from __future__ import annotations

from modules.labsheet_tools.tools._protocols import (
    PROTOCOLS,
    reagent_block,
    protocol_source_record,
)
from modules.labsheet_tools.tools.colony_calculator import colony_calculator
from modules.labsheet_tools.tools.verify_edit import verify_edit
from modules.crispr_tools.tools.citations import cites, format_citations

# predict_editing_efficiency is imported lazily inside run() — keeps
# lab_sheet usable in lightweight envs that don't ship the crispr scoring
# tool, and lets us fail soft (skip the prediction, render the sheet).
_PREDICT_EFFICIENCY_KNOWN = {"cas9", "cas12a"}

_NO_RESCUE_ANTIBIOTICS = {
    # Beta-lactams: act on cell-wall synthesis, slow enough that
    # resistance can be expressed during plating.
    "amp", "ampicillin",
    "carb", "carbenicillin",
    # Yeast auxotrophic markers: not antibiotics at all, plated on
    # dropout media (e.g. -URA, -LEU). No recovery step needed.
    "ura3", "leu2", "his3", "trp1", "lys2", "met15",
}


def _needs_rescue(antibiotic: str) -> bool:
    """Rescue is required unless the antibiotic is a beta-lactam (Amp/Carb)
    or a yeast auxotrophic marker (URA3/LEU2/HIS3/TRP1/...).

    Substring check so compound strings like
    'URA3 (yeast) / Ampicillin (bacteria)' still suppress rescue when
    they contain a no-rescue marker."""
    lowered = (antibiotic or "").strip().lower()
    if not lowered:
        return False
    return not any(marker in lowered for marker in _NO_RESCUE_ANTIBIOTICS)

# Map of CRISPRDelivery method name → registry protocol key. Lets the
# construction record specify {"step_type": "CRISPRDelivery",
# "parameters": {"method": "electroporation"}} and have the lab sheet
# render the appropriate canonical recipe.
_CRISPR_DELIVERY_PROTOCOLS = {
    "rnp_assembly": "crispr_rnp_assembly",
    "electroporation": "crispr_electroporation",
    "lipofection": "crispr_lipofection",
    "ivt_sgrna": "crispr_ivt_sgrna",
    "ivt": "crispr_ivt_sgrna",
}


def _protocols_io_search_url(query: str) -> str:
    """Build a deep link to a protocols.io keyword search. Safe URL-encoding
    via the stdlib so users can click through and see real, peer-shared
    protocols beyond the curated registry."""
    from urllib.parse import quote_plus
    return f"https://www.protocols.io/search?q={quote_plus(query)}"

_COMPLEMENT = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def _revcomp(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def _build_parts_map(parts: list[dict]) -> dict[str, dict]:
    """Index the parts list by name for sequence/template lookups."""
    out: dict[str, dict] = {}
    for p in parts or []:
        name = p.get("name")
        if name:
            out[name] = p
    return out


def _is_circular_part(part: dict) -> bool:
    """Plasmids/backbones are circular; everything else linear."""
    if not part:
        return False
    return part.get("part_type", "").lower() in {"plasmid", "circular_dna"}


def _locate_primer(primer_seq: str, template: str, circular: bool) -> int | None:
    """Find the position where the primer binds the template. Returns the
    0-based start of the whole primer string in `template`, or None.

    Handles both orientations:
    - Forward primer (5' overhang ... 3' anneal): anneal is at the END.
    - Reverse-complemented reverse primer (5' anneal ... 3' overhang): the
      anneal is at the START of the rev_rc string. We try the 3' end first
      (forward case) and fall back to the 5' end (rev_rc case)."""
    if not primer_seq or not template:
        return None
    primer_seq = primer_seq.upper()
    template = template.upper()
    search = template + (template if circular else "")
    # Pass 1: 3'-end anneal (forward primer with 5' overhang).
    for anneal_len in (len(primer_seq), 25, 22, 20, 18, 15):
        if anneal_len > len(primer_seq):
            continue
        anneal = primer_seq[-anneal_len:]
        idx = search.find(anneal)
        if idx != -1:
            return idx - (len(primer_seq) - anneal_len)
    # Pass 2: 5'-end anneal (rev_rc — the rev primer's 3' anneal becomes
    # the 5' end of the rev_rc string). Returns the 0-based start of the
    # rev_rc string on the template.
    for anneal_len in (len(primer_seq), 25, 22, 20, 18, 15):
        if anneal_len > len(primer_seq):
            continue
        anneal = primer_seq[:anneal_len]
        idx = search.find(anneal)
        if idx != -1:
            return idx
    return None


def _expected_amplicon_size(
    fwd_seq: str, rev_seq: str, template_seq: str, circular: bool
) -> int | None:
    """Compute the length of the PCR product given primer + template
    sequences. Returns None if either primer fails to locate."""
    if not (fwd_seq and rev_seq and template_seq):
        return None
    fwd_pos = _locate_primer(fwd_seq, template_seq, circular)
    rev_rc = _revcomp(rev_seq)
    rev_pos = _locate_primer(rev_rc, template_seq, circular)
    if fwd_pos is None or rev_pos is None:
        return None
    rev_end = rev_pos + len(rev_seq)
    n = len(template_seq)
    if circular:
        # All search indices live in [0, 2n). Normalize fwd to [0, n),
        # then advance rev_end forward until it sits past fwd_pos.
        fwd_pos %= n
        while rev_end <= fwd_pos:
            rev_end += n
        size = rev_end - fwd_pos
        # Vector PCR primers in Gibson assembly point OUTWARD from the
        # insertion site, so the small "gap" between them is the overlap
        # region — the actual amplicon wraps the whole plasmid the other
        # way. If the apparent size is shorter than the two primers
        # combined, interpret as wrap-around.
        if size < len(fwd_seq) + len(rev_seq):
            size = n - size + len(fwd_seq) + len(rev_seq)
        return size
    if rev_end <= fwd_pos:
        return None
    return rev_end - fwd_pos


def _gel_size_for_pcr(op: dict, parts_map: dict[str, dict]) -> int | None:
    """Compute expected gel-band size in bp for a PCR op, or None if
    sequences are unavailable."""
    p = op.get("parameters", {})
    fwd_name = p.get("forward_primer", "")
    rev_name = p.get("reverse_primer", "")
    tmpl_name = p.get("template", "")
    fwd_part = parts_map.get(fwd_name)
    rev_part = parts_map.get(rev_name)
    tmpl_part = parts_map.get(tmpl_name)
    if not (fwd_part and rev_part and tmpl_part):
        return None
    return _expected_amplicon_size(
        fwd_part.get("sequence", ""),
        rev_part.get("sequence", ""),
        tmpl_part.get("sequence", ""),
        circular=_is_circular_part(tmpl_part),
    )

# Default success rate when the user does not supply one. 0.8 reflects a
# typical 2-3 fragment Golden Gate / Gibson assembly into a competent E. coli
# strain like Mach1 or DH5a (NEB E1601 / E2611 product literature; Engler 2008).
# For CRISPR editing experiments, the caller should pass a real
# editing_efficiency or colony_preset so the screening burden is computed
# from a published benchmark instead.
_DEFAULT_CLONING_SUCCESS = 0.8

# Heuristic mapping from common cell/strain names to a colony_calculator
# preset. Keeps the lab sheet sensible when the user runs CRISPR workflows
# without supplying editing_efficiency by hand.
_MAMMALIAN_HINTS = ("HEK", "K562", "HELA", "U2OS", "IPSC", "JURKAT", "CHO", "293")
_YEAST_HINTS = ("SACCHAROMYCES", "YEAST", "CEREVISIAE", "S. CEREVISIAE", "S.CEREVISIAE", "BY4741", "BY4742", "W303")
_ECOLI_HINTS = ("E. COLI", "E.COLI", "ESCHERICHIA", "DH5A", "DH10B", "MACH1", "TOP10", "BL21", "JM109", "NEB10", "STBL3", "STBL4")


def _guess_preset(cells: str) -> str | None:
    """Map a strain/cell-line name to a colony_calculator preset, or None
    to fall back to _DEFAULT_CLONING_SUCCESS.

    Mammalian → cas9_plasmid_mammalian (20%)
    Yeast     → cas9_yeast (50%)
    E. coli   → None (the default 80% cloning rate is right for routine
                vector cloning; CRISPR editing in E. coli is uncommon and
                would need an explicit caller-supplied preset anyway.)
    """
    if not cells:
        return None
    up = cells.upper()
    for hint in _MAMMALIAN_HINTS:
        if hint in up:
            return "cas9_plasmid_mammalian"
    for hint in _YEAST_HINTS:
        if hint in up:
            return "cas9_yeast"
    return None  # E. coli / unknown → use cloning success rate default


def _colony_labels(n: int) -> list[str]:
    """Generate n short labels: A..Z, then AA..AZ, BA..ZZ. Used for tube
    labels like A1A, A1B, ..."""
    labels = []
    for i in range(n):
        if i < 26:
            labels.append(chr(65 + i))
        else:
            labels.append(chr(65 + (i // 26) - 1) + chr(65 + (i % 26)))
    return labels


def _build_colony_plan(
    op: dict,
    editing_efficiency: float | None,
    colony_preset: str | None,
    desired_clones: int,
    confidence: float,
    is_crispr: bool = False,
) -> dict:
    """Call colony_calculator to decide how many colonies to pick for a
    Transform step. Returns a dict with n, labels, the efficiency that
    was used, and the upstream citations."""
    params = op.get("parameters", {})
    cells = params.get("cells", "")

    # precedence: explicit editing_efficiency > explicit preset > strain hint
    # > default cloning success rate.
    if editing_efficiency is not None:
        result = colony_calculator(
            editing_efficiency=editing_efficiency,
            desired_clones=desired_clones,
            confidence=confidence,
        )
        used_preset = None
    elif colony_preset is not None:
        result = colony_calculator(
            preset=colony_preset,
            desired_clones=desired_clones,
            confidence=confidence,
        )
        used_preset = colony_preset
    else:
        guess = _guess_preset(cells)
        if guess is not None:
            result = colony_calculator(
                preset=guess,
                desired_clones=desired_clones,
                confidence=confidence,
            )
            used_preset = guess
        elif is_crispr:
            # E. coli CRISPR editing: use the published benchmark preset
            # rather than the generic cloning success rate.
            result = colony_calculator(
                preset="cas9_ecoli",
                desired_clones=desired_clones,
                confidence=confidence,
            )
            used_preset = "cas9_ecoli"
        else:
            result = colony_calculator(
                editing_efficiency=_DEFAULT_CLONING_SUCCESS,
                desired_clones=desired_clones,
                confidence=confidence,
            )
            used_preset = None

    n = result["colonies_to_pick"]
    return {
        "n": n,
        "labels": _colony_labels(n),
        "efficiency": result["editing_efficiency"],
        "preset": used_preset,
        "recommendation": result["recommendation"],
        "citations": result.get("citations", []),
    }

# Disclaimer included in every lab sheet output. The lab sheet is a
# convenience starting point — the user is expected to verify volumes
# against current manufacturer documentation before running the protocol.
_VERIFY_DISCLAIMER = (
    "Volumes and conditions in this lab sheet are transcribed from canonical "
    "sources (NEB, Zymo, original method papers) into the project's protocol "
    "registry. They reflect those sources as of the registry's last update — "
    "ALWAYS confirm against current manufacturer documentation and your lab's "
    "SOPs before benchwork. See the 'protocol_sources' field for citations."
)

_ENZYME_NOTE = (
    "note:\n"
    "Never let enzymes warm up! Only take the enzyme cooler out of the freezer when you\n"
    "are actively using it, and only take the tubes out of it when actively dispensing."
)

_TSV_HEADERS = [
    "section", "label",
    "primer1", "primer1_location",
    "primer2", "primer2_location",
    "template", "template_location",
    "fragments", "product", "product_location",
    "strain", "antibiotic", "incubation",
    "destination", "program", "notes",
]


def _col_table(headers: list[str], rows: list[list[str]]) -> str:
    widths = [len(h) for h in headers]
    for row in rows:
        for i, cell in enumerate(row):
            if i < len(widths):
                widths[i] = max(widths[i], len(str(cell)))
    lines = []
    lines.append("  ".join(h.ljust(widths[i]) for i, h in enumerate(headers)).rstrip())
    for row in rows:
        lines.append("  ".join(str(cell).ljust(widths[i]) for i, cell in enumerate(row) if i < len(widths)).rstrip())
    return "\n".join(lines)


def _tsv_row(values: dict) -> str:
    return "\t".join(str(values.get(h, "")) for h in _TSV_HEADERS)


def _format_pcr_section(ops: list[dict], thread: str) -> str:
    rows = []
    source_map: dict[str, str] = {}

    fwd_set: set[str] = set()
    rev_set: set[str] = set()
    tmpl_set: set[str] = set()

    for i, op in enumerate(ops, start=1):
        label = f"{thread}{i}"
        p = op.get("parameters", {})
        fwd = p.get("forward_primer", "")
        rev = p.get("reverse_primer", "")
        tmpl = p.get("template", "")
        product = op.get("output", "")
        rows.append([label, fwd, rev, tmpl, product])
        if fwd: fwd_set.add(fwd)
        if rev: rev_set.add(rev)
        if tmpl: tmpl_set.add(tmpl)

        oligo_col = 65 + (i - 1) * 2
        if fwd:
            source_map.setdefault(fwd, f"oligos1/{chr(oligo_col)}1")
        if rev:
            source_map.setdefault(rev, f"oligos1/{chr(oligo_col + 1)}1")
        if tmpl:
            source_map.setdefault(tmpl, f"templates/{thread}{i}")

    samples_table = _col_table(["label", "primer1", "primer2", "template", "product"], rows)
    src_rows = [[name, loc, ""] for name, loc in source_map.items()]
    source_table = _col_table(["label", "location", "note"], src_rows)

    # Reaction recipe / mastermix block.
    # Per-rxn (25 uL) volumes match standard Q5/Phusion setups:
    #   16.75 uL ddH2O, 5 uL 5x buffer, 2.5 uL 2mM dNTPs,
    #   0.5 uL each primer, 0.5 uL template, 0.25 uL polymerase.
    # For 4+ samples, build a mastermix that pools every component
    # shared across all reactions (10% dead volume, per Anderson lab).
    n = len(ops)
    fwd_shared = len(fwd_set) == 1
    rev_shared = len(rev_set) == 1
    tmpl_shared = len(tmpl_set) == 1
    fwd_name = next(iter(fwd_set)) if fwd_shared else None
    rev_name = next(iter(rev_set)) if rev_shared else None
    tmpl_name = next(iter(tmpl_set)) if tmpl_shared else None

    recipe_lines: list[str] = []
    if n >= 4:
        scale = n * 1.1  # 10% dead volume
        mm: list[tuple[str, float]] = [
            ("ddH2O", 16.75),
            ("5x Q5 buffer", 5.0),
            ("2 mM dNTPs", 2.5),
        ]
        per_rxn: list[tuple[str, float]] = []
        if fwd_shared:
            mm.append((f"{fwd_name} (10 uM)", 0.5))
        else:
            per_rxn.append(("forward primer (10 uM)", 0.5))
        if rev_shared:
            mm.append((f"{rev_name} (10 uM)", 0.5))
        else:
            per_rxn.append(("reverse primer (10 uM)", 0.5))
        if tmpl_shared:
            mm.append((f"{tmpl_name} (template)", 0.5))
        else:
            per_rxn.append(("template", 0.5))
        mm.append(("Q5 polymerase", 0.25))
        mm_total = sum(v for _, v in mm)
        per_rxn_mm_volume = mm_total
        recipe_lines.append(f"mastermix (covers {n} reactions + 10% dead volume):")
        for name, v in mm:
            recipe_lines.append(f"  {round(v * scale, 3)} uL  {name}")
        recipe_lines.append("")
        recipe_lines.append("reaction (per tube):")
        recipe_lines.append(f"  {per_rxn_mm_volume} uL mastermix")
        for name, v in per_rxn:
            recipe_lines.append(f"  {v} uL {name}")
    else:
        recipe_lines.append("reaction (per tube, 25 uL):")
        recipe_lines.append("  16.75 uL ddH2O")
        recipe_lines.append("  5 uL 5x Q5 buffer")
        recipe_lines.append("  2.5 uL 2 mM dNTPs")
        recipe_lines.append("  0.5 uL forward primer (10 uM)")
        recipe_lines.append("  0.5 uL reverse primer (10 uM)")
        recipe_lines.append("  0.5 uL template")
        recipe_lines.append("  0.25 uL Q5 polymerase")

    return "\n".join([
        f"{thread}: PCR",
        "",
        *recipe_lines,
        "",
        "samples:",
        samples_table,
        "",
        "source:",
        source_table,
        "",
        "destination: thermocycler1A",
        "program: Q5/Q5-4K",
        "",
        _ENZYME_NOTE,
    ])


def _format_gel_dpni_section(
    ops: list[dict], thread: str, parts_map: dict[str, dict] | None = None
) -> str:
    parts_map = parts_map or {}
    rows = []
    for i, op in enumerate(ops, start=1):
        label = f"{thread}{i}"
        product = op.get("output", "")
        size = _gel_size_for_pcr(op, parts_map)
        size_str = f"{size} bp" if size else "?"
        rows.append([label, size_str, product])

    samples_table = _col_table(["reaction", "size", "product"], rows)

    return "\n".join([
        f"{thread}: Gel and DpnI",
        "",
        "source: thermocycler1A",
        "",
        "samples:",
        samples_table,
        "",
        "protocol:",
        'In new PCR tubes, combine, mix and quick spin 6 uL of "1x Load" and 2 uL of PCR product.',
        "Run a gel with the full volume of the mixture in each well, and an additional well with",
        "5 uL of BstEII ladder.",
        "Add 0.5 uL DpnI to each PCR reaction, mix, quick spin, run thermocycler.",
        "",
        "destination: thermocycler1A",
        "program: main/SPE1",
        "",
        _ENZYME_NOTE,
    ])


def _format_zymo_section(ops: list[dict], thread: str) -> str:
    rows = []
    for i, op in enumerate(ops, start=1):
        label = f"{thread}{i}"
        label_p = f"{thread}{i}p"
        product = op.get("output", "")
        dest = f"box{thread}/{chr(64 + i)}1"
        rows.append([label, label_p, "50 uL", dest, f"{product}/pcrpdt"])

    samples_table = _col_table(["reaction", "label", "elution_volume", "destination", "product"], rows)

    return "\n".join([
        f"{thread}: Zymo",
        "",
        "source: thermocycler1A",
        "",
        "samples:",
        samples_table,
    ])


def _format_assemble_section(op: dict, thread: str) -> str:
    strategy = op.get("step_type", "GoldenGate")
    inputs = op.get("inputs", [])
    output = op.get("output", "")
    params = op.get("parameters", {})

    dna_mix = "\n".join(f"5 uL {inp}" for inp in inputs)

    # Recipe volumes come from the protocol registry, not hard-coded here.
    # See _protocols.PROTOCOLS for sources (Engler 2008, Gibson 2009, NEB).
    if strategy == "GoldenGate":
        proto_key = "goldengate_bsai"
        reaction = "reaction:\n" + reagent_block(proto_key)
        destination = "thermocycler1A"
        program = "main/GG1"
    else:
        proto_key = "gibson_neb_2x"
        reaction = "reaction:\n" + reagent_block(proto_key)
        destination = "thermocycler1B"
        program = "main/GIB2"

    src_rows = [[inp, f"box{thread}/{chr(65 + i)}1"] for i, inp in enumerate(inputs)]
    source_table = _col_table(["dna", "location"], src_rows)

    frags = ",".join(inputs)
    samples_table = _col_table(["label", "fragments", "product"], [[thread, frags, output]])

    return "\n".join([
        f"{thread}: Assemble",
        "",
        "DNA Mix:",
        dna_mix,
        "",
        reaction,
        "",
        "source:",
        source_table,
        "",
        "samples:",
        samples_table,
        "",
        f"destination: {destination}",
        f"program: {program}",
        "",
        _ENZYME_NOTE,
    ])


def _split_restriction_ligation_reagents() -> tuple[str, str]:
    """Split the `restriction_ligation` protocol's reagents into two
    blocks: one combined Digest recipe (vector + insert) and one Ligate
    recipe. The protocol entry uses '--- ... ---' marker rows to group
    them; this helper walks those groups."""
    proto = PROTOCOLS["restriction_ligation"]
    digest_lines: list[str] = []
    ligate_lines: list[str] = []
    bucket: list[str] | None = None
    for reagent, vol in proto.reagents:
        if reagent.startswith("---"):
            label = reagent.strip("- ").lower()
            if "digest" in label:
                bucket = digest_lines
                digest_lines.append(reagent)
            elif "ligation" in label or "ligate" in label:
                bucket = ligate_lines
            continue
        if bucket is not None:
            bucket.append(f"{vol} {reagent}")
    return "\n".join(digest_lines), "\n".join(ligate_lines)


def _format_restriction_ligation_section(op: dict, thread: str) -> str:
    """Render a RestrictionLigation step as TWO labsheet sub-sections —
    `{thread}: Digest` and `{thread}: Ligate` — matching the per-operation
    labpacket model from the LabPlanner lecture (one labsheet per physical
    bench session)."""
    proto = PROTOCOLS["restriction_ligation"]
    params = op.get("parameters", {})
    inputs = op.get("inputs", [])
    output = op.get("output", "")
    enzyme = params.get("enzyme", "EcoRI-HF + BamHI-HF")

    backbone = inputs[0] if len(inputs) > 0 else "vector"
    insert = inputs[1] if len(inputs) > 1 else "insert"

    digest_recipe, ligate_recipe = _split_restriction_ligation_reagents()

    digest_src_rows = [
        [backbone, f"box{thread}/A1"],
        [insert, f"box{thread}/B1"],
    ]
    digest_source_table = _col_table(["dna", "location"], digest_src_rows)
    digest_samples_table = _col_table(
        ["label", "dna", "source", "product"],
        [
            [f"{thread}1", backbone, f"box{thread}/A1", f"{backbone}/dig"],
            [f"{thread}2", insert,   f"box{thread}/B1", f"{insert}/dig"],
        ],
    )

    digest_section = "\n".join([
        f"{thread}: Digest",
        "",
        f"enzyme(s): {enzyme}",
        "",
        "reaction:",
        digest_recipe,
        "",
        "samples:",
        digest_samples_table,
        "",
        "source:",
        digest_source_table,
        "",
        "destination: thermocycler1A",
        "program: 37C 1 hr; 65C 20 min heat-inactivate; gel-purify",
        "",
        _ENZYME_NOTE,
    ])

    ligate_samples_table = _col_table(
        ["label", "digest", "source", "product"],
        [[f"{thread}1", f"{backbone}/dig + {insert}/dig", f"box{thread}/A1+B1", output]],
    )

    ligate_section = "\n".join([
        f"{thread}: Ligate",
        "",
        "reaction:",
        ligate_recipe,
        "",
        "samples:",
        ligate_samples_table,
        "",
        "destination: thermocycler1A",
        "program: 16C overnight (or RT 1 hr for sticky ends)",
        "",
        f"notes:\n{proto.notes}",
        "",
        _ENZYME_NOTE,
    ])

    return digest_section + "\n\n" + ligate_section


def _format_typeiis_oligo_section(op: dict, thread: str) -> str:
    """Render a TypeIISOligoCloning step. The construction record's
    inputs are [top_oligo, bottom_oligo, vector] and parameters carry
    the Type IIS enzyme name (BbsI / BsmBI). Without this section, the
    workflow is biologically incomplete — the lab sheet would jump
    straight to Transform with no annealed-and-ligated insert in hand.
    """
    proto_key = "typeiis_oligo_cloning"
    proto = PROTOCOLS[proto_key]
    params = op.get("parameters", {})
    inputs = op.get("inputs", [])
    output = op.get("output", "")
    enzyme = params.get("enzyme", "BbsI-HF")
    overhangs = params.get("overhangs", "vector-specific (e.g. CACC/AAAC)")

    top = inputs[0] if len(inputs) > 0 else "top_oligo"
    bot = inputs[1] if len(inputs) > 1 else "bottom_oligo"
    vec = inputs[2] if len(inputs) > 2 else "vector"

    src_rows = [
        [top, "oligos1/A1"],
        [bot, "oligos1/B1"],
        [vec, f"box{thread}/A1"],
    ]
    source_table = _col_table(["dna", "location"], src_rows)
    samples_table = _col_table(
        ["label", "fragments", "product"],
        [[thread, f"{top}+{bot} duplex / {vec}", output]],
    )

    return "\n".join([
        f"{thread}: TypeIISOligoCloning",
        "",
        f"protocol: {proto.name}",
        f"enzyme: {enzyme}",
        f"overhangs: {overhangs}",
        "",
        "reaction:",
        reagent_block(proto_key),
        "",
        f"program: {proto.program}",
        "",
        "source:",
        source_table,
        "",
        "samples:",
        samples_table,
        "",
        "destination: thermocycler1A",
        "",
        f"notes:\n{proto.notes}",
        "",
        _ENZYME_NOTE,
    ])


def _format_transform_section(op: dict, thread: str) -> str:
    params = op.get("parameters", {})
    inputs = op.get("inputs", [])

    strain = params.get("cells", "Mach1")
    antibiotic = params.get("selection", "Amp")
    temp = params.get("temperature_c", 37)
    rescue = _needs_rescue(antibiotic)
    rescue_str = "yes" if rescue else "no"

    product = inputs[0] if inputs else op.get("output", "")
    rows = [[f"{thread}1", product, strain, antibiotic, f"{temp}°C"]]
    samples_table = _col_table(["label", "product", "strain", "antibiotic", "incubate"], rows)

    lines = [
        f"{thread}: Transform",
        "",
        "source: thermocycler1A",
        "",
        "samples:",
        samples_table,
        "",
        f"rescue_required: {rescue_str}",
    ]

    if rescue:
        lines += [
            "",
            "note:",
            f"Because the antibiotic is {antibiotic}, there needs to be a rescue step. After the",
            "transformation completes, transfer the cells to a 1.5 mL eppendorf tube (labeled on lid",
            "with label) and 200 uL of 2YT medium. Tape the tubes to the inside of the shaker.",
            "Another experimentalist will perform plating.",
        ]

    return "\n".join(lines)


def _format_plate_section(op: dict, thread: str) -> str:
    """Render the post-transformation plating step. Same per-sample fields
    as Transform (label/product/strain/antibiotic/incubate) but framed as a
    distinct labsheet section so the per-step labpacket model matches the
    professor's LabPlanner format."""
    params = op.get("parameters", {})
    inputs = op.get("inputs", [])

    strain = params.get("cells", "Mach1")
    antibiotic = params.get("selection", "Amp")
    temp = params.get("temperature_c", 37)
    product = inputs[0] if inputs else op.get("output", "")

    rows = [[f"{thread}1", product, strain, antibiotic, f"{temp}°C"]]
    samples_table = _col_table(["label", "product", "strain", "antibiotic", "incubate"], rows)

    return "\n".join([
        f"{thread}: Plate",
        "",
        "samples:",
        samples_table,
        "",
        "protocol:",
        f"Warm up (or prepare) a plate containing {antibiotic}.",
        "Label the plate with date, product name, lab member name.",
        "Pipette 100 uL of transformed cells onto plate.",
        "Spread the cells with glass beads or a sterile wand.",
        f"Incubate the plates upside down, without shaking, in a {temp}°C incubator.",
    ])


def _format_pick_section(op: dict, thread: str, plan: dict) -> str:
    params = op.get("parameters", {})
    inputs = op.get("inputs", [])
    strain = params.get("cells", "Mach1")
    antibiotic = params.get("selection", "Amp")
    product = inputs[0] if inputs else op.get("output", "")

    n = plan["n"]
    labels = plan["labels"]
    label_str = ", ".join(f"{thread}1{ch}" for ch in labels)
    rows = [
        [f"{thread}1", product, strain, antibiotic, "37°C", str(n), label_str]
    ]
    samples_table = _col_table(
        ["source", "product", "strain", "antibiotic", "incubate", "number", "labels"], rows
    )

    rationale_line = (
        f"colonies to pick: {n} "
        f"(efficiency {plan['efficiency']:.0%}"
        + (f", preset={plan['preset']}" if plan["preset"] else "")
        + ")"
    )

    return "\n".join([
        f"{thread}: Pick",
        "",
        "samples:",
        samples_table,
        "",
        rationale_line,
        "",
        "rationale:",
        plan["recommendation"],
    ])


def _format_miniprep_section(op: dict, thread: str, plan: dict) -> str:
    params = op.get("parameters", {})
    inputs = op.get("inputs", [])
    product = inputs[0] if inputs else op.get("output", "")

    rows = []
    for i, ch in enumerate(plan["labels"]):
        loc = f"box{thread}/{chr(67 + i)}1" if i < 24 else f"box{thread}/overflow{i}"
        rows.append([f"{thread}1{ch}", f"{product}-{ch}", loc])
    samples_table = _col_table(["culture", "label", "location"], rows)

    return "\n".join([
        f"{thread}: Miniprep",
        "",
        "samples:",
        samples_table,
        "",
        "note:",
        f"Write the short label (e.g. {thread}1{plan['labels'][0]}) on the top of the eppendorf,",
        "put the full name (the label column) on the side of the tube.",
    ])


def _normalize_seq_primers(seq_primers: list[dict] | None) -> list[dict]:
    """Default to L4440 if no guide-specific primers are provided. When
    the caller passes verify_edit's output, those primers replace L4440
    so the lab sheet sequences across the actual edit site instead of a
    generic plasmid backbone landmark."""
    if not seq_primers:
        return [{"name": "L4440", "location": "oligos1/J1",
                 "note": "generic plasmid sequencing primer"}]
    out = []
    for i, p in enumerate(seq_primers):
        out.append({
            "name": p.get("name", f"verify_primer_{i+1}"),
            "location": p.get("location", f"oligos1/J{i+1}"),
            "sequence": p.get("sequence", ""),
            "note": p.get("note", "guide-specific verification primer (verify_edit)"),
        })
    return out


def _format_sequencing_section(
    op: dict,
    thread: str,
    plan: dict,
    seq_primers: list[dict] | None = None,
    protospacer: str | None = None,
) -> str:
    params = op.get("parameters", {})
    inputs = op.get("inputs", [])
    product = inputs[0] if inputs else op.get("output", "")
    seq_primers = _normalize_seq_primers(seq_primers)

    src_rows = []
    for i, ch in enumerate(plan["labels"]):
        loc = f"box{thread}/{chr(67 + i)}1" if i < 24 else f"box{thread}/overflow{i}"
        src_rows.append([f"{thread}1{ch}", loc, f"{product}-{ch}"])
    for sp in seq_primers:
        src_rows.append([sp["name"], sp["location"], sp["note"]])
    source_table = _col_table(["label", "location", "product"], src_rows)

    example_label = f"{thread}1{plan['labels'][0]}"
    instructions: list[str] = []
    for sp in seq_primers:
        save_loc = sp.get("save_location", sp["location"])
        instructions += [
            f"- Resuspend the oligo {sp['name']} to the appropriate volume for a 100 uM stock",
            f"- In an eppendorf, prepare a 2.66 uM dilute stock of {sp['name']} as:",
            "    o 487 uL ddH2O",
            "    o 13.3 uL of 100 uM oligo",
        ]
    instructions += [
        "- For each plasmid listed, mix the following sequencing reactions in an eppendorf tube:",
        "    o 6 uL ddH2O",
        "    o 4 uL miniprep DNA (undiluted)",
        "    o 3 uL oligo (2.66 uM)",
        f'- Label the tops of the tubes with the "label", ie "{example_label}"',
    ]
    for sp in seq_primers:
        save_loc = sp.get("save_location", sp["location"])
        instructions.append(f"- When done, save the {sp['name']} stock at: {save_loc}")
    instructions.append(
        "- Take the sequencing reactions and order form to: 237 Stanley Hall (second floor cold room)"
    )

    expected_block: list[str] = []
    if protospacer:
        expected_block = [
            "",
            "expected read should contain:",
            f"  {protospacer}  (20-nt protospacer; confirms guide cloning)",
            "  + sgRNA scaffold + terminator downstream",
        ]

    trailing_note = (
        []
        if any(sp["name"] == "L4440" for sp in seq_primers)
        else [
            "",
            "Note: using guide-specific verification primers (verify_edit) — read WILL span the edit site.",
        ]
    )

    return "\n".join([
        f"{thread}: Sequencing (plasmid — confirms vector clone, NOT the genomic edit)",
        "",
        "sources:",
        source_table,
        "",
        "Instructions:",
        *instructions,
        *expected_block,
        *trailing_note,
    ])


def _format_crispr_delivery_section(op: dict, thread: str) -> str:
    """Render a CRISPRDelivery step. The op's `parameters.method` selects
    which protocol from the registry to render."""
    params = op.get("parameters", {})
    method = (params.get("method") or "rnp_assembly").lower().strip()
    proto_key = _CRISPR_DELIVERY_PROTOCOLS.get(method)
    if proto_key is None:
        raise ValueError(
            f"Unknown CRISPRDelivery method '{method}'. "
            f"Choose one of: {sorted(_CRISPR_DELIVERY_PROTOCOLS)}."
        )
    proto = PROTOCOLS[proto_key]
    inputs = op.get("inputs", [])
    output = op.get("output", "")

    return "\n".join([
        f"{thread}: CRISPRDelivery ({method})",
        "",
        f"protocol: {proto.name}",
        "",
        "reaction:",
        reagent_block(proto_key),
        "",
        f"program: {proto.program}",
        "",
        "inputs: " + (", ".join(inputs) if inputs else "(none)"),
        f"product: {output}",
        ("\nnotes:\n" + proto.notes) if proto.notes else "",
        _ENZYME_NOTE,
    ]).rstrip()


def _format_edit_verification_placeholder(
    thread: str, protospacer: str, nuclease: str
) -> str:
    """Always-emit placeholder for the post-edit GENOMIC verification
    step when the design phase identified a protospacer but the caller
    didn't supply a genomic reference. Points the user at crispr_verify_edit
    so they can design specific primers when they're ready to verify."""
    return "\n".join([
        f"{thread}: EditVerification (post-edit genomic Sanger) — placeholder",
        "",
        "After your CRISPR edit succeeds, verify the genomic edit:",
        "  1. Extract genomic DNA from edited cells.",
        "  2. PCR-amplify the genomic locus around the cut site (~300-500 bp window).",
        "  3. Sanger-sequence the PCR product with one of the flanking primers.",
        "  4. Upload the .ab1 trace + amplicon sequence to Synthego ICE or TIDE.",
        "  5. Run crispr_interpret_ice_tide on the editing % + R-squared from ICE/TIDE.",
        "",
        f"To get specific verification primers for this guide ({protospacer}, {nuclease}),",
        "run crispr_verify_edit with the genomic reference for your locus. Without a",
        "genomic reference, this section is a placeholder — primer design needs the",
        "actual sequence around the edit site.",
    ])


def _format_edit_verification_section(
    thread: str,
    ve: dict,
    nuclease: str,
    delivery: str,
) -> str:
    """Render the post-edit GENOMIC verification protocol from a
    verify_edit result. This is distinct from the plasmid Sequencing
    step — those Sanger reads confirm the vector was cloned correctly,
    while this section confirms the actual genomic edit (e.g. mouse
    Tyr locus) AFTER transfection / nucleofection / electroporation.

    Workflow rendered:
      1. Extract genomic DNA from edited cells (post-recovery).
      2. PCR-amplify the locus using verify_edit's flanking primers.
      3. Run a confirmation gel for amplicon size.
      4. Sanger-sequence the amplicon (3 uL of one of the same primers).
      5. Upload the .ab1 trace + amplicon sequence + cut offset to
         Synthego ICE or run TIDE locally.
    """
    fwd = ve.get("forward_primer", "")
    rev = ve.get("reverse_primer", "")
    fwd_tm = ve.get("forward_primer_tm", "?")
    rev_tm = ve.get("reverse_primer_tm", "?")
    cut_pos = ve.get("cut_position", "?")
    amp_len = ve.get("amplicon_length", "?")
    cut_off = ve.get("cut_offset_in_amplicon", "?")
    primer_warnings = ve.get("primer_warnings", []) or []

    primer_table = _col_table(
        ["primer", "sequence", "Tm (C)", "location"],
        [
            ["verify_F (genomic)", fwd or "(verify_edit could not design)", str(fwd_tm), "oligos1/M1"],
            ["verify_R (genomic)", rev or "(verify_edit could not design)", str(rev_tm), "oligos1/M2"],
        ],
    )

    gdna_step = (
        "Day 0 (post-edit): harvest edited cells (or genomic-DNA template) "
        f"~{2 if delivery in ('rnp','electroporation','lipofection') else 3} "
        "days after delivery to allow time for indel formation."
    )
    if delivery in ("rnp", "electroporation"):
        gdna_step += " For RNP/electroporation, 48 hr is typical."
    elif delivery in ("plasmid", "lentivirus", "aav"):
        gdna_step += (
            " For plasmid transfection/lentivirus/AAV, wait 72 hr or longer "
            "for full Cas9 expression and editing."
        )
    gdna_step += (
        " Extract genomic DNA via column kit (e.g. Qiagen DNeasy, Zymo "
        "Quick-DNA Miniprep), elute in 50 uL ddH2O, normalize to ~50 ng/uL."
    )

    ice_url = "https://ice.synthego.com/"
    interpretation = ve.get("interpretation_guide", "")

    pcr_step = (
        f"PCR (50 uL Q5 reaction): 25 uL Q5 2X master mix + 1.25 uL each "
        f"primer (10 uM stock) + 1 uL gDNA (~50 ng) + 21.5 uL ddH2O. "
        f"Program: 98C 30s; [98C 10s, 65C 30s, 72C 30s/kb] x 35; 72C 2 min; "
        f"4C hold. Expected amplicon: ~{amp_len} bp."
    )
    gel_step = (
        f"Run 5 uL of PCR product on 1.5% agarose to confirm a single "
        f"~{amp_len} bp band before Sanger submission."
    )
    sanger_step = (
        "Sanger premix (1 reaction per primer, 10 uL total): 4 uL PCR "
        "product (100-300 ng) + 3 uL primer (2.66 uM) + 3 uL ddH2O. "
        "Submit to your sequencing facility."
    )
    ice_step = (
        f"Once .ab1 traces are back: upload the EDITED trace + the UNEDITED "
        f"reference trace (or paste the amplicon sequence) to {ice_url} along "
        f"with the protospacer ({ve.get('protospacer','?')}). ICE reports KO "
        f"score + indel distribution. The cut should fall at offset {cut_off} "
        f"in the amplicon (position {cut_pos} in the reference). Run TIDE "
        f"(Brinkman 2014) for a second-opinion deconvolution."
    )

    warning_block = ""
    if primer_warnings:
        warning_block = "\nprimer_warnings:\n" + "\n".join(f"- {w}" for w in primer_warnings)

    handoff = (
        "After running ICE/TIDE, paste your result back. Required: "
        "the editing percentage (KO score for ICE, total indel % for TIDE) "
        "AND the R^2 fit value (both are reported on the same ICE/TIDE "
        "results page). Optional but useful: the indel size distribution "
        "(e.g. '+1: 45%, -3: 12%') so the dominant indel can be flagged, "
        "and which tool you used (ICE vs TIDE — default is ICE).\n"
        "Example prompt: 'My ICE score is 45% with R^2 = 0.93' or 'TIDE: "
        "62% indels, R^2 0.91, dominant indel +1 at 38%'.\n"
        "After interpret_ice_tide returns, the workflow continues with: "
        "(a) more single-clone Sanger if you need a homozygous edit "
        "(call interpret_ice_tide again per clone), "
        "(b) colony_calculator to size additional screening, or "
        "(c) functional validation of the edited line (out of scope here)."
    )

    return "\n".join([
        f"{thread}: EditVerification (post-edit genomic Sanger)",
        "",
        f"target nuclease: {nuclease}",
        f"protospacer: {ve.get('protospacer','?')}",
        f"cut position (reference): {cut_pos}",
        f"amplicon length: {amp_len} bp",
        f"cut offset within amplicon: {cut_off}",
        "",
        "primers (order separately from cloning oligos):",
        primer_table,
        warning_block,
        "",
        "step 1 — gDNA prep:",
        gdna_step,
        "",
        "step 2 — locus PCR:",
        pcr_step,
        "",
        "step 3 — confirmation gel:",
        gel_step,
        "",
        "step 4 — Sanger:",
        sanger_step,
        "",
        "step 5 — ICE/TIDE upload:",
        ice_step,
        "",
        "handoff:",
        handoff,
    ]).rstrip() + ("\n\n" + interpretation if interpretation else "")


def _format_unknown_step_section(op: dict, thread: str) -> str:
    """Last-resort renderer for a step_type lab_sheet doesn't recognize.
    Better than silently dropping it — the user sees the inputs/outputs
    and a flag that no canonical recipe is on file."""
    step_type = op.get("step_type", "Unknown")
    inputs = op.get("inputs", [])
    output = op.get("output", "")
    params = op.get("parameters", {})
    return "\n".join([
        f"{thread}: {step_type}",
        "",
        "WARNING: lab_sheet does not have a canonical recipe for this step.",
        "The construction file specified the following — fill in volumes and",
        "thermocycler settings against the manufacturer's protocol manually.",
        "",
        f"inputs: {', '.join(inputs) if inputs else '(none)'}",
        f"output: {output}",
        f"parameters: {params}",
    ])


def _maybe_predict_efficiency(
    protospacer: str | None,
    pam: str | None,
    nuclease: str,
    delivery: str,
    outcome: str,
) -> dict | None:
    """Best-effort call into predict_editing_efficiency. Returns the
    full prediction dict on success, or None if the inputs are
    insufficient or the tool errors out (we don't want lab_sheet to
    crash because the predictor was strict about a PAM)."""
    if not protospacer or not pam:
        return None
    if nuclease.lower() not in _PREDICT_EFFICIENCY_KNOWN:
        return None
    # Normalize delivery — predict_editing_efficiency only knows the
    # 5 categories below, but users (and the CRISPRDelivery step type)
    # may pass synonyms like "lipofection" or "nucleofection".
    _DELIVERY_ALIAS = {
        "lipofection": "plasmid",       # lipofected plasmid
        "transfection": "plasmid",
        "nucleofection": "electroporation",
        "rnp_assembly": "rnp",
    }
    delivery_norm = _DELIVERY_ALIAS.get(delivery.lower(), delivery.lower())
    try:
        from modules.crispr_tools.tools.predict_editing_efficiency import (
            predict_editing_efficiency,
        )
        return predict_editing_efficiency(
            protospacer=protospacer,
            pam=pam,
            nuclease=nuclease.lower(),
            delivery=delivery_norm,
            outcome=outcome,
        )
    except Exception:
        return None


def _format_direct_synthesis_section(op: dict, thread: str) -> str:
    inputs = op.get("inputs", [])
    output = op.get("output", "")
    construct = inputs[0] if inputs else output

    return "\n".join([
        f"{thread}: DirectSynthesis",
        "",
        "samples:",
        _col_table(["label", "construct"], [[f"{thread}1", output]]),
        "",
        "note:",
        f"Order {construct} as a synthetic gene fragment from a vendor (e.g. Twist, IDT).",
        "Resuspend per vendor instructions before use.",
    ])


def _build_tsv(
    operations: list[dict],
    thread: str,
    notes_str: str,
    transform_plans: list[dict] | None = None,
    parts_map: dict[str, dict] | None = None,
    sequencing_primers: list[dict] | None = None,
    verify_edit_result: dict | None = None,
    protospacer: str | None = None,
) -> str:
    transform_plans = transform_plans or []
    parts_map = parts_map or {}
    seq_primers = _normalize_seq_primers(sequencing_primers)
    rows = ["\t".join(_TSV_HEADERS)]

    pcr_ops = [op for op in operations if op.get("step_type") == "PCR"]
    assemble_ops = [op for op in operations if op.get("step_type") in ("GoldenGate", "Gibson")]
    transform_ops = [op for op in operations if op.get("step_type") == "Transform"]
    direct_ops = [op for op in operations if op.get("step_type") == "DirectSynthesis"]
    typeiis_ops = [op for op in operations if op.get("step_type") == "TypeIISOligoCloning"]
    rxn_lig_ops = [op for op in operations if op.get("step_type") == "RestrictionLigation"]

    source_map: dict[str, str] = {}

    # PCR rows
    for i, op in enumerate(pcr_ops, start=1):
        p = op.get("parameters", {})
        fwd = p.get("forward_primer", "")
        rev = p.get("reverse_primer", "")
        tmpl = p.get("template", "")
        product = op.get("output", "")
        oligo_col = 65 + (i - 1) * 2
        loc_fwd = f"oligos1/{chr(oligo_col)}1"
        loc_rev = f"oligos1/{chr(oligo_col + 1)}1"
        loc_tmpl = f"templates/{thread}{i}"
        source_map.setdefault(fwd, loc_fwd)
        source_map.setdefault(rev, loc_rev)
        source_map.setdefault(tmpl, loc_tmpl)
        rows.append(_tsv_row({
            "section": f"{thread}: PCR",
            "label": f"{thread}{i}",
            "primer1": fwd,
            "primer1_location": loc_fwd,
            "primer2": rev,
            "primer2_location": loc_rev,
            "template": tmpl,
            "template_location": loc_tmpl,
            "product": product,
            "destination": "thermocycler1A",
            "program": "Q5/Q5-4K",
            "notes": "Never let enzymes warm up!",
        }))

    # Gel/DpnI rows
    for i, op in enumerate(pcr_ops, start=1):
        size = _gel_size_for_pcr(op, parts_map)
        size_note = f"Expected band: {size} bp. " if size else ""
        rows.append(_tsv_row({
            "section": f"{thread}: Gel and DpnI",
            "label": f"{thread}{i}",
            "product": op.get("output", ""),
            "destination": "thermocycler1A",
            "program": "main/SPE1",
            "notes": size_note + '6 uL "1x Load" + 2 uL PCR product per well; BstEII ladder; 0.5 uL DpnI per reaction. Never let enzymes warm up!',
        }))

    # Zymo rows
    for i, op in enumerate(pcr_ops, start=1):
        product = op.get("output", "")
        rows.append(_tsv_row({
            "section": f"{thread}: Zymo",
            "label": f"{thread}{i}",
            "product": f"{product}/pcrpdt",
            "product_location": f"box{thread}/{chr(64 + i)}1",
            "notes": "Elute in 50 uL",
        }))

    # Assemble rows
    for op in assemble_ops:
        strategy = op.get("step_type", "GoldenGate")
        inputs = op.get("inputs", [])
        output = op.get("output", "")
        params = op.get("parameters", {})
        if strategy == "GoldenGate":
            enzyme = params.get("enzyme", "BsaI")
            dest, prog = "thermocycler1A", "main/GG1"
            note = (
                f"DNA Mix: 5 uL per fragment. "
                f"Reaction: 7 uL ddH2O / 1 uL T4 ligase buffer / 1 uL DNA Mix / "
                f"0.5 uL T4 ligase / 0.5 uL {enzyme}. Never let enzymes warm up!"
            )
        else:
            dest, prog = "thermocycler1B", "main/GIB2"
            note = (
                "DNA Mix: 5 uL per fragment. "
                "Reaction: 4 uL ddH2O / 1 uL DNA Mix / 5 uL 2X Gibson Mix. "
                "Never let enzymes warm up!"
            )
        rows.append(_tsv_row({
            "section": f"{thread}: Assemble",
            "label": thread,
            "fragments": ", ".join(inputs),
            "product": output,
            "destination": dest,
            "program": prog,
            "notes": note,
        }))

    # RestrictionLigation rows
    for op in rxn_lig_ops:
        params = op.get("parameters", {})
        inputs = op.get("inputs", [])
        enzyme = params.get("enzyme", "EcoRI-HF + BamHI-HF")
        rows.append(_tsv_row({
            "section": f"{thread}: RestrictionLigation",
            "label": thread,
            "fragments": ", ".join(inputs),
            "product": op.get("output", ""),
            "destination": "thermocycler1A -> bench",
            "program": f"digest 37C 1hr ({enzyme}); gel-purify; T4 ligate 16C O/N",
            "notes": (
                f"Digest vector + insert in parallel with {enzyme} (NEB rCutSmart, "
                f"20 uL each); heat-inactivate 65C 20min; gel-purify both; ligate "
                f"3:1 insert:vector with T4 ligase 16C overnight. Optionally rSAP-"
                f"treat vector to reduce self-ligation. Never let enzymes warm up!"
            ),
        }))

    # TypeIIS oligo cloning rows
    for op in typeiis_ops:
        params = op.get("parameters", {})
        inputs = op.get("inputs", [])
        enzyme = params.get("enzyme", "BbsI-HF")
        rows.append(_tsv_row({
            "section": f"{thread}: TypeIISOligoCloning",
            "label": thread,
            "fragments": ", ".join(inputs),
            "product": op.get("output", ""),
            "destination": "thermocycler1A",
            "program": "anneal -> [37C/16C]x30 -> 60C 5min -> 80C 5min",
            "notes": (
                f"Phase 1: anneal top+bottom oligos (95C 5min, slow cool). "
                f"Phase 2: 1 uL diluted duplex + vector + {enzyme} + T4 ligase "
                f"in 10 uL. Never let enzymes warm up!"
            ),
        }))

    # Transform rows
    for op in transform_ops:
        params = op.get("parameters", {})
        inputs = op.get("inputs", [])
        product = inputs[0] if inputs else op.get("output", "")
        strain = params.get("cells", "Mach1")
        antibiotic = params.get("selection", "Amp")
        temp = params.get("temperature_c", 37)
        rescue = _needs_rescue(antibiotic)
        note = (
            f"Rescue required: after transformation add 200 uL 2YT, shake before plating."
            if rescue else ""
        )
        rows.append(_tsv_row({
            "section": f"{thread}: Transform",
            "label": f"{thread}1",
            "product": product,
            "strain": strain,
            "antibiotic": antibiotic,
            "incubation": f"{temp}°C",
            "notes": note,
        }))
        rows.append(_tsv_row({
            "section": f"{thread}: Plate",
            "label": f"{thread}1",
            "product": product,
            "strain": strain,
            "antibiotic": antibiotic,
            "incubation": f"{temp}°C",
            "notes": (
                f"Plate 100 uL on {antibiotic} plate; spread; incubate upside down at "
                f"{temp}°C, no shaking."
            ),
        }))

    # Pick rows
    for idx, op in enumerate(transform_ops):
        params = op.get("parameters", {})
        inputs = op.get("inputs", [])
        product = inputs[0] if inputs else op.get("output", "")
        strain = params.get("cells", "Mach1")
        antibiotic = params.get("selection", "Amp")
        plan = transform_plans[idx] if idx < len(transform_plans) else {"n": 2, "labels": ["A", "B"]}
        n = plan["n"]
        labels = plan["labels"]
        label_str = ", ".join(f"{thread}1{ch}" for ch in labels)
        rows.append(_tsv_row({
            "section": f"{thread}: Pick",
            "label": f"{thread}1",
            "product": product,
            "strain": strain,
            "antibiotic": antibiotic,
            "incubation": "37°C",
            "notes": f"Pick {n} colonies. Label as {label_str}.",
        }))

    # Miniprep rows
    for idx, op in enumerate(transform_ops):
        inputs = op.get("inputs", [])
        product = inputs[0] if inputs else op.get("output", "")
        plan = transform_plans[idx] if idx < len(transform_plans) else {"labels": ["A", "B"]}
        labels = plan["labels"]
        label_str = " / ".join(f"{thread}1{ch}" for ch in labels)
        product_str = " / ".join(f"{product}-{ch}" for ch in labels)
        loc_str = " / ".join(
            (f"box{thread}/{chr(67 + i)}1" if i < 24 else f"box{thread}/overflow{i}")
            for i in range(len(labels))
        )
        rows.append(_tsv_row({
            "section": f"{thread}: Miniprep",
            "label": label_str,
            "product": product_str,
            "product_location": loc_str,
            "notes": "Write short label on top of tube; full name on side.",
        }))

    # Sequencing rows
    for idx, op in enumerate(transform_ops):
        inputs = op.get("inputs", [])
        product = inputs[0] if inputs else op.get("output", "")
        plan = transform_plans[idx] if idx < len(transform_plans) else {"labels": ["A", "B"]}
        labels = plan["labels"]
        label_str = " / ".join(f"{thread}1{ch}" for ch in labels)
        product_str = " / ".join(f"{product}-{ch}" for ch in labels)
        primer_names = " / ".join(sp["name"] for sp in seq_primers)
        primer_locs = " / ".join(sp["location"] for sp in seq_primers)
        seq_note = (
            f"6 uL ddH2O / 4 uL miniprep DNA / 3 uL primer (2.66 uM). "
            f"Submit to sequencing facility. (Confirms plasmid; does NOT "
            f"confirm the genomic edit — see EditVerification section.)"
        )
        if protospacer:
            seq_note += f" Expected read should contain: {protospacer} + scaffold + terminator."
        rows.append(_tsv_row({
            "section": f"{thread}: Sequencing (plasmid)",
            "label": label_str,
            "product": product_str,
            "primer1": primer_names,
            "primer1_location": primer_locs,
            "notes": seq_note,
        }))

    # EditVerification row (genomic Sanger after the experiment)
    if verify_edit_result:
        ve = verify_edit_result
        rows.append(_tsv_row({
            "section": f"{thread}: EditVerification (genomic)",
            "label": f"{thread}_genomic",
            "primer1": "verify_F (genomic)",
            "primer1_location": "oligos1/M1",
            "primer2": "verify_R (genomic)",
            "primer2_location": "oligos1/M2",
            "template": "edited gDNA",
            "product": f"{ve.get('amplicon_length','?')} bp amplicon",
            "program": "Q5 35 cycles, 65C anneal",
            "notes": (
                f"PCR genomic locus, gel-confirm ~{ve.get('amplicon_length','?')} bp band, "
                f"Sanger sequence, upload .ab1 + amplicon to Synthego ICE "
                f"(cut at offset {ve.get('cut_offset_in_amplicon','?')} in amplicon). "
                f"Run interpret_ice_tide on the result."
            ),
        }))

    # DirectSynthesis rows
    for op in direct_ops:
        inputs = op.get("inputs", [])
        output = op.get("output", "")
        construct = inputs[0] if inputs else output
        rows.append(_tsv_row({
            "section": f"{thread}: DirectSynthesis",
            "label": f"{thread}1",
            "product": output,
            "notes": f"Order {construct} from Twist/IDT. Resuspend per vendor instructions.",
        }))

    return "\n".join(rows)


class LabSheet:
    """
    Description:
        Renders a bench-ready lab sheet from a structured construction file in two
        formats: a plain-text step-by-step protocol (matching the professor's
        LabPlanner format) and a tab-separated (TSV) table suitable for Excel.

        Both formats cover the same information: sample labels, primer 1 & 2 names
        and plate locations, template, fragment list, final product, destination
        machine, thermocycler program, and any relevant notes.

        Supports PCR, GoldenGate, Gibson, Transform, and DirectSynthesis steps.
        PCR steps automatically generate Gel+DpnI and Zymo cleanup sections.

    Input:
        construction_record (dict): the structured construction file dict produced
                                    by create_construction_file. Must contain at
                                    minimum: construct_name, assembly_strategy,
                                    and operations (list of step dicts).
        thread (str): single-letter thread identifier (default "A"). Matches the
                      letter prefix used in the lab sheet (A1, A2, etc.).
        include_notes (bool): if True, append any notes at the bottom. default True.

    Output:
        dict with keys:
            - construct_name: name of the construct
            - assembly_strategy: assembly method used
            - lab_sheet_text: full plain-text bench protocol (step-by-step)
            - lab_sheet_tsv: tab-separated table of every step (paste into Excel)
            - step_count: total number of sections rendered

    Tests:
        - Case:
            Input: construction_record with one PCR step
            Expected Output: "lab_sheet_text" contains "PCR" and "samples:" and "source:"
            Description: PCR step renders samples and source tables.
        - Case:
            Input: construction_record with GoldenGate step
            Expected Output: "lab_sheet_text" contains "Assemble" and "DNA Mix:" and "source:" and "samples:"
            Description: GoldenGate step renders DNA mix, reaction volumes, source and samples with colons.
        - Case:
            Input: construction_record with Gibson step
            Expected Output: "lab_sheet_tsv" contains "thermocycler1B" and "main/GIB2" and "Gibson Mix"
            Description: TSV format captures destination, program, and reaction notes for Gibson.
        - Case:
            Input: construction_record with PCR step
            Expected Output: "lab_sheet_tsv" first line equals tab-joined _TSV_HEADERS
            Description: TSV starts with the correct header row.
        - Case:
            Input: construction_record={}
            Expected Exception: ValueError
            Description: empty record raises ValueError.
    """

    def initiate(self) -> None:
        pass

    def run(
        self,
        construction_record: dict,
        thread: str = "A",
        include_notes: bool = True,
        editing_efficiency: float | None = None,
        colony_preset: str | None = None,
        desired_clones: int = 1,
        confidence: float = 0.95,
        sequencing_primers: list[dict] | None = None,
        protospacer: str | None = None,
        verification_reference: str | None = None,
        nuclease: str = "cas9",
        pam: str | None = None,
        delivery: str = "plasmid",
        outcome: str = "nhej",
        location_overrides: dict | None = None,
    ) -> dict:
        if not construction_record:
            raise ValueError("construction_record must not be empty.")

        construct_name = construction_record.get("construct_name", "Unknown")
        assembly_strategy = construction_record.get("assembly_strategy", "Unknown")
        operations = construction_record.get("operations", [])
        notes = construction_record.get("notes", "")
        parts_map = _build_parts_map(construction_record.get("parts", []))

        # Auto-sniff the protospacer from the construction record's parts
        # if the caller didn't supply one. The full-workflow tool builds a
        # part named like "<vector>_annealed_guide_insert" whose sequence
        # is "<20bp protospacer><sgRNA scaffold tail>". When that pattern
        # is present we lift out the first 20 nt as the protospacer so
        # predict_editing_efficiency runs even when the upstream tool
        # didn't thread protospacer through to lab_sheet explicitly.
        if protospacer is None:
            for part in (construction_record.get("parts") or []):
                pname = (part.get("name") or "").lower()
                if "guide" in pname and "insert" in pname:
                    seq = (part.get("sequence") or "").upper()
                    # Sequence must look like real DNA and be at least 20 nt
                    if len(seq) >= 20 and set(seq[:20]) <= set("ACGT"):
                        protospacer = seq[:20]
                        break
        # Same for PAM — agent rarely passes it, but the parts list
        # sometimes carries the reference around the cut site. We don't
        # try to extract that automatically (high error rate); we just
        # default pam to "NGG" so predict_editing_efficiency can fire on
        # Cas9 guides without the user spelling out the PAM.
        if protospacer and pam is None and nuclease.lower() == "cas9":
            pam = "NGG"
        if protospacer and pam is None and nuclease.lower() == "cas12a":
            # Use TTTA (most common Cas12a PAM, real V base) so the
            # predictor's PAM bonus actually fires.
            pam = "TTTA"

        # Best-effort pre-experiment editing efficiency prediction.
        # Surfaced in the output so the user knows what to expect at the
        # bench BEFORE running ICE/TIDE on real data.
        predicted_eff = _maybe_predict_efficiency(
            protospacer, pam, nuclease, delivery, outcome,
        )
        # If we got a prediction, also feed it into the colony plan so
        # colony counts are guide-specific instead of preset-generic
        # (only when caller didn't already pin editing_efficiency).
        if predicted_eff and editing_efficiency is None and colony_preset is None:
            editing_efficiency = predicted_eff["on_target_efficiency_pct"] / 100.0

        # If the caller passes a genomic reference + protospacer, run
        # verify_edit. The result is rendered as a dedicated post-edit
        # GENOMIC verification section — distinct from the plasmid
        # Sequencing step (which sequences the assembled vector with
        # L4440 or a user-supplied vector primer).
        verify_edit_result: dict | None = None
        if protospacer and verification_reference:
            try:
                verify_edit_result = verify_edit(
                    protospacer=protospacer,
                    reference=verification_reference,
                    nuclease=nuclease,
                )
            except Exception:
                verify_edit_result = None

        pcr_ops = [op for op in operations if op.get("step_type") == "PCR"]
        assemble_ops = [op for op in operations if op.get("step_type") in ("GoldenGate", "Gibson")]
        transform_ops = [op for op in operations if op.get("step_type") == "Transform"]
        direct_ops = [op for op in operations if op.get("step_type") == "DirectSynthesis"]
        crispr_ops = [op for op in operations if op.get("step_type") == "CRISPRDelivery"]
        typeiis_ops = [op for op in operations if op.get("step_type") == "TypeIISOligoCloning"]
        rxn_lig_ops = [op for op in operations if op.get("step_type") == "RestrictionLigation"]

        # Safety net: any step_type not in a known bucket gets explicit
        # placeholder rendering instead of being silently dropped.
        # We split this into two classes:
        #   - "generic_assembly_ops": agent-flexible synonyms ("Assemble",
        #     "Ligate", "Anneal") that LOOK like an assembly step. Render
        #     these BEFORE Transform so the biology sequence stays intact.
        #   - "unknown_ops": truly unrecognized step types (no obvious slot).
        #     Render last with a WARNING so the user spots the gap.
        _known = {
            "PCR", "GoldenGate", "Gibson", "Transform", "DirectSynthesis",
            "TypeIISOligoCloning", "RestrictionLigation", "CRISPRDelivery",
        }
        _GENERIC_ASSEMBLY_SYNONYMS = {
            "assemble", "assembly", "ligate", "ligation", "anneal",
            "annealedoligocloning", "oligocloning", "annealed_oligo_cloning",
        }

        def _is_generic_assembly(op: dict) -> bool:
            t = (op.get("step_type") or "").lower().replace(" ", "").replace("-", "_")
            return t in _GENERIC_ASSEMBLY_SYNONYMS

        generic_assembly_ops = [
            op for op in operations
            if op.get("step_type") not in _known and _is_generic_assembly(op)
        ]
        unknown_ops = [
            op for op in operations
            if op.get("step_type") not in _known and not _is_generic_assembly(op)
        ]

        # Decide colony-screening burden once per Transform op via the
        # colony_calculator tool, instead of hardcoding "pick 2 colonies".
        transform_plans = [
            _build_colony_plan(
                op, editing_efficiency, colony_preset, desired_clones, confidence,
                is_crispr=bool(crispr_ops),
            )
            for op in transform_ops
        ]

        sections = []

        if pcr_ops:
            sections.append(_format_pcr_section(pcr_ops, thread))
            sections.append(_format_gel_dpni_section(pcr_ops, thread, parts_map))
            sections.append(_format_zymo_section(pcr_ops, thread))

        for op in assemble_ops:
            sections.append(_format_assemble_section(op, thread))

        # TypeIISOligoCloning is the assembly step for sgRNA cloning into
        # BbsI/BsmBI-cut vectors (pX330, Lenti-Guide, etc.). Without
        # rendering it the workflow would jump straight to Transform with
        # no annealed-and-ligated insert in hand — biologically invalid.
        for op in typeiis_ops:
            sections.append(_format_typeiis_oligo_section(op, thread))

        for op in rxn_lig_ops:
            sections.append(_format_restriction_ligation_section(op, thread))

        # Generic assembly steps with a synonym name (e.g. agent passes
        # step_type="Assemble" instead of "TypeIISOligoCloning"). These
        # MUST render before Transform — otherwise the lab sheet reads
        # transform-first which is biologically wrong.
        for op in generic_assembly_ops:
            # Pick the closest registry recipe based on parameters.enzyme
            params = op.get("parameters", {})
            enzyme = (params.get("enzyme") or "").lower()
            if any(e in enzyme for e in ("bsai", "bsmbi", "bbsi", "swai", "bcli")):
                # Type IIS-class — render via the typeiis section
                sections.append(_format_typeiis_oligo_section(op, thread))
            elif "gibson" in enzyme:
                op_copy = {**op, "step_type": "Gibson"}
                sections.append(_format_assemble_section(op_copy, thread))
            else:
                # Fall back to GoldenGate-shape rendering — most agent-
                # generated "Assemble" steps are some form of restriction-
                # plus-ligase reaction.
                op_copy = {**op, "step_type": "GoldenGate"}
                sections.append(_format_assemble_section(op_copy, thread))

        for op, plan in zip(transform_ops, transform_plans):
            sections.append(_format_transform_section(op, thread))
            sections.append(_format_plate_section(op, thread))
            sections.append(_format_pick_section(op, thread, plan))
            sections.append(_format_miniprep_section(op, thread, plan))
            sections.append(_format_sequencing_section(op, thread, plan, sequencing_primers, protospacer))

        for op in crispr_ops:
            sections.append(_format_crispr_delivery_section(op, thread))

        for op in direct_ops:
            sections.append(_format_direct_synthesis_section(op, thread))

        for op in unknown_ops:
            sections.append(_format_unknown_step_section(op, thread))

        # Post-edit GENOMIC verification — appended after all cloning
        # steps because it happens after the experiment, on harvested
        # cells, not on the assembled plasmid.
        if verify_edit_result:
            sections.append(
                _format_edit_verification_section(
                    thread, verify_edit_result, nuclease, delivery,
                )
            )
        elif protospacer:
            # No verification_reference was provided, but we know there's a
            # protospacer in this design (CRISPR workflow). Emit a placeholder
            # so the user knows that post-edit genomic Sanger is the next
            # step after their plasmid is confirmed — and how to get specific
            # primers when they're ready.
            sections.append(
                _format_edit_verification_placeholder(
                    thread, protospacer, nuclease,
                )
            )

        if include_notes and notes:
            sections.append(f"note:\n{notes}")

        lab_sheet_text = "\n\n".join(sections)
        lab_sheet_tsv = _build_tsv(
            operations, thread, notes, transform_plans, parts_map,
            sequencing_primers, verify_edit_result, protospacer,
        )

        # Apply user-supplied location overrides as a final pass over both
        # the text and TSV outputs. Examples:
        #   {"oligos1": "Freezer 3 / Rack B", "boxA": "Bench drawer 2"}
        # Replacements are longest-token-first so "boxA1" beats "boxA".
        if location_overrides:
            sorted_keys = sorted(location_overrides.keys(), key=len, reverse=True)
            for token in sorted_keys:
                replacement = str(location_overrides[token])
                lab_sheet_text = lab_sheet_text.replace(token, replacement)
                lab_sheet_tsv = lab_sheet_tsv.replace(token, replacement)

        # Build a list of canonical protocol sources (one entry per step type
        # that fired). Each entry includes the manufacturer/method paper its
        # volumes came from, so the output is fully auditable.
        protocol_keys: list[str] = []
        if pcr_ops:
            protocol_keys.append("pcr_q5")
        for op in assemble_ops:
            if op.get("step_type") == "GoldenGate":
                protocol_keys.append("goldengate_bsai")
            else:
                protocol_keys.append("gibson_neb_2x")
        if typeiis_ops:
            protocol_keys.append("typeiis_oligo_cloning")
        if rxn_lig_ops:
            protocol_keys.append("restriction_ligation")
        if transform_ops:
            protocol_keys.extend([
                "transformation_chemical",
                "miniprep_zymo",
                "sanger_sequencing",
            ])
        for op in crispr_ops:
            method = (op.get("parameters", {}).get("method") or "rnp_assembly").lower().strip()
            key = _CRISPR_DELIVERY_PROTOCOLS.get(method)
            if key:
                protocol_keys.append(key)
        # de-duplicate while preserving order
        seen: set[str] = set()
        protocol_keys = [k for k in protocol_keys if not (k in seen or seen.add(k))]

        protocol_sources = [protocol_source_record(k) for k in protocol_keys]

        # Surface every underlying citation as a flat list too, matching the
        # citations field shape used by the other tools (rank_guides,
        # predict_offtargets, colony_calculator, etc.).
        citation_keys: list[str] = []
        for k in protocol_keys:
            proto = PROTOCOLS[k]
            citation_keys.append(proto.source)
            citation_keys.extend(proto.extra_sources)
        # de-dupe citations too
        seen_cit: set[str] = set()
        citation_keys = [c for c in citation_keys if not (c in seen_cit or seen_cit.add(c))]

        # Merge in colony_calculator citations (Kim 2014 / Paquet 2016 /
        # Jiang 2013 depending on the preset/efficiency used).
        protocol_citations = format_citations(cites(*citation_keys)) if citation_keys else []
        seen_labels = {c["label"] for c in protocol_citations}
        for plan in transform_plans:
            for c in plan.get("citations", []):
                if c["label"] not in seen_labels:
                    protocol_citations.append(c)
                    seen_labels.add(c["label"])

        # Surface the colony plan so downstream tooling can audit it.
        colony_plan_summary = [
            {
                "transform_index": i,
                "colonies_to_pick": p["n"],
                "labels": [f"{thread}1{ch}" for ch in p["labels"]],
                "efficiency_used": p["efficiency"],
                "preset": p["preset"],
                "rationale": p["recommendation"],
            }
            for i, p in enumerate(transform_plans)
        ]

        # Build protocols.io enrichment links — one per protocol used.
        # We don't fetch from protocols.io (free-text protocol prose is
        # easy for an LLM to misread), but the user can click through to
        # see real, peer-shared variants beyond the curated registry.
        protocols_io_links = [
            {
                "protocol": PROTOCOLS[k].name,
                "search_url": _protocols_io_search_url(PROTOCOLS[k].name),
            }
            for k in protocol_keys
        ]

        return {
            "construct_name": construct_name,
            "assembly_strategy": assembly_strategy,
            "lab_sheet_text": lab_sheet_text,
            "lab_sheet_tsv": lab_sheet_tsv,
            "step_count": len(sections),
            "protocol_sources": protocol_sources,
            "citations": protocol_citations,
            "verify_before_use": _VERIFY_DISCLAIMER,
            "colony_plan": colony_plan_summary,
            "protocols_io_links": protocols_io_links,
            "predicted_editing_efficiency": (
                {
                    "on_target_efficiency_pct": predicted_eff["on_target_efficiency_pct"],
                    "confidence_range": predicted_eff["confidence_range"],
                    "interpretation": predicted_eff["interpretation"],
                    "warnings": predicted_eff.get("warnings", []),
                    "delivery": predicted_eff.get("delivery"),
                    "outcome": predicted_eff.get("outcome"),
                }
                if predicted_eff
                else None
            ),
            "verify_edit_summary": (
                {
                    "protospacer": verify_edit_result["protospacer"],
                    "cut_position": verify_edit_result["cut_position"],
                    "amplicon_length": verify_edit_result["amplicon_length"],
                    "forward_primer": verify_edit_result["forward_primer"],
                    "reverse_primer": verify_edit_result["reverse_primer"],
                    "forward_primer_tm": verify_edit_result["forward_primer_tm"],
                    "reverse_primer_tm": verify_edit_result["reverse_primer_tm"],
                }
                if verify_edit_result
                else None
            ),
        }


_instance = LabSheet()
_instance.initiate()
lab_sheet = _instance.run
