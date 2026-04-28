from __future__ import annotations

from modules.labsheet_tools.tools._protocols import (
    PROTOCOLS,
    reagent_block,
    protocol_source_record,
)
from modules.labsheet_tools.tools.colony_calculator import colony_calculator
from modules.labsheet_tools.tools.verify_edit import verify_edit
from modules.crispr_tools.tools.citations import cites, format_citations

_RESCUE_ANTIBIOTICS = {"Spec", "spectinomycin", "Spc"}

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
    """Find the position where the primer (or its 3'-end annealing region)
    binds the template. Returns the 0-based start in `template`, or None."""
    if not primer_seq or not template:
        return None
    primer_seq = primer_seq.upper()
    template = template.upper()
    search = template + (template if circular else "")
    # Try progressively shorter 3' anneal lengths. Real primers often have
    # 5' overhangs (Gibson, Golden Gate BsaI sites) that don't anneal.
    for anneal_len in (len(primer_seq), 25, 22, 20, 18, 15):
        if anneal_len > len(primer_seq):
            continue
        anneal = primer_seq[-anneal_len:]
        idx = search.find(anneal)
        if idx != -1:
            # `idx` is the start of the anneal region; the primer's 5' end
            # would be at idx - (len(primer_seq) - anneal_len).
            return idx - (len(primer_seq) - anneal_len)
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
        # Add primer 5' overhang contribution beyond template-anneal:
        # caller expects total amplified length, which already equals
        # (rev_end - fwd_pos) since fwd_pos was anchored to the primer's
        # 5' end in `_locate_primer`.
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


def _guess_preset(cells: str) -> str | None:
    """Map a strain/cell-line name to a colony_calculator preset, or None
    to fall back to _DEFAULT_CLONING_SUCCESS."""
    if not cells:
        return None
    up = cells.upper()
    for hint in _MAMMALIAN_HINTS:
        if hint in up:
            return "cas9_plasmid_mammalian"
    return None  # E. coli / yeast default → use cloning success rate


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

    for i, op in enumerate(ops, start=1):
        label = f"{thread}{i}"
        p = op.get("parameters", {})
        fwd = p.get("forward_primer", "")
        rev = p.get("reverse_primer", "")
        tmpl = p.get("template", "")
        product = op.get("output", "")
        rows.append([label, fwd, rev, tmpl, product])

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

    return "\n".join([
        f"{thread}: PCR",
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


def _format_restriction_ligation_section(op: dict, thread: str) -> str:
    """Render a RestrictionLigation step (vector + insert, one or two
    restriction enzymes, then T4 ligation). Without this section the
    workflow would jump straight from PCR to Transform — the same
    biological-completeness gap that TypeIISOligoCloning had."""
    proto_key = "restriction_ligation"
    proto = PROTOCOLS[proto_key]
    params = op.get("parameters", {})
    inputs = op.get("inputs", [])
    output = op.get("output", "")
    enzyme = params.get("enzyme", "EcoRI-HF + BamHI-HF")

    backbone = inputs[0] if len(inputs) > 0 else "vector"
    insert = inputs[1] if len(inputs) > 1 else "insert"

    src_rows = [
        [backbone, f"box{thread}/A1"],
        [insert, f"box{thread}/B1"],
    ]
    source_table = _col_table(["dna", "location"], src_rows)
    samples_table = _col_table(
        ["label", "fragments", "product"],
        [[thread, f"{backbone} / {insert}", output]],
    )

    return "\n".join([
        f"{thread}: RestrictionLigation",
        "",
        f"protocol: {proto.name}",
        f"enzyme(s): {enzyme}",
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
        "destination: thermocycler1A (digest) -> bench (ligate)",
        "",
        f"notes:\n{proto.notes}",
        "",
        _ENZYME_NOTE,
    ])


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
    rescue = antibiotic in _RESCUE_ANTIBIOTICS
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
    op: dict, thread: str, plan: dict, seq_primers: list[dict] | None = None,
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

    return "\n".join([
        f"{thread}: Sequencing",
        "",
        "sources:",
        source_table,
        "",
        "Instructions:",
    ] + [
        f"For each plasmid + primer combination, mix the following in an eppendorf tube:"
    ] + [
        f"Primers in this run: {', '.join(sp['name'] for sp in seq_primers)}.",
        "Resuspend each primer to a 100 uM stock; dilute 1:37.6 in ddH2O for a 2.66 uM working stock.",
        "  6 uL ddH2O",
        "  4 uL miniprep DNA (undiluted)",
        "  3 uL primer (2.66 uM)",
        f'Label the tops of the tubes with the label (e.g. "{thread}1{plan["labels"][0]}").',
        "Take the sequencing reactions and order form to: 237 Stanley Hall (second floor cold room).",
        "" if any(sp["name"] == "L4440" for sp in seq_primers) else
        "Note: using guide-specific verification primers (verify_edit) — sequence WILL span the edit site.",
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
        rescue = antibiotic in _RESCUE_ANTIBIOTICS
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
        primer_note = (
            "(guide-specific verification primers from verify_edit) "
            if not any(sp["name"] == "L4440" for sp in seq_primers) else ""
        )
        rows.append(_tsv_row({
            "section": f"{thread}: Sequencing",
            "label": label_str,
            "product": product_str,
            "primer1": primer_names,
            "primer1_location": primer_locs,
            "notes": (
                f"6 uL ddH2O / 4 uL miniprep DNA / 3 uL primer (2.66 uM). "
                f"{primer_note}Submit to 237 Stanley Hall."
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
    ) -> dict:
        if not construction_record:
            raise ValueError("construction_record must not be empty.")

        construct_name = construction_record.get("construct_name", "Unknown")
        assembly_strategy = construction_record.get("assembly_strategy", "Unknown")
        operations = construction_record.get("operations", [])
        notes = construction_record.get("notes", "")
        parts_map = _build_parts_map(construction_record.get("parts", []))

        # If the caller passes a protospacer + reference, auto-call
        # verify_edit to design guide-specific Sanger primers that span
        # the cut site. Explicit `sequencing_primers` still wins (explicit
        # > auto > L4440 fallback). This closes the loop the SKILL.md
        # workflow describes: rank_guides -> verify_edit -> lab_sheet.
        verify_edit_result: dict | None = None
        if sequencing_primers is None and protospacer and verification_reference:
            try:
                verify_edit_result = verify_edit(
                    protospacer=protospacer,
                    reference=verification_reference,
                    nuclease=nuclease,
                )
                sequencing_primers = [
                    {
                        "name": "verify_F",
                        "sequence": verify_edit_result["forward_primer"],
                        "location": "oligos1/M1",
                        "note": (
                            f"verify_edit forward primer; "
                            f"Tm={verify_edit_result['forward_primer_tm']}C"
                        ),
                    },
                    {
                        "name": "verify_R",
                        "sequence": verify_edit_result["reverse_primer"],
                        "location": "oligos1/M2",
                        "note": (
                            f"verify_edit reverse primer; "
                            f"Tm={verify_edit_result['reverse_primer_tm']}C"
                        ),
                    },
                ]
            except Exception:
                # If verify_edit can't design primers (e.g. protospacer
                # not found in reference), fall through to L4440 so the
                # lab sheet still renders. The user gets a generic primer
                # and can re-run with a corrected reference.
                verify_edit_result = None

        pcr_ops = [op for op in operations if op.get("step_type") == "PCR"]
        assemble_ops = [op for op in operations if op.get("step_type") in ("GoldenGate", "Gibson")]
        transform_ops = [op for op in operations if op.get("step_type") == "Transform"]
        direct_ops = [op for op in operations if op.get("step_type") == "DirectSynthesis"]
        crispr_ops = [op for op in operations if op.get("step_type") == "CRISPRDelivery"]
        typeiis_ops = [op for op in operations if op.get("step_type") == "TypeIISOligoCloning"]
        rxn_lig_ops = [op for op in operations if op.get("step_type") == "RestrictionLigation"]

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

        for op, plan in zip(transform_ops, transform_plans):
            sections.append(_format_transform_section(op, thread))
            sections.append(_format_pick_section(op, thread, plan))
            sections.append(_format_miniprep_section(op, thread, plan))
            sections.append(_format_sequencing_section(op, thread, plan, sequencing_primers))

        for op in crispr_ops:
            sections.append(_format_crispr_delivery_section(op, thread))

        for op in direct_ops:
            sections.append(_format_direct_synthesis_section(op, thread))

        if include_notes and notes:
            sections.append(f"note:\n{notes}")

        lab_sheet_text = "\n\n".join(sections)
        lab_sheet_tsv = _build_tsv(
            operations, thread, notes, transform_plans, parts_map, sequencing_primers
        )

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
