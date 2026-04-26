from __future__ import annotations

_RESCUE_ANTIBIOTICS = {"Spec", "spectinomycin", "Spc"}

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


def _format_gel_dpni_section(ops: list[dict], thread: str) -> str:
    rows = []
    for i, op in enumerate(ops, start=1):
        label = f"{thread}{i}"
        product = op.get("output", "")
        rows.append([label, "?", product])

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

    if strategy == "GoldenGate":
        enzyme = params.get("enzyme", "BsaI")
        reaction = (
            "reaction:\n"
            "7 uL ddH2O\n"
            "1 uL T4 DNA ligase buffer\n"
            "1 uL DNA Mix\n"
            "0.5 uL T4 DNA ligase\n"
            f"0.5 uL {enzyme}"
        )
        destination = "thermocycler1A"
        program = "main/GG1"
    else:
        reaction = (
            "reaction:\n"
            "4 uL ddH2O\n"
            "1 uL DNA Mix\n"
            "5 uL 2X Gibson Mix"
        )
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


def _format_pick_section(op: dict, thread: str) -> str:
    params = op.get("parameters", {})
    inputs = op.get("inputs", [])
    strain = params.get("cells", "Mach1")
    antibiotic = params.get("selection", "Amp")
    product = inputs[0] if inputs else op.get("output", "")

    rows = [
        [f"{thread}1", product, strain, antibiotic, "37°C", "2", f"{thread}1A, {thread}1B"]
    ]
    samples_table = _col_table(
        ["source", "product", "strain", "antibiotic", "incubate", "number", "labels"], rows
    )

    return "\n".join([
        f"{thread}: Pick",
        "",
        "samples:",
        samples_table,
    ])


def _format_miniprep_section(op: dict, thread: str) -> str:
    params = op.get("parameters", {})
    inputs = op.get("inputs", [])
    product = inputs[0] if inputs else op.get("output", "")

    rows = [
        [f"{thread}1A", f"{product}-A", f"box{thread}/C1"],
        [f"{thread}1B", f"{product}-B", f"box{thread}/D1"],
    ]
    samples_table = _col_table(["culture", "label", "location"], rows)

    return "\n".join([
        f"{thread}: Miniprep",
        "",
        "samples:",
        samples_table,
        "",
        "note:",
        f"Write the short label (e.g. {thread}1A) on the top of the eppendorf,",
        "put the full name (the label column) on the side of the tube.",
    ])


def _format_sequencing_section(op: dict, thread: str) -> str:
    params = op.get("parameters", {})
    inputs = op.get("inputs", [])
    product = inputs[0] if inputs else op.get("output", "")

    src_rows = [
        [f"{thread}1A", f"box{thread}/C1", f"{product}-A"],
        [f"{thread}1B", f"box{thread}/D1", f"{product}-B"],
        ["L4440", "oligos1/J1", "sequencing primer"],
    ]
    source_table = _col_table(["label", "location", "product"], src_rows)

    return "\n".join([
        f"{thread}: Sequencing",
        "",
        "sources:",
        source_table,
        "",
        "Instructions:",
        "Resuspend the oligo L4440 to the appropriate volume for a 100 uM stock.",
        "In an eppendorf, prepare a 2.66 uM dilute stock of L4440 as:",
        "  487 uL ddH2O",
        "  13.3 uL of 100 uM oligo",
        "For each plasmid listed, mix the following sequencing reactions in an eppendorf tube:",
        "  6 uL ddH2O",
        "  4 uL miniprep DNA (undiluted)",
        "  3 uL oligo (2.66 uM)",
        f'Label the tops of the tubes with the label (e.g. "{thread}1A").',
        "When done, save the L4440 stock at: oligos1/A5",
        "Take the sequencing reactions and order form to: 237 Stanley Hall (second floor cold room)",
    ])


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


def _build_tsv(operations: list[dict], thread: str, notes_str: str) -> str:
    rows = ["\t".join(_TSV_HEADERS)]

    pcr_ops = [op for op in operations if op.get("step_type") == "PCR"]
    assemble_ops = [op for op in operations if op.get("step_type") in ("GoldenGate", "Gibson")]
    transform_ops = [op for op in operations if op.get("step_type") == "Transform"]
    direct_ops = [op for op in operations if op.get("step_type") == "DirectSynthesis"]

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
        rows.append(_tsv_row({
            "section": f"{thread}: Gel and DpnI",
            "label": f"{thread}{i}",
            "product": op.get("output", ""),
            "destination": "thermocycler1A",
            "program": "main/SPE1",
            "notes": '6 uL "1x Load" + 2 uL PCR product per well; BstEII ladder; 0.5 uL DpnI per reaction. Never let enzymes warm up!',
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
    for op in transform_ops:
        params = op.get("parameters", {})
        inputs = op.get("inputs", [])
        product = inputs[0] if inputs else op.get("output", "")
        strain = params.get("cells", "Mach1")
        antibiotic = params.get("selection", "Amp")
        rows.append(_tsv_row({
            "section": f"{thread}: Pick",
            "label": f"{thread}1",
            "product": product,
            "strain": strain,
            "antibiotic": antibiotic,
            "incubation": "37°C",
            "notes": f"Pick 2 colonies. Label as {thread}1A and {thread}1B.",
        }))

    # Miniprep rows
    for op in transform_ops:
        inputs = op.get("inputs", [])
        product = inputs[0] if inputs else op.get("output", "")
        rows.append(_tsv_row({
            "section": f"{thread}: Miniprep",
            "label": f"{thread}1A / {thread}1B",
            "product": f"{product}-A / {product}-B",
            "product_location": f"box{thread}/C1 / box{thread}/D1",
            "notes": "Write short label on top of tube; full name on side.",
        }))

    # Sequencing rows
    for op in transform_ops:
        inputs = op.get("inputs", [])
        product = inputs[0] if inputs else op.get("output", "")
        rows.append(_tsv_row({
            "section": f"{thread}: Sequencing",
            "label": f"{thread}1A / {thread}1B",
            "product": f"{product}-A / {product}-B",
            "primer1": "L4440",
            "primer1_location": "oligos1/J1",
            "notes": (
                "6 uL ddH2O / 4 uL miniprep DNA / 3 uL L4440 2.66 uM. "
                "Submit to 237 Stanley Hall."
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

    def run(self, construction_record: dict, thread: str = "A", include_notes: bool = True) -> dict:
        if not construction_record:
            raise ValueError("construction_record must not be empty.")

        construct_name = construction_record.get("construct_name", "Unknown")
        assembly_strategy = construction_record.get("assembly_strategy", "Unknown")
        operations = construction_record.get("operations", [])
        notes = construction_record.get("notes", "")

        pcr_ops = [op for op in operations if op.get("step_type") == "PCR"]
        assemble_ops = [op for op in operations if op.get("step_type") in ("GoldenGate", "Gibson")]
        transform_ops = [op for op in operations if op.get("step_type") == "Transform"]
        direct_ops = [op for op in operations if op.get("step_type") == "DirectSynthesis"]

        sections = []

        if pcr_ops:
            sections.append(_format_pcr_section(pcr_ops, thread))
            sections.append(_format_gel_dpni_section(pcr_ops, thread))
            sections.append(_format_zymo_section(pcr_ops, thread))

        for op in assemble_ops:
            sections.append(_format_assemble_section(op, thread))

        for op in transform_ops:
            sections.append(_format_transform_section(op, thread))
            sections.append(_format_pick_section(op, thread))
            sections.append(_format_miniprep_section(op, thread))
            sections.append(_format_sequencing_section(op, thread))

        for op in direct_ops:
            sections.append(_format_direct_synthesis_section(op, thread))

        if include_notes and notes:
            sections.append(f"note:\n{notes}")

        lab_sheet_text = "\n\n".join(sections)
        lab_sheet_tsv = _build_tsv(operations, thread, notes)

        return {
            "construct_name": construct_name,
            "assembly_strategy": assembly_strategy,
            "lab_sheet_text": lab_sheet_text,
            "lab_sheet_tsv": lab_sheet_tsv,
            "step_count": len(sections),
        }


_instance = LabSheet()
_instance.initiate()
lab_sheet = _instance.run
