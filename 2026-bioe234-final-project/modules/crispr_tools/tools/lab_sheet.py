from __future__ import annotations

_RESCUE_ANTIBIOTICS = {"Spec", "spectinomycin", "Spc"}

_ENZYME_NOTE = (
    "note:\n"
    "Never let enzymes warm up! Only take the enzyme cooler out of the freezer when you\n"
    "are actively using it, and only take the tubes out of it when actively dispensing."
)


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
        "source",
        source_table,
        "",
        "samples",
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


class LabSheet:
    """
    Description:
        Renders a bench-ready plain-text lab sheet from a structured construction
        file, matching the professor's LabPlanner format. Output is organized into
        lettered sections (A: PCR, A: Assemble, A: Transform) that a researcher
        can print and follow at the bench.

        Supports PCR, GoldenGate, Gibson, Transform, and DirectSynthesis steps.
        PCR steps automatically generate Gel+DpnI and Zymo cleanup sections.
        Source/location fields use placeholder box positions since no live
        inventory is connected.

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
            - lab_sheet_text: full plain-text bench protocol
            - step_count: total number of sections rendered

    Tests:
        - Case:
            Input: construction_record with one PCR step
            Expected Output: "lab_sheet_text" contains "PCR" and "samples:" and "source:"
            Description: PCR step renders samples and source tables.
        - Case:
            Input: construction_record with GoldenGate step
            Expected Output: "lab_sheet_text" contains "Assemble" and "DNA Mix:"
            Description: GoldenGate step renders DNA mix and reaction volumes.
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

        for op in direct_ops:
            sections.append(_format_direct_synthesis_section(op, thread))

        if include_notes and notes:
            sections.append(f"note:\n{notes}")

        lab_sheet_text = "\n\n".join(sections)

        return {
            "construct_name": construct_name,
            "assembly_strategy": assembly_strategy,
            "lab_sheet_text": lab_sheet_text,
            "step_count": len(sections),
        }


_instance = LabSheet()
_instance.initiate()
lab_sheet = _instance.run
