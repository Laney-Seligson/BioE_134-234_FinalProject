from __future__ import annotations


def _format_pcr_step(op: dict, idx: int) -> str:
    primers = op.get("primers", [])
    template = op.get("template", "unknown template")
    output = op.get("output", f"product_{idx}")
    lines = [
        f"### Step {idx}: PCR — {output}",
        f"- Template: {template}",
        f"- Primers: {', '.join(primers) if primers else 'see construction file'}",
        f"- Output: {output}",
        "",
        "**Reagents (25 µL reaction):**",
        "| Reagent | Volume |",
        "|---------|--------|",
        "| 2× Q5 Master Mix | 12.5 µL |",
        "| Forward primer (10 µM) | 1.25 µL |",
        "| Reverse primer (10 µM) | 1.25 µL |",
        f"| Template ({template}) | 1 µL |",
        "| Nuclease-free water | to 25 µL |",
        "",
        "**Thermocycler program:** 98°C 30 s → [98°C 10 s / 60°C 30 s / 72°C 30 s/kb] × 30 → 72°C 2 min",
        "",
    ]
    return "\n".join(lines)


def _format_goldengate_step(op: dict, idx: int) -> str:
    enzyme = op.get("enzyme", "BsaI")
    inputs = op.get("inputs", [])
    output = op.get("output", "assembled_plasmid")
    lines = [
        f"### Step {idx}: Golden Gate Assembly — {output}",
        f"- Enzyme: {enzyme}",
        f"- Inputs: {', '.join(inputs) if inputs else 'see construction file'}",
        f"- Output: {output}",
        "",
        "**Reagents (20 µL reaction):**",
        "| Reagent | Volume |",
        "|---------|--------|",
        f"| {enzyme}-HF (NEB) | 1 µL |",
        "| T4 DNA Ligase (NEB) | 1 µL |",
        "| T4 DNA Ligase Buffer (10×) | 2 µL |",
        "| Each PCR product (~100 ng) | variable |",
        "| Nuclease-free water | to 20 µL |",
        "",
        "**Thermocycler program:** [37°C 5 min / 16°C 5 min] × 25 → 60°C 5 min → 4°C hold",
        "",
    ]
    return "\n".join(lines)


def _format_gibson_step(op: dict, idx: int) -> str:
    inputs = op.get("inputs", [])
    output = op.get("output", "assembled_plasmid")
    lines = [
        f"### Step {idx}: Gibson Assembly — {output}",
        f"- Inputs: {', '.join(inputs) if inputs else 'see construction file'}",
        f"- Output: {output}",
        "",
        "**Reagents (20 µL reaction):**",
        "| Reagent | Volume |",
        "|---------|--------|",
        "| 2× Gibson Assembly Master Mix (NEB) | 10 µL |",
        "| Each fragment (~100 ng each) | variable |",
        "| Nuclease-free water | to 20 µL |",
        "",
        "**Incubation:** 50°C for 60 minutes, then place on ice.",
        "",
    ]
    return "\n".join(lines)


def _format_transform_step(op: dict, idx: int) -> str:
    strain = op.get("strain", "DH5alpha")
    selection = op.get("selection", "antibiotic selection")
    temperature = op.get("temperature_c", 37)
    lines = [
        f"### Step {idx}: Transformation — {strain}",
        f"- Strain: {strain}",
        f"- Selection: {selection}",
        f"- Growth temperature: {temperature}°C",
        "",
        "**Protocol:**",
        "1. Thaw 50 µL competent cells on ice.",
        "2. Add 2 µL of assembly reaction; mix gently.",
        "3. Incubate on ice 30 min.",
        "4. Heat shock 42°C for 30 s.",
        "5. Return to ice 2 min.",
        "6. Add 950 µL SOC media; recover 37°C 1 h shaking.",
        f"7. Plate on LB + {selection}; incubate {temperature}°C overnight.",
        "",
    ]
    return "\n".join(lines)


_STEP_FORMATTERS = {
    "PCR": _format_pcr_step,
    "GoldenGate": _format_goldengate_step,
    "Gibson": _format_gibson_step,
    "Transform": _format_transform_step,
}


class LabSheet:
    """
    Description:
        Renders a bench-ready markdown lab sheet from a structured construction
        file. The output is a human-readable document a researcher can print and
        follow at the bench — listing each cloning step with reagent tables,
        volumes, and thermocycler programs.

        Currently fully implemented for PCR, GoldenGate, Gibson, and Transform
        steps. DirectSynthesis and Digest steps return a placeholder noting that
        the fragment should be ordered or the protocol should be followed manually.

        Input comes from the dict returned by create_construction_file
        (the "structured_construction_file" key).

    Input:
        construction_record (dict): the structured construction file dict produced
                                    by create_construction_file. Must contain at
                                    minimum: construct_name, assembly_strategy,
                                    and operations (list of step dicts).
        include_notes (bool): if True, append any notes from the construction
                              file at the bottom. default True.

    Output:
        dict with keys:
            - construct_name: name of the construct
            - assembly_strategy: assembly method used
            - lab_sheet_md: full markdown text of the bench protocol
            - step_count: total number of steps rendered

    Tests:
        - Case:
            Input: minimal construction_record with one PCR step
            Expected Output: "lab_sheet_md" contains "PCR" and reagent table
            Description: PCR step renders reagent table and thermocycler program.
        - Case:
            Input: construction_record with GoldenGate step
            Expected Output: "lab_sheet_md" contains "Golden Gate" and enzyme name
            Description: GoldenGate step renders correct enzyme and protocol.
        - Case:
            Input: construction_record={}
            Expected Exception: ValueError
            Description: empty record raises ValueError.
    """

    def initiate(self) -> None:
        pass

    def run(self, construction_record: dict, include_notes: bool = True) -> dict:
        if not construction_record:
            raise ValueError("construction_record must not be empty.")

        construct_name = construction_record.get("construct_name", "Unknown")
        assembly_strategy = construction_record.get("assembly_strategy", "Unknown")
        operations = construction_record.get("operations", [])
        notes = construction_record.get("notes", "")

        sections = [
            f"# Lab Sheet: {construct_name}",
            f"**Assembly strategy:** {assembly_strategy}",
            f"**Steps:** {len(operations)}",
            "",
            "---",
            "",
        ]

        for idx, op in enumerate(operations, start=1):
            step_type = op.get("type", "Unknown")
            formatter = _STEP_FORMATTERS.get(step_type)
            if formatter:
                sections.append(formatter(op, idx))
            else:
                # TODO: add formatters for Digest and DirectSynthesis steps
                sections.append(
                    f"### Step {idx}: {step_type}\n"
                    f"_Protocol not yet implemented for {step_type} steps. "
                    f"Refer to the construction file for details._\n"
                )

        if include_notes and notes:
            sections += ["---", "", f"**Notes:** {notes}", ""]

        lab_sheet_md = "\n".join(sections)

        return {
            "construct_name": construct_name,
            "assembly_strategy": assembly_strategy,
            "lab_sheet_md": lab_sheet_md,
            "step_count": len(operations),
        }


_instance = LabSheet()
_instance.initiate()
lab_sheet = _instance.run
