from dataclasses import dataclass
from pathlib import Path
from typing import Optional


@dataclass(frozen=True)
class VectorConfig:
    assembly_method: str
    enzyme: str
    top_overhang: str
    bottom_overhang: str
    u6_requires_5prime_g: bool
    backbone_resource: Optional[str] = None
    cell_strain: str = ""
    selection: str = ""
    notes: str = ""
    source: str = ""


@dataclass(frozen=True)
class LocalReference:
    organism: str
    resource_name: str
    path: Path
    description: str


TOOL_DIR = Path(__file__).resolve().parent
CRISPR_TOOLS_DIR = TOOL_DIR.parent
BUNDLED_DATA_DIR = CRISPR_TOOLS_DIR / "data"


VECTOR_CONFIG = {
    "px330": VectorConfig(
        assembly_method="annealed_oligos",
        enzyme="BbsI",
        top_overhang="CACC",
        bottom_overhang="AAAC",
        u6_requires_5prime_g=True,
        notes="Standard mammalian Cas9 cloning vector.",
        source=(
            "pX330/Addgene #42230; Cong et al. Science 2013; "
            "Ran et al. Nat Protoc 2013 sgRNA cloning protocol; "
            "BbsI-compatible CACC/AAAC-style oligo overhangs."
        ),
    ),
    "lenticrispr_v2": VectorConfig(
        assembly_method="annealed_oligos",
        enzyme="BsmBI",
        top_overhang="CACC",
        bottom_overhang="AAAC",
        u6_requires_5prime_g=True,
        notes="LentiCRISPR v2-style Cas9 cloning vector.",
        source=(
            "lentiCRISPR v2/Addgene #52961; Sanjana et al. Nat Methods 2014; "
            "Zhang lab lentiCRISPR target-guide cloning protocol uses BsmBI "
            "with CACC/AAAC-style guide oligos."
        ),
    ),
    "pdr274": VectorConfig(
        assembly_method="golden_gate",
        enzyme="BsaI",
        top_overhang="TAGG",
        bottom_overhang="AAAC",
        u6_requires_5prime_g=True,
        notes="Common zebrafish sgRNA cloning vector.",
        source=(
            "pDR274/Addgene #42250; Hwang et al. PLoS One 2013; "
            "Auer et al. Genome Res 2014 and Varshney et al. Nat Protoc 2016 "
            "describe BsaI cloning using TAGG/AAAC guide oligo overhangs."
        ),
    ),
    "pcrispr": VectorConfig(
        assembly_method="golden_gate",
        enzyme="BsaI",
        top_overhang="AAAC",
        bottom_overhang="AAAAC",
        u6_requires_5prime_g=False,
        backbone_resource="pCRISPR_rpsL",
        cell_strain="HME63 or MG1655 carrying pCas9",
        selection="Kan",
        notes=(
            "E. coli pCRISPR::rpsL guide-array plasmid used with a separate "
            "pCas9 plasmid carrying tracrRNA and Cas9. The local backbone "
            "sequence is Addgene plasmid #44505, pCRISPR::rpsL. Jiang et al. "
            "used this two-plasmid system with Lambda Red recombineering for "
            "efficient genome editing."
        ),
        source=(
            "Jiang et al. Nat Biotechnol 2013, doi:10.1038/nbt.2508; "
            "Addgene plasmid #44505 pCRISPR::rpsL; pCas9 carries tracrRNA/Cas9 "
            "and pCRISPR carries the spacer array. Supplementary Fig. 9 and "
            "Supplementary Table 2 describe BsaI insertion of spacers into "
            "pCRISPR using AAAC/AAAAC-style oligos."
        ),
    ),
    "pet28a": VectorConfig(
        assembly_method="golden_gate",
        enzyme="BsaI",
        top_overhang="TAGT",
        bottom_overhang="AAAC",
        u6_requires_5prime_g=False,
        backbone_resource="pET28a",
        cell_strain="Mach1",
        selection="Kan",
        notes="Bundled pET28a backbone GenBank file for local verification.",
        source=(
            "Local project GenBank file modules/crispr_tools/data/pET28a.gb. "
            "This is a construction-file test backbone, not a published CRISPR "
            "guide-cloning vector configuration."
        ),
    ),
    "pbr322": VectorConfig(
        assembly_method="annealed_oligos",
        enzyme="N/A",
        top_overhang="CACC",
        bottom_overhang="AAAC",
        u6_requires_5prime_g=False,
        backbone_resource="pBR322",
        cell_strain="DH5alpha",
        selection="Amp",
        notes="Bundled pBR322 backbone GenBank file for local verification.",
        source=(
            "Local project GenBank file modules/crispr_tools/data/pBR322.gb. "
            "This is a construction-file test backbone, not a published CRISPR "
            "guide-cloning vector configuration."
        ),
    ),
}


BACKBONE_RESOURCES = {
    "pET28a": BUNDLED_DATA_DIR / "pET28a.gb",
    "pBR322": BUNDLED_DATA_DIR / "pBR322.gb",
    "pCRISPR_rpsL": BUNDLED_DATA_DIR / "pCRISPR_rpsL.gbk",
}


LOCAL_REFERENCES = {
    "ecoli_rpsl": LocalReference(
        organism="Escherichia coli str. K-12 substr. MG1655",
        resource_name="ecoli_rpsl",
        path=(
            BUNDLED_DATA_DIR
            / "references"
            / "ecoli_rpsl"
            / "ncbi_dataset"
            / "data"
            / "gene.fna"
        ),
        description="Downloaded NCBI rpsL gene FASTA for E. coli K-12 MG1655.",
    ),
}


ORGANISM_ALIASES = {
    "e coli": "ecoli_rpsl",
    "e. coli": "ecoli_rpsl",
    "e_coli": "ecoli_rpsl",
    "ecoli": "ecoli_rpsl",
    "escherichia coli": "ecoli_rpsl",
    "rpsl": "ecoli_rpsl",
    "stra": "ecoli_rpsl",
}


DEFAULT_VECTOR_BY_REFERENCE = {
    "ecoli_rpsl": "pcrispr",
}


class DesignCloningOligos:
    """
    Design top and bottom strand oligos for cloning a CRISPR protospacer into a
    restriction-digested vector. The returned construction_file_inputs dict is
    intentionally shaped so it can be passed to create_construction_file.
    """

    def initiate(self) -> None:
        self._complement = {"A": "T", "T": "A", "G": "C", "C": "G"}

    def _reverse_complement(self, seq: str) -> str:
        return "".join(self._complement[b] for b in reversed(seq))

    def _validate_dna(self, seq: str, label: str) -> str:
        seq = seq.upper().strip()

        if not seq:
            raise ValueError(f"{label} must not be empty.")

        invalid = sorted(set(seq) - set("ATGC"))
        if invalid:
            raise ValueError(
                f"Invalid base(s) in {label}: {invalid}. "
                "Only A, T, G, C are accepted."
            )

        return seq

    def _read_sequence_file(self, path: Path) -> str:
        if not path.exists():
            raise ValueError(f"Required local sequence file is missing: {path}")

        text = path.read_text()
        if path.suffix.lower() in {".fa", ".fasta", ".fna"}:
            sequence = "".join(
                line.strip()
                for line in text.splitlines()
                if line.strip() and not line.startswith(">")
            )
            return self._validate_dna(sequence, path.name)

        if path.suffix.lower() in {".gb", ".gbk", ".genbank"}:
            in_origin = False
            sequence_parts = []
            for line in text.splitlines():
                if line.startswith("ORIGIN"):
                    in_origin = True
                    continue
                if in_origin and line.startswith("//"):
                    break
                if in_origin:
                    sequence_parts.extend(char for char in line if char.isalpha())
            return self._validate_dna("".join(sequence_parts), path.name)

        raise ValueError(f"Unsupported local sequence file type: {path}")

    def _resolve_local_reference(
        self,
        organism: Optional[str],
        target_reference: Optional[str],
    ) -> Optional[LocalReference]:
        if target_reference:
            key = target_reference.strip().lower().replace(" ", "_")
            key = ORGANISM_ALIASES.get(key, key)
        elif organism:
            organism_key = organism.strip().lower().replace("-", " ")
            key = ORGANISM_ALIASES.get(organism_key)
        else:
            return None

        if not key or key not in LOCAL_REFERENCES:
            valid = ", ".join(sorted(LOCAL_REFERENCES))
            raise ValueError(
                "This tool only verifies targets against downloaded local "
                f"references. Use one of: {valid}."
            )

        reference = LOCAL_REFERENCES[key]
        if not reference.path.exists():
            raise ValueError(
                f"Local reference '{key}' is configured but missing at {reference.path}."
            )
        return reference

    def _verify_protospacer_in_reference(
        self,
        protospacer: str,
        reference: Optional[LocalReference],
    ) -> dict:
        if reference is None:
            return {
                "verified": False,
                "note": (
                    "No local organism/reference was supplied, so the protospacer "
                    "was not checked against a downloaded sequence."
                ),
            }

        reference_sequence = self._read_sequence_file(reference.path)
        reverse_complement = self._reverse_complement(protospacer)
        forward_index = reference_sequence.find(protospacer)
        reverse_index = reference_sequence.find(reverse_complement)

        if forward_index == -1 and reverse_index == -1:
            raise ValueError(
                f"Protospacer was not found in local reference '{reference.resource_name}'. "
                "Choose a protospacer from the downloaded target sequence before "
                "generating cloning oligos."
            )

        if forward_index != -1:
            strand = "+"
            position = forward_index
        else:
            strand = "-"
            position = reverse_index

        return {
            "verified": True,
            "organism": reference.organism,
            "reference": reference.resource_name,
            "reference_file": str(reference.path),
            "description": reference.description,
            "strand": strand,
            "zero_based_position": position,
        }

    def _resolve_vector(
        self,
        vector: Optional[str],
        top_overhang: Optional[str],
        bottom_overhang: Optional[str],
        prepend_g: Optional[bool],
    ):
        if vector is not None:
            key = vector.lower().strip()
            if key not in VECTOR_CONFIG:
                valid = ", ".join(sorted(VECTOR_CONFIG))
                raise ValueError(
                    f"Unknown vector '{vector}'. Known vectors: {valid}."
                )
            cfg = VECTOR_CONFIG[key]

            resolved_top = cfg.top_overhang if top_overhang is None else top_overhang
            resolved_bottom = (
                cfg.bottom_overhang if bottom_overhang is None else bottom_overhang
            )
            resolved_prepend_g = (
                cfg.u6_requires_5prime_g if prepend_g is None else prepend_g
            )
            vector_notes = f"Vector '{key}' uses {cfg.enzyme} cloning. {cfg.notes}"
            return key, cfg, resolved_top, resolved_bottom, resolved_prepend_g, vector_notes

        cfg = VectorConfig(
            assembly_method="annealed_oligos",
            enzyme="N/A",
            top_overhang="CACC",
            bottom_overhang="AAAC",
            u6_requires_5prime_g=False,
            notes="No vector specified. Using default overhangs CACC / AAAC.",
        )
        resolved_top = "CACC" if top_overhang is None else top_overhang
        resolved_bottom = "AAAC" if bottom_overhang is None else bottom_overhang
        resolved_prepend_g = False if prepend_g is None else prepend_g
        return None, cfg, resolved_top, resolved_bottom, resolved_prepend_g, cfg.notes

    def _default_vector_for_reference(
        self,
        reference: Optional[LocalReference],
    ) -> Optional[str]:
        if reference is None:
            return None
        return DEFAULT_VECTOR_BY_REFERENCE.get(reference.resource_name)

    def _resolve_backbone(self, vector_key: Optional[str], cfg: VectorConfig) -> tuple[str, str]:
        backbone_resource = cfg.backbone_resource
        if backbone_resource is None and vector_key in BACKBONE_RESOURCES:
            backbone_resource = vector_key

        if backbone_resource is None:
            return vector_key or "guide_backbone", "N"

        if backbone_resource not in BACKBONE_RESOURCES:
            raise ValueError(f"No local backbone resource configured for {backbone_resource}.")

        return backbone_resource, self._read_sequence_file(BACKBONE_RESOURCES[backbone_resource])

    def run(
        self,
        protospacer: str,
        vector: Optional[str] = None,
        top_overhang: Optional[str] = None,
        bottom_overhang: Optional[str] = None,
        prepend_g: Optional[bool] = None,
        organism: Optional[str] = None,
        target_reference: Optional[str] = None,
        construct_name: Optional[str] = None,
    ) -> dict:
        """
        Design top and bottom cloning oligos.

        If organism or target_reference is supplied, the protospacer must be
        present in a downloaded local sequence. This keeps MCP runs anchored to
        organisms/data that the project can actually verify.
        """
        protospacer = self._validate_dna(protospacer, "protospacer")
        local_reference = self._resolve_local_reference(organism, target_reference)
        target_verification = self._verify_protospacer_in_reference(
            protospacer=protospacer,
            reference=local_reference,
        )
        if vector is None:
            vector = self._default_vector_for_reference(local_reference)

        (
            vector_key,
            vector_cfg,
            top_overhang,
            bottom_overhang,
            prepend_g,
            vector_notes,
        ) = self._resolve_vector(
            vector=vector,
            top_overhang=top_overhang,
            bottom_overhang=bottom_overhang,
            prepend_g=prepend_g,
        )

        top_overhang = self._validate_dna(top_overhang, "top_overhang")
        bottom_overhang = self._validate_dna(bottom_overhang, "bottom_overhang")

        g_prepended = False
        original_protospacer = protospacer

        if prepend_g and protospacer[0] != "G":
            protospacer = "G" + protospacer
            g_prepended = True

        rc = self._reverse_complement(protospacer)
        top_oligo = top_overhang + protospacer
        bottom_oligo = bottom_overhang + rc

        guide_label = vector_key or "guide"
        top_oligo_name = f"{guide_label}_top_oligo"
        bottom_oligo_name = f"{guide_label}_bottom_oligo"
        backbone_name, backbone_sequence = self._resolve_backbone(vector_key, vector_cfg)
        insert_name = f"{guide_label}_annealed_guide_insert"
        resolved_construct_name = construct_name or f"{guide_label}_guide_construct"

        notes = [vector_notes]
        if target_verification["verified"]:
            notes.append(
                "Target verified against local reference "
                f"{target_verification['reference']} "
                f"({target_verification['organism']}, strand "
                f"{target_verification['strand']}, zero-based position "
                f"{target_verification['zero_based_position']})."
            )
        else:
            notes.append(target_verification["note"])

        if g_prepended:
            notes.append(
                "A 5' G was prepended for promoter compatibility, so the cloned "
                "guide sequence is one base longer than the original protospacer."
            )
        else:
            notes.append("No 5' G was added.")

        if protospacer != original_protospacer:
            notes.append(f"Original protospacer: {original_protospacer}")
            notes.append(f"Modified protospacer: {protospacer}")

        construction_file_inputs = {
            "construct_name": resolved_construct_name,
            "assembly_strategy": "DirectSynthesis",
            "backbone_name": backbone_name,
            "backbone_sequence": backbone_sequence,
            "insert_name": insert_name,
            "insert_sequence": protospacer,
            "insert_forward_primer_name": top_oligo_name,
            "insert_forward_primer_sequence": top_oligo,
            "insert_reverse_primer_name": bottom_oligo_name,
            "insert_reverse_primer_sequence": bottom_oligo,
            "vector_forward_primer_name": "",
            "vector_forward_primer_sequence": "",
            "vector_reverse_primer_name": "",
            "vector_reverse_primer_sequence": "",
            "enzyme": "" if vector_cfg.enzyme == "N/A" else vector_cfg.enzyme,
            "cell_strain": vector_cfg.cell_strain,
            "selection": vector_cfg.selection,
            "temperature_c": 37,
            "notes": " ".join(notes),
        }

        return {
            "vector": vector_key,
            "enzyme": vector_cfg.enzyme,
            "source": vector_cfg.source,
            "assembly_method": vector_cfg.assembly_method,
            "top_overhang": top_overhang,
            "bottom_overhang": bottom_overhang,
            "top_oligo_name": top_oligo_name,
            "bottom_oligo_name": bottom_oligo_name,
            "top_oligo": top_oligo,
            "bottom_oligo": bottom_oligo,
            "g_prepended": g_prepended,
            "final_protospacer": protospacer,
            "target_verification": target_verification,
            "construction_file_inputs": construction_file_inputs,
        }


_instance = DesignCloningOligos()
_instance.initiate()
design_cloning_oligos = _instance.run
