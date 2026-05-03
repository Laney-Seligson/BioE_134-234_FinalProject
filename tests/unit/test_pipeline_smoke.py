from modules.semantic_tools.semantic_wrapper import SemanticGeneWrapper
from modules.annotation_tools.go_term_to_genes import GOTermGeneLookup
from modules.locus_tools.gene_locus_lookup import GeneLocusLookup


def test_oxidative_stress_pipeline_smoke():
    semantic = SemanticGeneWrapper()
    sem_result = semantic.run("oxidative stress in yeast").to_dict()

    go_terms = sem_result["go_terms"]
    assert len(go_terms) > 0

    top_term = go_terms[0]
    assert top_term["go_id"] == "GO:0006979"

    annot = GOTermGeneLookup()
    annot_result = annot.run(
        go_id=top_term["go_id"],
        go_label=top_term["label"],
        organism="Saccharomyces cerevisiae",
        max_genes=15,
    ).to_dict()

    assert len(annot_result["genes"]) > 0

    symbols = {g["symbol"] for g in annot_result["genes"]}
    assert "YAP1" in symbols or "SOD1" in symbols

    locus = GeneLocusLookup()
    locus_result = locus.run(
        gene_symbol="YAP1",
        organism="Saccharomyces cerevisiae",
        max_loci=1,
        include_fasta=False,
    ).to_dict()

    assert locus_result["resolved_symbol"] == "YAP1"
    assert len(locus_result["loci"]) == 1