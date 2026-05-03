from modules.sequence_tools.gene_sequence_lookup import GeneSequenceLookup


def test_gene_sequence_lookup_resolves_yap1_correctly():
    lookup = GeneSequenceLookup()
    result = lookup.run(
        gene_symbol="YAP1",
        organism="Saccharomyces cerevisiae",
        max_nucleotide_records=2,
        include_fasta=False,
    ).to_dict()

    assert result["resolved_symbol"] == "YAP1"
    assert result["resolved_gene_id"] is not None
    assert isinstance(result["nucleotide_records"], list)
    assert len(result["nucleotide_records"]) >= 1