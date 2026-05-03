from modules.locus_tools.gene_locus_lookup import GeneLocusLookup


def test_gene_locus_lookup_resolves_yap1_and_returns_coordinates():
    lookup = GeneLocusLookup()
    result = lookup.run(
        gene_symbol="YAP1",
        organism="Saccharomyces cerevisiae",
        max_loci=1,
        include_fasta=False,
    ).to_dict()

    assert result["resolved_symbol"] == "YAP1"
    assert len(result["loci"]) == 1

    locus = result["loci"][0]
    assert locus["chr_accession"] is not None
    assert locus["start_1_based"] < locus["stop_1_based"]
    assert locus["strand"] in {"plus", "minus"}


def test_gene_locus_lookup_can_fetch_fasta():
    lookup = GeneLocusLookup()
    result = lookup.run(
        gene_symbol="YAP1",
        organism="Saccharomyces cerevisiae",
        max_loci=1,
        include_fasta=True,
    ).to_dict()

    locus = result["loci"][0]
    fasta = locus["fasta"]

    assert fasta is not None
    assert fasta.startswith(">")
    assert "NC_" in fasta or "BK" in fasta