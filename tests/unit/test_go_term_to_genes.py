from modules.annotation_tools.go_term_to_genes import GOTermGeneLookup


def test_go_term_gene_lookup_returns_expected_shape():
    lookup = GOTermGeneLookup()
    result = lookup.run(
        go_id="GO:0006979",
        go_label="response to oxidative stress",
        organism="Saccharomyces cerevisiae",
        max_genes=10,
    ).to_dict()

    assert result["go_id"] == "GO:0006979"
    assert result["organism"] == "Saccharomyces cerevisiae"
    assert isinstance(result["genes"], list)
    assert len(result["genes"]) > 0


def test_go_term_gene_lookup_contains_plausible_oxidative_stress_genes():
    lookup = GOTermGeneLookup()
    result = lookup.run(
        go_id="GO:0006979",
        go_label="response to oxidative stress",
        organism="Saccharomyces cerevisiae",
        max_genes=15,
    ).to_dict()

    symbols = {g["symbol"] for g in result["genes"]}
    assert "YAP1" in symbols or "SOD1" in symbols