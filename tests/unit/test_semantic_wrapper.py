from modules.semantic_tools.semantic_wrapper import SemanticGeneWrapper


def test_parse_oxidative_stress_query_returns_expected_go_term():
    wrapper = SemanticGeneWrapper()
    result = wrapper.run("oxidative stress in yeast").to_dict()

    assert result["parsed_query"]["organism"] == "Saccharomyces cerevisiae"
    assert "response to oxidative stress" in result["parsed_query"]["ontology_terms"]

    go_ids = [term["go_id"] for term in result["go_terms"]]
    assert "GO:0006979" in go_ids