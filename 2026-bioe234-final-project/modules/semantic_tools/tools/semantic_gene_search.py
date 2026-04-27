"""
Semantic Gene Search MCP tool.

Natural language biology query -> GO / ontology terms

This tool only performs semantic concept discovery.
Gene lookup should be handled by go_term_gene_lookup.
"""

from modules.semantic_tools.semantic_wrapper import SemanticGeneWrapper


class SemanticGeneSearch:
    """
    Description:
        Search for biology concepts using a natural-language query by mapping
        the query to Gene Ontology / ontology terms.

    Input:
        query (str): Natural-language biology query, e.g.
                     "oxidative stress in yeast"

    Output:
        dict: Structured result containing:
              - parsed_query
              - go_terms
    """

    def initiate(self) -> None:
        self.wrapper = SemanticGeneWrapper()

    def run(
        self,
        query: str,
    ) -> dict:
        result = self.wrapper.run(query=query)
        return result.to_dict()


_instance = SemanticGeneSearch()
_instance.initiate()
semantic_gene_search = _instance.run