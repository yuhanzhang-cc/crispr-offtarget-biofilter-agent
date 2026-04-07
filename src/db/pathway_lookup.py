"""Pathway importance database adapter."""

from typing import Dict, Set


class PathwayImportanceDB:
    """
    Local pathway importance database.
    
    Placeholder implementation. In production, would connect to:
    - KEGG pathways
    - Reactome
    - Gene Ontology
    - Disease pathway databases
    """
    
    # Mock data: pathway -> set of genes
    _PATHWAYS: Dict[str, Set[str]] = {
        "p53_signaling": {"TP53", "MDM2", "CDKN1A", "BAX"},
        "PI3K_AKT_mTOR": {"PIK3CA", "AKT1", "MTOR", "PTEN", "TSC1", "TSC2"},
        "RAS_RAF_MEK_ERK": {"KRAS", "BRAF", "MAPK1", "MAPK3", "EGFR"},
        "cell_cycle": {"CDK4", "CDK6", "CCND1", "RB1", "TP53"},
        "DNA_repair": {"BRCA1", "BRCA2", "ATM", "ATR", "CHEK2"},
        "apoptosis": {"BCL2", "BAX", "BAK1", "CASP9", "CASP3"},
        "stem_cell_pluripotency": {"SOX2", "OCT4", "NANOG", "MYC"},
    }
    
    # Pathway importance scores
    _PATHWAY_IMPORTANCE: Dict[str, float] = {
        "p53_signaling": 1.0,
        "PI3K_AKT_mTOR": 0.95,
        "RAS_RAF_MEK_ERK": 0.95,
        "cell_cycle": 0.9,
        "DNA_repair": 0.9,
        "apoptosis": 0.85,
        "stem_cell_pluripotency": 0.9,
    }
    
    def __init__(self):
        """Initialize database."""
        self._gene_to_pathways: Dict[str, Set[str]] = {}
        self._build_gene_index()
    
    def _build_gene_index(self):
        """Build reverse index from gene to pathways."""
        for pathway, genes in self._PATHWAYS.items():
            for gene in genes:
                gene = gene.upper()
                if gene not in self._gene_to_pathways:
                    self._gene_to_pathways[gene] = set()
                self._gene_to_pathways[gene].add(pathway)
    
    def get_pathways_for_gene(self, gene_name: str) -> Set[str]:
        """
        Get all pathways a gene belongs to.
        
        Args:
            gene_name: HGNC gene symbol
            
        Returns:
            Set of pathway names
        """
        gene_name = gene_name.upper().strip()
        return self._gene_to_pathways.get(gene_name, set())
    
    def get_importance(self, gene_name: str) -> float:
        """
        Get pathway-based importance for a gene.
        
        Returns max importance across all pathways the gene belongs to.
        Returns 0.0 if gene is not in any known pathway.
        """
        gene_name = gene_name.upper().strip()
        pathways = self.get_pathways_for_gene(gene_name)
        
        if not pathways:
            return 0.0
        
        # Return max importance across all pathways
        return max(
            self._PATHWAY_IMPORTANCE.get(p, 0.5)
            for p in pathways
        )
    
    def get_pathway_details(self, gene_name: str) -> Dict[str, float]:
        """
        Get detailed pathway information for a gene.
        
        Returns:
            Dict mapping pathway name to importance score
        """
        gene_name = gene_name.upper().strip()
        pathways = self.get_pathways_for_gene(gene_name)
        
        return {
            p: self._PATHWAY_IMPORTANCE.get(p, 0.5)
            for p in pathways
        }
    
    def is_in_pathway(self, gene_name: str, pathway: str) -> bool:
        """Check if gene is in specific pathway."""
        gene_name = gene_name.upper().strip()
        pathway = pathway.lower()
        pathways = self.get_pathways_for_gene(gene_name)
        return any(p.lower() == pathway for p in pathways)
