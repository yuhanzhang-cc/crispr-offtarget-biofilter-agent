"""Gene importance database adapter."""

from typing import Optional, Dict


class GeneImportanceDB:
    """
    Local gene importance database.
    
    This is a placeholder implementation using a simple dictionary.
    In production, this would connect to external databases like:
    - OMIM (morbid genes)
    - Cancer Gene Census
    - Essential gene databases
    """
    
    # Mock data: gene -> importance score (0-1)
    # Higher = more biologically important/concerning if edited
    _IMPORTANCE_DATA: Dict[str, float] = {
        # Tumor suppressors (high importance)
        "TP53": 1.0,
        "PTEN": 0.95,
        "BRCA1": 0.95,
        "BRCA2": 0.95,
        "RB1": 0.9,
        "APC": 0.9,
        
        # Oncogenes (high importance)
        "MYC": 0.95,
        "KRAS": 0.95,
        "EGFR": 0.9,
        "BRAF": 0.9,
        "PIK3CA": 0.85,
        
        # Signaling pathway genes
        "AKT1": 0.8,
        "MTOR": 0.8,
        "PTK2": 0.7,
        "SRC": 0.75,
        
        # Housekeeping genes (moderate importance)
        "ACTB": 0.5,
        "GAPDH": 0.5,
        "TUBB": 0.5,
        
        # Developmental regulators
        "SOX2": 0.85,
        "OCT4": 0.85,
        "NANOG": 0.85,
    }
    
    def __init__(self):
        """Initialize database connection."""
        self._cache = {}
    
    def get_importance(self, gene_name: str) -> float:
        """
        Get importance score for a gene.
        
        Args:
            gene_name: HGNC gene symbol
            
        Returns:
            Importance score (0-1), default 0.5 if unknown
        """
        # Normalize gene name
        gene_name = gene_name.upper().strip()
        
        # Check cache first
        if gene_name in self._cache:
            return self._cache[gene_name]
        
        # Look up in database
        score = self._IMPORTANCE_DATA.get(gene_name, 0.5)
        
        # Cache result
        self._cache[gene_name] = score
        
        return score
    
    def is_essential(self, gene_name: str, threshold: float = 0.8) -> bool:
        """Check if gene is considered essential/important."""
        return self.get_importance(gene_name) >= threshold
    
    def get_all_genes(self) -> list:
        """Return list of all genes in database."""
        return list(self._IMPORTANCE_DATA.keys())
