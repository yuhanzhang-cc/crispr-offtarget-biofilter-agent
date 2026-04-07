"""LLM interpreter for generating natural language reports."""

import os
from typing import List, Optional
from dataclasses import dataclass

from src.models.schema import ScoredSite, RiskLevel


@dataclass
class LLMConfig:
    """Configuration for LLM client."""
    api_key: Optional[str] = None
    base_url: str = "https://api.openai.com/v1"
    model: str = "gpt-4"
    temperature: float = 0.3
    max_tokens: int = 1000


class MockLLMClient:
    """Mock LLM client for testing without API access."""
    
    def generate(self, prompt: str) -> str:
        """Generate a mock response."""
        return f"""[MOCK LLM RESPONSE]

Based on the provided off-target analysis data, here is a summary:

{prompt[:500]}...

[This is a mock response. Configure OPENAI_API_KEY for real LLM generation.]
"""


class OpenAIClient:
    """Real OpenAI API client."""
    
    def __init__(self, config: LLMConfig):
        self.config = config
        self._client = None
    
    def _get_client(self):
        """Lazy import and initialize client."""
        if self._client is None:
            try:
                from openai import OpenAI
            except ImportError:
                raise ImportError("openai package required. Install with: pip install openai")
            
            self._client = OpenAI(
                api_key=self.config.api_key,
                base_url=self.config.base_url,
            )
        return self._client
    
    def generate(self, prompt: str) -> str:
        """Generate text using OpenAI API."""
        client = self._get_client()
        
        try:
            response = client.chat.completions.create(
                model=self.config.model,
                messages=[
                    {"role": "system", "content": "You are a bioinformatics expert analyzing CRISPR off-target risk."},
                    {"role": "user", "content": prompt},
                ],
                temperature=self.config.temperature,
                max_tokens=self.config.max_tokens,
            )
            return response.choices[0].message.content
        except Exception as e:
            return f"[LLM Error: {e}]"


class LLMInterpreter:
    """Generates natural language interpretations of risk scores."""
    
    def __init__(self, config: Optional[LLMConfig] = None):
        """
        Initialize interpreter.
        
        Args:
            config: LLM configuration. If None, uses environment variables.
        """
        self.config = config or self._load_config_from_env()
        self.client = self._create_client()
    
    def _load_config_from_env(self) -> LLMConfig:
        """Load configuration from environment variables."""
        return LLMConfig(
            api_key=os.getenv("OPENAI_API_KEY"),
            base_url=os.getenv("OPENAI_BASE_URL", "https://api.openai.com/v1"),
            model=os.getenv("OPENAI_MODEL", "gpt-4"),
        )
    
    def _create_client(self):
        """Create appropriate LLM client."""
        if self.config.api_key:
            return OpenAIClient(self.config)
        else:
            return MockLLMClient()
    
    def interpret(
        self,
        sgRNA: str,
        top_hits: List[ScoredSite],
        max_sites: int = 10,
    ) -> str:
        """
        Generate natural language summary of top risk sites.
        
        Args:
            sgRNA: The sgRNA sequence
            top_hits: List of scored sites (should be pre-sorted by score)
            max_sites: Maximum number of sites to include in interpretation
            
        Returns:
            Natural language summary
        """
        # Select top sites
        sites_to_interpret = top_hits[:max_sites]
        
        # Build prompt
        prompt = self._build_prompt(sgRNA, sites_to_interpret)
        
        # Generate interpretation
        return self.client.generate(prompt)
    
    def _build_prompt(self, sgRNA: str, sites: List[ScoredSite]) -> str:
        """Build prompt for LLM."""
        lines = [
            f"CRISPR Off-Target Risk Analysis for sgRNA: {sgRNA}",
            "",
            f"Total high-risk sites identified: {len(sites)}",
            "",
            "Top off-target sites by risk score:",
            "",
        ]
        
        for i, site in enumerate(sites, 1):
            lines.append(f"{i}. {site.name}")
            lines.append(f"   Location: {site.chrom}:{site.start}-{site.end} ({site.strand})")
            lines.append(f"   Region: {site.region_type.value}")
            lines.append(f"   Gene: {site.gene_name or 'None'}")
            lines.append(f"   Mismatches: {site.mismatch_count}")
            lines.append(f"   Risk Score: {site.total_score:.3f} ({site.risk_level.value.upper()})")
            lines.append(f"   Score Breakdown:")
            lines.append(f"     - Region component: {site.components.region_score:.3f}")
            lines.append(f"     - Gene importance: {site.components.gene_score:.3f}")
            lines.append(f"     - Pathway involvement: {site.components.pathway_score:.3f}")
            lines.append("")
        
        lines.extend([
            "Please provide:",
            "1. A concise summary of the overall off-target risk profile",
            "2. Key concerns based on the high-risk sites identified",
            "3. Recommendations for experimental validation",
        ])
        
        return "\n".join(lines)
    
    def generate_summary_text(
        self,
        sgRNA: str,
        total_sites: int,
        high_risk: int,
        medium_risk: int,
        low_risk: int,
        top_hits: List[ScoredSite],
    ) -> str:
        """
        Generate a structured summary without LLM (deterministic).
        
        This is useful when LLM is not available or for reproducible reports.
        """
        lines = [
            "=" * 60,
            f"CRISPR Off-Target Risk Assessment Report",
            "=" * 60,
            "",
            f"sgRNA Sequence: {sgRNA}",
            f"Total Off-Target Sites Analyzed: {total_sites}",
            "",
            "Risk Distribution:",
            f"  High Risk:   {high_risk} sites",
            f"  Medium Risk: {medium_risk} sites",
            f"  Low Risk:    {low_risk} sites",
            "",
            "-" * 60,
            "Top High-Risk Off-Target Sites:",
            "-" * 60,
            "",
        ]
        
        for i, site in enumerate(top_hits[:10], 1):
            lines.append(f"{i}. {site.name}")
            lines.append(f"   Location: {site.chrom}:{site.start}-{site.end}")
            lines.append(f"   Gene: {site.gene_name or 'N/A'} | Region: {site.region_type.value}")
            lines.append(f"   Risk Score: {site.total_score:.3f}")
            lines.append("")
        
        lines.extend([
            "-" * 60,
            "Scoring Method:",
            "  Risk = 0.4 * Region + 0.35 * Gene + 0.25 * Pathway",
            "  Region priority: exon > promoter > intron > intergenic",
            "-" * 60,
        ])
        
        return "\n".join(lines)
