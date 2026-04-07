"""Scoring utilities."""


def normalize_score(value: float, min_val: float = 0.0, max_val: float = 1.0) -> float:
    """
    Normalize a score to [0, 1] range.
    
    Args:
        value: Input value
        min_val: Minimum expected value
        max_val: Maximum expected value
        
    Returns:
        Normalized value in [0, 1]
    """
    if max_val == min_val:
        return 0.5
    
    normalized = (value - min_val) / (max_val - min_val)
    return max(0.0, min(1.0, normalized))


def clamp(value: float, min_val: float = 0.0, max_val: float = 1.0) -> float:
    """Clamp value to range."""
    return max(min_val, min(max_val, value))


def weighted_average(values: list, weights: list) -> float:
    """
    Compute weighted average.
    
    Args:
        values: List of values
        weights: List of weights (same length as values)
        
    Returns:
        Weighted average
    """
    if len(values) != len(weights):
        raise ValueError("values and weights must have same length")
    
    if not values:
        return 0.0
    
    total_weight = sum(weights)
    if total_weight == 0:
        return sum(values) / len(values)
    
    return sum(v * w for v, w in zip(values, weights)) / total_weight
