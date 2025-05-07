# core/state.py
import numpy as np
import torch
from typing import Dict, Any, Optional, List, Tuple

class SystemState:
    """Lightweight container for system state that can be efficiently passed between samplers"""
    
    def __init__(self):
        # Core data
        self.positions: Dict[str, np.ndarray] = {}
        self.sigma: Dict[str, float] = {}
        self.sigma_range: Dict[str, Tuple[float, float]] = {}
        self.box_size: float = 0.0
        
        # Cached structures (clear when positions change)
        self._tetramers: Optional[List[Tuple]] = None
        self._octets: Optional[List[Tuple]] = None
        
        # Flags
        self.use_sigma_distribution: bool = False
        
        # Additional metadata
        self.metadata: Dict[str, Any] = {}
        
    def update_positions(self, new_positions: Dict[str, np.ndarray]) -> None:
        """Update positions and clear cached structures"""
        self.positions = {k: v.copy() for k, v in new_positions.items()}
        self._tetramers = None
        self._octets = None
        
    def update_sigma(self, new_sigma: Dict[str, float]) -> None:
        """Update sigma values"""
        self.sigma = {k: v for k, v in new_sigma.items()}
    
    def copy(self) -> 'SystemState':
        """Create a minimal copy with only the essential state"""
        state_copy = SystemState()
        state_copy.positions = {k: v.copy() for k, v in self.positions.items()}
        state_copy.sigma = {k: v for k, v in self.sigma.items()}
        state_copy.sigma_range = {k: v for k, v in self.sigma_range.items()}
        state_copy.box_size = self.box_size
        state_copy.use_sigma_distribution = self.use_sigma_distribution
        return state_copy
    
    @property
    def tetramers(self) -> List[Tuple]:
        """Get tetramers, calculating if needed"""
        if self._tetramers is None:
            from samplers.tetramer import get_tetramers
            self._tetramers = get_tetramers(self.positions)
        return self._tetramers
    
    @property
    def octets(self) -> List[Tuple]:
        """Get octets, calculating if needed"""
        if self._octets is None:
            from samplers.octet import get_octets
            self._octets = get_octets(self.positions, self.tetramers)
        return self._octets