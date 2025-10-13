# core/state.py
import numpy as np
import torch
from typing import Dict, Any, Optional, List, Tuple, Sequence

class SystemState:
    """Lightweight container for system state that can be efficiently passed between samplers"""

    def __init__(
        self,
        sampler_sequence: Optional[Sequence[str]] = None,
        current_sampler: Optional[str] = None
    ):
        # Core data
        self.positions: Dict[str, np.ndarray] = {}
        self.sigma_range: Dict[str, Tuple[float, float]] = {}
        self.sigma_range = {'AA': (0.2, 12.0), 'AB': (0.1, 10.0), 'BC': (0.15, 11.0)}  # Default sigma ranges
        self.sigma: Dict[str, float] = {}
        # choose a random number within the range as default
        self.sigma = {k: np.random.uniform(v[0], v[1]) for k, v in self.sigma_range.items()}
        self.box_size: float = 0.0

        # Sampler sequencing
        # User-defined ordered list of samplers (e.g. ["pair","tetramer","octet", "full"])
        self.sampler_sequence: List[str] = list(sampler_sequence) if sampler_sequence is not None else []
        # Name of the sampler currently operating on this state
        self.current_sampler: Optional[str] = current_sampler

        # Cached structures (clear when positions or sampler change)
        self._tetramers: Optional[List[Tuple]] = None
        self._octets: Optional[List[Tuple]] = None

        # Flags
        self.use_sigma_distribution: bool = False

        # Additional metadata
        self.metadata: Dict[str, Any] = {}
        
    def __getstate__(self):
        """Custom pickle support - explicitly save all attributes"""
        state = self.__dict__.copy()
        # Explicitly include sigma and sigma_range
        return state
    
    def __setstate__(self, state):
        """Custom unpickle support - restore all attributes"""
        self.__dict__.update(state)
        # Ensure sigma exists
        if not hasattr(self, 'sigma') or self.sigma is None:
            self.sigma = {}
        if not hasattr(self, 'sigma_range') or self.sigma_range is None:
            self.sigma_range = {}
            
    def update_positions(self, new_positions: Dict[str, np.ndarray]) -> None:
        """Update positions and clear cached structures"""
        self.positions = {k: v.copy() for k, v in new_positions.items()}
        # Clear caches when positions change
        self._tetramers = None
        self._octets = None

    def update_sigma(self, new_sigma: Dict[str, float]) -> None:
        """Update sigma values"""
        self.sigma = {k: v for k, v in new_sigma.items()}

    def copy(self) -> 'SystemState':
        """Create a minimal copy with only the essential state"""
        state_copy = SystemState(self.sampler_sequence, self.current_sampler)
        state_copy.positions = {k: v.copy() for k, v in self.positions.items()}
        state_copy.sigma = {k: v for k, v in self.sigma.items()}
        state_copy.sigma_range = {k: v for k, v in self.sigma_range.items()}
        state_copy.box_size = self.box_size
        state_copy.use_sigma_distribution = self.use_sigma_distribution
        state_copy.metadata = dict(self.metadata)
        
        # Copy sigma_prior if it exists
        if hasattr(self, 'sigma_prior'):
            state_copy.sigma_prior = self.sigma_prior
        
        return state_copy

    @property
    def tetramers(self) -> List[Tuple]:
        """Get tetramers, calculating if needed and if sampler sequence requires it"""
        # Only compute tetramers if "tetramer" appears before or at current sampler
        if self._tetramers is None:
            if self._should_compute('tetramer'):
                from samplers.tetramer import get_tetramers
                self._tetramers = get_tetramers(self.positions)
            else:
                self._tetramers = []
        return self._tetramers

    @property
    def octets(self) -> List[Tuple]:
        """Always compute tetramers and octets if octet sampler is active"""
        if self._octets is None:
            from samplers.octet import get_octets

            # Always recompute tetramers when octets are needed
            self._octets, self._tetramers = get_octets(self.positions)
        return self._octets


    def _should_compute(self, target: str) -> bool:
        """
        Determine whether to compute a derived structure (e.g. tetramers or octets)
        based on the user-defined sampler sequence and the current sampler.
        """
        if not self.sampler_sequence or self.current_sampler is None:
            # no sequence defined: default to always compute
            return True
        try:
            target_idx = self.sampler_sequence.index(target)
            current_idx = self.sampler_sequence.index(self.current_sampler)
        except ValueError:
            # target or current sampler not in sequence: do not compute
            return False
        # compute if target stage is at or before current stage
        return target_idx <= current_idx