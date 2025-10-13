#============================================================================
# Tetrameric scoring function
#============================================================================
import numpy as np
from core.parameters import SystemParameters
from core.movers import get_tetramers

class TetramerNLL:
    def __init__(self, state):
        """
        Initialize tetramer scoring with system state.
        
        Args:
            state: SystemState object containing positions and sigma values
        """
        self.state = state
        self.coordinates = state.positions
        self.nuisance_parameters = state.sigma
        self.params = SystemParameters()

    def calculate_tetramer_scores_batch(self, positions, tetramers, sig, debug_logging=False):
        """Calculate scores for all tetramers using vectorized operations."""
        if not tetramers:
            return np.array([], dtype=np.float32)
        
        # Extract indices
        a_indices = np.array([t[0] for t in tetramers], dtype=np.int32)
        b_indices = np.array([t[1] for t in tetramers], dtype=np.int32)
        c1_indices = np.array([t[2] for t in tetramers], dtype=np.int32)
        c2_indices = np.array([t[3] for t in tetramers], dtype=np.int32)
        
        # Get positions
        pos_a = positions['A'][a_indices]
        pos_b = positions['B'][b_indices]
        pos_c1 = positions['C'][c1_indices]
        pos_c2 = positions['C'][c2_indices]
        
        # Calculate distances
        ab_dists = np.linalg.norm(pos_a - pos_b, axis=1)
        bc1_dists = np.linalg.norm(pos_b - pos_c1, axis=1)
        bc2_dists = np.linalg.norm(pos_b - pos_c2, axis=1)
        
        # Target distances
        ab_target = self.params.pair_distances['AB']
        bc_target = self.params.pair_distances['BC']
        
        # Calculate negative log-likelihood scores
        # NLL = (d - d0)²/(2σ²) + 0.5 * log(2πσ²)
        ab_scores = ((ab_dists - ab_target)**2) / (2 * sig['AB']**2) + 0.5 * np.log(2 * np.pi * sig['AB']**2)
        bc1_scores = ((bc1_dists - bc_target)**2) / (2 * sig['BC']**2) + 0.5 * np.log(2 * np.pi * sig['BC']**2)
        bc2_scores = ((bc2_dists - bc_target)**2) / (2 * sig['BC']**2) + 0.5 * np.log(2 * np.pi * sig['BC']**2)
        
        # Total score per tetramer (AB + BC1 + BC2)
        scores = ab_scores + bc1_scores + bc2_scores
        
        if debug_logging:
            print(f"Tetramer batch scores: mean={np.mean(scores):.3f}, std={np.std(scores):.3f}")
        
        return scores

    def compute_score(self) -> float:
        """
        Tetrameric scoring function. A tetramer is ABCC, so 
        score will be sum of pair scores for AB, BC1 and BC2.
        This has to be added to the pair score for the whole system to 
        obtain the total negative log likelihood and then adding the 
        prior score to get the posterior, which is done in the tetramer_sampler
        code.
        
        Returns:
            float: Total negative log-likelihood score for all tetramers
        """
        # Get tetramers from the current state
        tetramers = get_tetramers(self.state)
        
        if not tetramers:
            return 0.0
        
        # Calculate scores for all tetramers using batch method
        tetramer_scores = self.calculate_tetramer_scores_batch(
            positions=self.coordinates,
            tetramers=tetramers,
            sig=self.nuisance_parameters,
            debug_logging=False
        )
        
        # Sum all tetramer scores
        total_score = np.sum(tetramer_scores)
        
        return float(total_score)
#============================================================================