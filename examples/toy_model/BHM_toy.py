import numpy as np
import pickle
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import Circle
import mpl_toolkits.mplot3d.art3d as art3d
from scipy.spatial.distance import cdist
from dataclasses import dataclass
from typing import Dict, List, Tuple
import random

# Define the 8-fold ideal symmetric system to evaluate the likelihood
array_A = np.array([[0, 100, 0],
                   [-70.71, 70.71, 0],
                   [-100, 0, 0],
                   [-70.71, -70.71, 0],
                   [0, -100, 0],
                   [70.71, -70.71, 0],
                   [100, 0, 0],
                   [70.71, 70.71, 0]])
array_B = np.array([[0, 100, -50],
                   [-70.71, 70.71, -50],
                   [-100, 0, -50],
                   [-70.71, -70.71, -50],
                   [0, -100, -50],
                   [70.71, -70.71, -50],
                   [100, 0, -50],
                   [70.71, 70.71, -50]])
array_C = np.array([[11.25, 100, -68],    # pair 1
                    [-11.25, 100, -68],
                    [-59.46, 70.71, -68], # pair 2
                    [-81.96, 70.71, -68],
                    [-100, 11.25, -68],   # pair 3
                    [-100, -11.25, -68],
                    [-81.96, -70.71, -68], # pair 4
                    [-59.46, -70.71, -68],
                    [11.25, -100, -68],    # pair 5
                    [-11.25, -100, -68],
                    [81.96, -70.71, -68],  # pair 6
                    [59.46, -70.71, -68],
                    [100, 11.25, -68],     # pair 7
                    [100, -11.25, -68],
                    [81.96, 70.71, -68],   # pair 8
                    [59.46, 70.71, -68]])

ideal_coord = {'A': array_A,
               'B': array_B,
               'C': array_C}

class PairSampler:
    def __init__(self):
        self.box_size = 300
        self.radii = {'A': 40, 'B': 10, 'C': 16}
        self.pair_distances = {
            'AA': 80.5,
            'AB': 50.5,
            'CC': 31.5,
            'BC': 36.5
        }
        # Tighter tolerances for more focused sampling
        self.distance_tolerances = {
            'AA': 2.0,
            'AB': 2.0,
            'CC': 2.0,
            'BC': 2.0
        }
        self.counts = {'A': 8, 'B': 8, 'C': 16}
        # Minimum likelihood value before taking log
        self.min_likelihood = 1e-10

    def initialize_positions(self):
        """Initialize positions for all components"""
        positions = {}
        for type_name, count in self.counts.items():
            # Create array of shape (count, 3) for x,y,z coordinates
            positions[type_name] = np.random.uniform(0, self.box_size, (count, 3))
        return positions
        
    def calculate_pair_neg_log_likelihood_matrix(self, pos1, pos2, target_dist, tolerance):
        """Calculate negative log likelihood matrix for pairs"""
        distances = cdist(pos1, pos2)
        # Gaussian likelihood with minimum cutoff
        likelihood = np.exp(-(distances - target_dist)**2 / (2 * tolerance**2))
        likelihood = np.maximum(likelihood, self.min_likelihood)
        return -np.log(likelihood)
    
    def calculate_exclusion_neg_log_likelihood(self, positions):
        """Calculate negative log likelihood for excluded volume"""
        neg_log_likelihood = 0
        k_exclude = 100  # Strong exclusion penalty
        
        for type1 in positions:
            for type2 in positions:
                if type1 <= type2:
                    pos1, pos2 = positions[type1], positions[type2]
                    min_dist = self.radii[type1] + self.radii[type2]
                    distances = cdist(pos1, pos2)
                    
                    if type1 == type2:
                        mask = np.triu(np.ones_like(distances), k=1)
                        violations = np.where((distances < min_dist) & (mask > 0))
                    else:
                        violations = np.where(distances < min_dist)
                    
                    if violations[0].size > 0:
                        penalty = k_exclude * np.sum((min_dist - distances[violations])**2)
                        neg_log_likelihood += penalty
        
        return neg_log_likelihood
    
    def calculate_system_neg_log_likelihood(self, positions):
        """Calculate total system negative log likelihood"""
        neg_log_likelihood = 0
        
        # Add exclusion term
        neg_log_likelihood += self.calculate_exclusion_neg_log_likelihood(positions)
        
        # Calculate pair likelihoods
        pair_types = [('A', 'A'), ('A', 'B'), ('B', 'C'), ('C', 'C')]
        for type1, type2 in pair_types:
            if f"{type1}{type2}" in self.pair_distances:
                target_dist = self.pair_distances[f"{type1}{type2}"]
                tolerance = self.distance_tolerances[f"{type1}{type2}"]
                
                neg_log_likelihood_matrix = self.calculate_pair_neg_log_likelihood_matrix(
                    positions[type1], 
                    positions[type2] if type1 != type2 else positions[type2],
                    target_dist,
                    tolerance
                )
                
                if type1 == type2:
                    mask = np.triu(np.ones_like(neg_log_likelihood_matrix), k=1)
                    neg_log_likelihood_matrix *= mask
                
                # Get indices of minimum values (best pairs)
                row_mins = np.min(neg_log_likelihood_matrix, axis=1)
                col_mins = np.min(neg_log_likelihood_matrix, axis=0)
                
                # Add unique contributions
                unique_pairs = set()
                for i, j in enumerate(np.argmin(neg_log_likelihood_matrix, axis=1)):
                    unique_pairs.add((i, j))
                for i in np.argmin(neg_log_likelihood_matrix, axis=0):
                    unique_pairs.add((i, j))
                
                for i, j in unique_pairs:
                    neg_log_likelihood += neg_log_likelihood_matrix[i, j]
        
        return neg_log_likelihood
    
    def propose_move(self, positions):
        """Propose a move with adaptive step size"""
        new_positions = {k: v.copy() for k, v in positions.items()}
    
        # Select random component type and index
        type_name = np.random.choice(list(self.counts.keys()))
        idx = np.random.randint(self.counts[type_name])
    
        # Adaptive step size based on component type
        step_size = self.radii[type_name] * 0.1  # 10% of radius
    
        # Apply random displacement to one particle's position
        new_positions[type_name][idx, :] += np.random.normal(0, step_size, 3)
    
        # Ensure within box with soft boundary
        new_positions[type_name][idx, :] = np.clip(
            new_positions[type_name][idx, :], -10, self.box_size + 10)
    
        return new_positions

    def save_state(self, positions, step):
        """Save current state information"""
        state = {
            'step': step,
            'positions': positions,
            'sigmas': {
                'sigma_AA': self.distance_tolerances['AA'],
                'sigma_AB': self.distance_tolerances['AB'],
                'sigma_CC': self.distance_tolerances['CC'],
                'sigma_BC': self.distance_tolerances['BC']
            },
            'types': {i: 'A' for i in range(8)},
            'bead_numbers': {i: i+1 for i in range(32)}
        }
        # Add type information for all beads
        state['types'].update({i: 'B' for i in range(8, 16)})
        state['types'].update({i: 'C' for i in range(16, 32)})
        return state
    
    def run_mc(self, n_steps=50000, initial_temperature=10.0, save_freq=10000, traj_filename='trajectory'):
        """Run Monte Carlo sampling with temperature schedule"""
        positions = self.initialize_positions()
        current_neg_log_likelihood = self.calculate_system_neg_log_likelihood(positions)
        
        best_positions = {k: v.copy() for k, v in positions.items()}
        best_neg_log_likelihood = current_neg_log_likelihood
        
        neg_log_likelihoods = [current_neg_log_likelihood]
        
        # Temperature schedule
        temperature = initial_temperature
        cooling_rate = 0.99
        min_temperature = 0.1

        # Initialize trajectory
        self.trajectory = []
        
        accepted = 0
        for step in range(n_steps):
            if step % save_freq == 0:
                self.trajectory.append(self.save_state(positions, step))
            # Adjust temperature
            if step % 100 == 0 and temperature > min_temperature:
                temperature *= cooling_rate
            
            proposed_positions = self.propose_move(positions)
            proposed_neg_log_likelihood = self.calculate_system_neg_log_likelihood(
                proposed_positions)
            
            # Metropolis criterion for negative log likelihood
            delta = proposed_neg_log_likelihood - current_neg_log_likelihood
            if delta < 0 or np.random.random() < np.exp(-delta / temperature):
                positions = proposed_positions
                current_neg_log_likelihood = proposed_neg_log_likelihood
                accepted += 1
                
                if current_neg_log_likelihood < best_neg_log_likelihood:
                    best_neg_log_likelihood = current_neg_log_likelihood
                    best_positions = {k: v.copy() for k, v in positions.items()}
            
            if step % 100000 == 0:
                acceptance_rate = accepted / (step + 1)
                print(f"Step {step}, Neg Log Likelihood: {current_neg_log_likelihood:.2f}, "
                      f"Temperature: {temperature:.2f}, "
                      f"Acceptance Rate: {acceptance_rate:.2f}")
            
            neg_log_likelihoods.append(current_neg_log_likelihood)

        # Save trajectory to file
        filename = f"{traj_filename}.pkl"
        with open(filename, 'wb') as f:
            pickle.dump(self.trajectory, f)
        
        return best_positions, neg_log_likelihoods

    def visualize_configuration(self, positions, title):
        """Visualize current configuration"""
        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(111, projection='3d')
        
        colors = {'A': 'red', 'B': 'blue', 'C': 'green'}
        
        for type_name, pos in positions.items():
            ax.scatter(pos[:,0], pos[:,1], pos[:,2], 
                      c=colors[type_name], 
                      s=self.radii[type_name]*5,
                      alpha=0.6,
                      label=f'Type {type_name}')
        
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_title(title)
        ax.legend()
        plt.show()

sampler = PairSampler()
    
# Show initial configuration
initial_positions = sampler.initialize_positions()
sampler.visualize_configuration(initial_positions, "Initial Configuration")

sampler.calculate_system_neg_log_likelihood(ideal_coord)
sampler.calculate_system_neg_log_likelihood(initial_positions)

# Run sampling
final_positions, likelihoods = sampler.run_mc(n_steps=10000000, save_freq=10000, traj_filename='trajectory')

sampler.visualize_configuration(final_positions, "Initial Configuration")

sampler.visualize_configuration(ideal_coord, "Initial Configuration")

@dataclass
class SystemParameters:
    """Parameters for the molecular system"""
    box_size: float = 300.0
    radii: Dict[str, float] = None
    pair_distances: Dict[str, float] = None
    component_counts: Dict[str, int] = None
    
    def __post_init__(self):
        if self.radii is None:
            self.radii = {'A': 40.0, 'B': 10.0, 'C': 16.0}
        if self.pair_distances is None:
            self.pair_distances = {
                'AA': 80.5,
                'AB': 50.5,
                'CC': 31.5,
                'BC': 36.5
            }
        if self.component_counts is None:
            self.component_counts = {'A': 8, 'B': 8, 'C': 16}

class BaseMCSampler:
    def __init__(self, params: SystemParameters):
        self.params = params
        self.positions = {}
        self.trajectory = []
        
        # Initialize sigma parameters with Jeffreys prior
        self.sigma = {}
        for pair_type in self.params.pair_distances.keys():
            # Sample from log-uniform distribution between 0.5 and 20 Å
            log_min, log_max = np.log(0.5), np.log(20.0)
            self.sigma[pair_type] = np.exp(np.random.uniform(log_min, log_max))
            
        # Minimum likelihood value to avoid numerical issues
        self.min_likelihood = 1e-10
    
    def initialize_positions(self) -> Dict[str, np.ndarray]:
        """Initialize random positions for all components"""
        positions = {}
        for type_name, count in self.params.component_counts.items():
            positions[type_name] = np.random.uniform(
                0, self.params.box_size, (count, 3))
        return positions

    def calculate_exclusion_score(self, positions: Dict[str, np.ndarray]) -> float:
        """Calculate score term for excluded volume"""
        score = 0.0
        k_exclude = 100.0  # Strong exclusion penalty
        
        for type1 in positions:
            for type2 in positions:
                if type1 <= type2:
                    pos1, pos2 = positions[type1], positions[type2]
                    min_dist = self.params.radii[type1] + self.params.radii[type2]
                    distances = cdist(pos1, pos2)
                    
                    if type1 == type2:
                        # Only consider upper triangle for same type
                        mask = np.triu(np.ones_like(distances), k=1)
                        violations = np.where((distances < min_dist) & (mask > 0))
                    else:
                        violations = np.where(distances < min_dist)
                    
                    if violations[0].size > 0:
                        score += k_exclude * np.sum((min_dist - distances[violations])**2)
        
        return score

    def calculate_pair_score(self, pos1: np.ndarray, pos2: np.ndarray, 
                            target_dist: float, sigma: float) -> float:
        """Calculate score for a pair of positions"""
        distance = np.sqrt(np.sum((pos1 - pos2)**2))
        score = ((distance - target_dist)**2) / (2 * sigma**2)
        return score
    
    def calculate_pair_scores_matrix(self, pos1: np.ndarray, pos2: np.ndarray, 
                                     target_dist: float, sigma: float) -> np.ndarray:
        """Calculate score matrix for all pairs between two sets of positions"""
        distances = cdist(pos1, pos2)
        scores = ((distances - target_dist)**2) / (2 * sigma**2)
        return scores
    
    def propose_sigma_move(self) -> Dict[str, float]:
        """Propose new sigma values"""
        new_sigma = self.sigma.copy()
        # Pick random pair type
        pair_type = np.random.choice(list(self.sigma.keys()))
        
        # Propose new value in log space
        log_sigma = np.log(new_sigma[pair_type])
        step_size = 0.1
        new_log_sigma = log_sigma + np.random.normal(0, step_size)
        
        # Ensure within bounds [0.5, 20]
        new_log_sigma = np.clip(new_log_sigma, np.log(0.5), np.log(20.0))
        new_sigma[pair_type] = np.exp(new_log_sigma)
        
        return new_sigma
    
    def propose_position_move(self, positions: Dict[str, np.ndarray]) -> Dict[str, np.ndarray]:
        """Propose new positions with adaptive step size"""
        new_positions = {k: v.copy() for k, v in positions.items()}
        
        # Select random component type and index
        type_name = np.random.choice(list(self.params.component_counts.keys()))
        idx = np.random.randint(self.params.component_counts[type_name])
        
        # Adaptive step size based on component type
        step_size = self.params.radii[type_name] * 0.1
        
        # Apply random displacement
        new_positions[type_name][idx] += np.random.normal(0, step_size, 3)
        
        # Ensure within box
        new_positions[type_name][idx] = np.clip(
            new_positions[type_name][idx], 0, self.params.box_size)
        
        return new_positions
    
    def save_state(self, step: int, positions: Dict[str, np.ndarray], 
                  score: float) -> Dict:
        """Save current state information"""
        state = {
            'step': step,
            'positions': positions,
            'sigma': self.sigma.copy(),
            'score': score,
            'types': {},
            'bead_numbers': {}
        }
        
        # Add type and bead number information
        current_idx = 0
        for type_name, count in self.params.component_counts.items():
            for i in range(count):
                state['types'][current_idx] = type_name
                state['bead_numbers'][current_idx] = current_idx + 1
                current_idx += 1
        
        return state
    
    def save_trajectory(self, filename: str = None):
        """Save trajectory to file"""
        if filename is None:
            filename = f"trajectory.pkl"
        
        with open(filename, 'wb') as f:
            pickle.dump(self.trajectory, f)
        
        return filename
    
    @staticmethod
    def load_trajectory(filename: str) -> List[Dict]:
        """Load trajectory from file"""
        with open(filename, 'rb') as f:
            return pickle.load(f)
            
    def run_mc(self, n_steps: int = 50000, save_freq: int = 100) -> Tuple:
        """Base Monte Carlo sampling with position and sigma moves"""
        positions = self.initialize_positions()
        current_score = self.calculate_score(positions)
        
        best_positions = {k: v.copy() for k, v in positions.items()}
        best_score = current_score
        
        self.trajectory = []
        accepted = 0
        
        # Temperature schedule
        initial_temp = 10.0
        final_temp = 0.1
        
        for step in range(n_steps):
            # Update temperature
            temp = initial_temp * (final_temp / initial_temp)**(step / n_steps)
            
            # Choose between position (90%) and sigma (10%) moves
            if np.random.random() < 0.9:
                proposed_positions = self.propose_position_move(positions)
                proposed_sigma = self.sigma.copy()
            else:
                proposed_positions = positions.copy()
                proposed_sigma = self.propose_sigma_move()
            
            # Store current state
            old_sigma = self.sigma.copy()
            self.sigma = proposed_sigma
            
            # Calculate scores
            proposed_score = self.calculate_score(proposed_positions)
            delta_e = proposed_score - current_score
            
            # Metropolis criterion
            if delta_e < 0 or np.random.random() < np.exp(-delta_e/temp):
                positions = proposed_positions
                current_score = proposed_score
                accepted += 1
                
                if current_score < best_score:
                    best_score = current_score
                    best_positions = {k: v.copy() for k, v in positions.items()}
            else:
                self.sigma = old_sigma
            
            # Save state and report progress
            if step % save_freq == 0:
                self.trajectory.append(self.save_state(step, positions, current_score))
                print(f"Step {step}, Score: {current_score:.2f}, "
                      f"T: {temp:.4f}, Accept: {accepted/(step+1):.2f}")
        
        return best_positions, self.trajectory, self.save_trajectory()

class PairSampler(BaseMCSampler):
    def __init__(self, params: SystemParameters):
        super().__init__(params)
    
    def calculate_score(self, positions: Dict[str, np.ndarray]) -> float:
        """Calculate score for pair-level interactions"""
        score = 0.0
        
        # Add excluded volume score
        score += self.calculate_exclusion_score(positions)
        
        # Calculate pair scores
        pair_types = [('A', 'A'), ('A', 'B'), ('B', 'C'), ('C', 'C')]
        for type1, type2 in pair_types:
            if f"{type1}{type2}" in self.params.pair_distances:
                target_dist = self.params.pair_distances[f"{type1}{type2}"]
                sigma_value = self.sigma[f"{type1}{type2}"]
                
                score_matrix = self.calculate_pair_scores_matrix(
                    positions[type1],
                    positions[type2] if type1 != type2 else positions[type2],
                    target_dist,
                    sigma_value
                )
                
                if type1 == type2:
                    mask = np.triu(np.ones_like(score_matrix), k=1)
                    score_matrix *= mask
                
                # Get unique pairs and their scores
                row_indices = np.argmin(score_matrix, axis=1)
                col_indices = np.argmin(score_matrix, axis=0)
                
                row_pairs = set((i, row_indices[i]) for i in range(len(row_indices)))
                col_pairs = set((col_indices[j], j) for j in range(len(col_indices)))
                
                unique_pairs = row_pairs.union(col_pairs)
                
                for i, j in unique_pairs:
                    score += score_matrix[i, j]
        
        return score
        
class TetramerSampler(PairSampler):
    def __init__(self, params: SystemParameters):
        super().__init__(params)
        
    def run_mc(self, n_steps: int = 50000, save_freq: int = 100) -> Tuple:
        """Monte Carlo sampling with tetramer moves"""
        positions = self.initialize_positions()
        current_score = self.calculate_score(positions)
        
        best_positions = {k: v.copy() for k, v in positions.items()}
        best_score = current_score
        
        self.trajectory = []
        accepted = 0
        
        # Temperature schedule
        initial_temp = 10.0
        final_temp = 0.1
        
        for step in range(n_steps):
            temp = initial_temp * (final_temp / initial_temp)**(step / n_steps)
            
            # Choose move type: position (40%), sigma (10%), tetramer (50%)
            move_type = np.random.choice(['position', 'sigma', 'tetramer'], 
                                       p=[0.4, 0.1, 0.5])
            
            if move_type == 'position':
                proposed_positions = self.propose_position_move(positions)
                proposed_sigma = self.sigma.copy()
            elif move_type == 'sigma':
                proposed_positions = positions.copy()
                proposed_sigma = self.propose_sigma_move()
            else:  # tetramer move
                proposed_positions = self.propose_tetramer_move(positions)
                proposed_sigma = self.sigma.copy()
            
            old_sigma = self.sigma.copy()
            self.sigma = proposed_sigma
            
            proposed_score = self.calculate_score(proposed_positions)
            delta_e = proposed_score - current_score
            
            if delta_e < 0 or np.random.random() < np.exp(-delta_e/temp):
                positions = proposed_positions
                current_score = proposed_score
                accepted += 1
                
                if current_score < best_score:
                    best_score = current_score
                    best_positions = {k: v.copy() for k, v in positions.items()}
            else:
                self.sigma = old_sigma
            
            if step % save_freq == 0:
                self.trajectory.append(self.save_state(step, positions, current_score))
                print(f"Step {step}, Score: {current_score:.2f}, "
                      f"Temp: {temp:.4f}, Accept: {accepted/(step+1):.2f}")
        
        return best_positions, self.trajectory, self.save_trajectory()

    def get_tetramers(self, positions: Dict[str, np.ndarray]) -> List[Tuple[int, ...]]:
        """Identify ABCC tetramers based on distances"""
        tetramers = []
        
        # For each A component
        for i in range(len(positions['A'])):
            # Find closest B
            dist_AB = cdist([positions['A'][i]], positions['B'])
            b_idx = np.argmin(dist_AB)
            
            # Find two closest C's to this B
            dist_BC = cdist([positions['B'][b_idx]], positions['C'])
            c_indices = np.argsort(dist_BC[0])[:2]
            
            tetramers.append((i, b_idx, c_indices[0], c_indices[1]))
        
        return tetramers
    
    def calculate_tetramer_score(self, positions: Dict[str, np.ndarray], 
                                tetramer: Tuple[int, ...]) -> float:
        """Calculate score for a single tetramer"""
        a_idx, b_idx, c1_idx, c2_idx = tetramer
        score = 0.0
        
        # AB distance
        score += self.calculate_pair_score(
            positions['A'][a_idx],
            positions['B'][b_idx],
            self.params.pair_distances['AB'],
            self.sigma['AB']
        )
        
        # BC distances
        for c_idx in [c1_idx, c2_idx]:
            score += self.calculate_pair_score(
                positions['B'][b_idx],
                positions['C'][c_idx],
                self.params.pair_distances['BC'],
                self.sigma['BC']
            )
            
        # CC distance
        score += self.calculate_pair_score(
            positions['C'][c1_idx],
            positions['C'][c2_idx],
            self.params.pair_distances['CC'],
            self.sigma['CC']
        )
        
        return score
    
    def calculate_score(self, positions: Dict[str, np.ndarray]) -> float:
        """Calculate total score """
        # Start with pair-level score
        score = super().calculate_score(positions)
        
        # Add tetramer score
        tetramers = self.get_tetramers(positions)
        for tetramer in tetramers:
            score += self.calculate_tetramer_score(positions, tetramer)
        
        return score

    def propose_tetramer_move(self, positions: Dict[str, np.ndarray]) -> Dict[str, np.ndarray]:
        """Propose move for tetramer components"""
        new_positions = {k: v.copy() for k, v in positions.items()}
        
        # Select random A component
        a_idx = np.random.randint(len(positions['A']))
        
        # Find closest B
        dist_AB = cdist([positions['A'][a_idx]], positions['B'])
        b_idx = np.argmin(dist_AB)
        
        # Find two closest C's to this B
        dist_BC = cdist([positions['B'][b_idx]], positions['C'])
        c_indices = np.argsort(dist_BC[0])[:2]
        
        # Generate random displacement for whole tetramer
        displacement = np.random.normal(0, 5, 3)
        
        # Apply to all components in tetramer
        new_positions['A'][a_idx] += displacement
        new_positions['B'][b_idx] += displacement
        new_positions['C'][c_indices[0]] += displacement
        new_positions['C'][c_indices[1]] += displacement
        
        return new_positions
        
class OctetSampler(TetramerSampler):
    def __init__(self, params: SystemParameters):
        super().__init__(params)

    def propose_octet_move(self, positions: Dict[str, np.ndarray]) -> Dict[str, np.ndarray]:
        """Propose coordinated move of octet components"""
        new_positions = {k: v.copy() for k, v in positions.items()}
        
        # Select random starting A component
        a_idx = np.random.randint(len(positions['A']) - 1)  # Leave room for second A
        
        # Find closest second A
        dist_AA = cdist([positions['A'][a_idx]], positions['A'])
        dist_AA[0, a_idx] = np.inf  # Exclude self
        second_a_idx = np.argmin(dist_AA)
        
        # For each A, find its tetramer components
        for a_idx in [a_idx, second_a_idx]:
            # Find closest B
            dist_AB = cdist([positions['A'][a_idx]], positions['B'])
            b_idx = np.argmin(dist_AB)
            
            # Find two closest C's to this B
            dist_BC = cdist([positions['B'][b_idx]], positions['C'])
            c_indices = np.argsort(dist_BC[0])[:2]
            
            # Apply displacement to this tetramer
            displacement = np.random.normal(0, 3, 3)  # Smaller displacement for stability
            
            new_positions['A'][a_idx] += displacement
            new_positions['B'][b_idx] += displacement
            new_positions['C'][c_indices[0]] += displacement
            new_positions['C'][c_indices[1]] += displacement
        
        return new_positions
        
    def run_mc(self, n_steps: int = 50000, save_freq: int = 1000) -> Tuple:
        """Monte Carlo sampling with octet moves"""
        positions = self.initialize_positions()
        current_score = self.calculate_score(positions)
        
        best_positions = {k: v.copy() for k, v in positions.items()}
        best_score = current_score
        
        self.trajectory = []
        accepted = 0
        
        # Temperature schedule
        initial_temp = 10.0
        final_temp = 0.1
        
        for step in range(n_steps):
            temp = initial_temp * (final_temp / initial_temp)**(step / n_steps)
            
            # Choose move type: position (30%), sigma (10%), tetramer (30%), octet (30%)
            move_type = np.random.choice(
                ['position', 'sigma', 'tetramer', 'octet'], 
                p=[0.3, 0.1, 0.3, 0.3]
            )
            
            if move_type == 'position':
                proposed_positions = self.propose_position_move(positions)
                proposed_sigma = self.sigma.copy()
            elif move_type == 'sigma':
                proposed_positions = positions.copy()
                proposed_sigma = self.propose_sigma_move()
            elif move_type == 'tetramer':
                proposed_positions = self.propose_tetramer_move(positions)
                proposed_sigma = self.sigma.copy()
            else:  # octet move
                proposed_positions = self.propose_octet_move(positions)
                proposed_sigma = self.sigma.copy()
            
            old_sigma = self.sigma.copy()
            self.sigma = proposed_sigma
            
            proposed_score = self.calculate_score(proposed_positions)
            delta_e = proposed_score - current_score
            
            if delta_e < 0 or np.random.random() < np.exp(-delta_e/temp):
                positions = proposed_positions
                current_score = proposed_score
                accepted += 1
                
                if current_score < best_score:
                    best_score = current_score
                    best_positions = {k: v.copy() for k, v in positions.items()}
            else:
                self.sigma = old_sigma
            
            if step % save_freq == 0:
                self.trajectory.append(self.save_state(step, positions, current_score))
                print(f"Step {step}, Score: {current_score:.2f}, "
                      f"T: {temp:.4f}, Accept: {accepted/(step+1):.2f}")
        
        return best_positions, self.trajectory, self.save_trajectory()

class FullSystemSampler(OctetSampler):
    def __init__(self, params: SystemParameters):
        super().__init__(params)
        
    def run_mc(self, n_steps: int = 50000, save_freq: int = 1000) -> Tuple:
        """Monte Carlo sampling with all move types"""
        positions = self.initialize_positions()
        current_score = self.calculate_score(positions)
        
        best_positions = {k: v.copy() for k, v in positions.items()}
        best_score = current_score
        
        self.trajectory = []
        accepted = 0
        
        # Temperature schedule with slower cooling
        initial_temp = 20.0  # Higher initial temperature
        final_temp = 0.01   # Lower final temperature
        
        for step in range(n_steps):
            temp = initial_temp * (final_temp / initial_temp)**(step / n_steps)
            
            # Choose move type with emphasis on larger scale moves
            move_type = np.random.choice(
                ['position', 'sigma', 'tetramer', 'octet'], 
                p=[0.2, 0.1, 0.3, 0.4]  # More weight on larger moves
            )
            
            if move_type == 'position':
                proposed_positions = self.propose_position_move(positions)
                proposed_sigma = self.sigma.copy()
            elif move_type == 'sigma':
                proposed_positions = positions.copy()
                proposed_sigma = self.propose_sigma_move()
            elif move_type == 'tetramer':
                proposed_positions = self.propose_tetramer_move(positions)
                proposed_sigma = self.sigma.copy()
            else:  # octet move
                proposed_positions = self.propose_octet_move(positions)
                proposed_sigma = self.sigma.copy()
            
            old_sigma = self.sigma.copy()
            self.sigma = proposed_sigma
            
            proposed_score = self.calculate_score(proposed_positions)
            delta_e = proposed_score - current_score
            
            if delta_e < 0 or np.random.random() < np.exp(-delta_e/temp):
                positions = proposed_positions
                current_score = proposed_score
                accepted += 1
                
                if current_score < best_score:
                    best_score = current_score
                    best_positions = {k: v.copy() for k, v in positions.items()}
            else:
                self.sigma = old_sigma
            
            if step % save_freq == 0:
                self.trajectory.append(self.save_state(step, positions, current_score))
                print(f"Step {step}, score: {current_score:.2f}, "
                      f"Temp: {temp:.4f}, Accept: {accepted/(step+1):.2f}, "
                      f"Move: {move_type}")
        
        return best_positions, self.trajectory, self.save_trajectory()

def visualize_3d_configuration(positions, title="System Configuration"):
    """Visualize the current state of the system"""
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection='3d')
    
    # Plot each type with different colors and sizes
    colors = {'A': 'red', 'B': 'blue', 'C': 'green'}
    for type_name, pos in positions.items():
        ax.scatter(pos[:, 0], pos[:, 1], pos[:, 2], 
                  c=colors[type_name], 
                  label=f'Type {type_name}',
                  alpha=0.6)
    
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title(title)
    ax.legend()
    plt.show()

def run_hierarchical_sampling():
    # Initialize system parameters
    params = SystemParameters(
        box_size=300.0,
        radii={'A': 40.0, 'B': 10.0, 'C': 16.0},
        pair_distances={
            'AA': 80.5,
            'AB': 50.5,
            'CC': 31.5,
            'BC': 36.5
        },
        component_counts={'A': 8, 'B': 8, 'C': 16}
    )
    
    # Create samplers for each level
    samplers = {
        'pair': PairSampler(params),
        'tetramer': TetramerSampler(params),
        'octet': OctetSampler(params),
        'full': FullSystemSampler(params)
    }
    
    # Run sampling at each level with visualization
    results = {}
    current_positions = None
    for level, sampler in samplers.items():
        print(f"\nRunning {level}-level sampling...")

        # Initialize positions: use previous best positions if available, otherwise random
        if current_positions is None:
            initial_positions = sampler.initialize_positions()
        else:
            initial_positions = current_positions
            sampler.positions = initial_positions 
        
        # Visualize initial configuration
        visualize_3d_configuration(initial_positions, f"Initial {level} configuration")
        
        # Run sampling
        best_positions, trajectory, traj_file = sampler.run_mc(
            n_steps=100000,  
            save_freq=1000
        )
        
        # Visualize final configuration
        visualize_3d_configuration(best_positions, f"Final {level} configuration")
        
        # Plot scores
        scores = [state['score'] for state in trajectory]
        plt.figure(figsize=(10, 5))
        plt.plot(Scores)
        plt.xlabel('Step')
        plt.ylabel('score')
        plt.title(f'{level} Sampling')
        plt.show()
        
        # Save results
        results[level] = {
            'best_positions': best_positions,
            'trajectory_file': traj_file,
            'final_sigma': sampler.sigma
        }
        
        print(f"{level} sampling complete:")
        print(f"Final sigma values: {sampler.sigma}")
        print(f"Trajectory saved to: {traj_file}")
        
    return results

def analyze_results(results):
    """Analyze sampling results from all levels"""
    for level, result in results.items():
        print(f"\nAnalysis for {level} sampling:")
        print("Final sigma values:")
        for pair_type, sigma in result['final_sigma'].items():
            print(f"  {pair_type}: {sigma:.3f}")
        
        # Load and analyze trajectory
        trajectory = BaseMCSampler.load_trajectory(result['trajectory_file'])
        
        # Calculate acceptance rates
        n_steps = len(trajectory)
        score_changes = np.diff([state['score'] for state in trajectory])
        acceptance_rate = np.sum(score_changes < 0) / len(score_changes)
        print(f"Acceptance rate: {acceptance_rate:.2f}")

# Run hierarchical sampling
print("Starting hierarchical sampling...")
results = run_hierarchical_sampling()
    
