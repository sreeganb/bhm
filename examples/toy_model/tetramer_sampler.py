

class TetramerSampler(PairSampler):
    def __init__(
        self,
        params: SystemParameters,
        use_sigma_distribution=True
    ):
        super().__init__(params)
        self.use_sigma_distribution = use_sigma_distribution

        # Additional initialization for tetramer-specific features
        self.tetramer_trans_step = 0.5
        self.tetramer_rot_step = 0.15
        self.sigma_prior_dist = {}
        #-------------------------------------
        # tracking acceptance rates
        self.tet_trans_acc_rate = 0.5
        self.tet_rot_acc_rate = 0.5
        self.target_acceptance = 0.4
        #-------------------------------------
        self.sigma_history_tetramer = {
            "AA": [],
            "AB": [],
            "BC": [],
            "CC": []
        }

        required_pairs = {"AA", "AB", "BC", "CC"}
        if not required_pairs.issubset(self.params.pair_distances.keys()):
            raise ValueError("Missing required pair types for tetramer sampling.")

        if use_sigma_distribution:
            # Load pre-fit GMM parameters from the JSON file
            gmm_file = "output_analysis/gmm_parameters.json"
            print(f"Loading GMM parameters from: {gmm_file}")
            try:
                with open(gmm_file, "r") as f:
                    self.gmm_data = json.load(f)
                print("GMM parameters loaded successfully.")
            except FileNotFoundError:
                print(f"Error: GMM parameters file not found: {gmm_file}")
                raise
            except json.JSONDecodeError:
                print(f"Error: Could not decode JSON in {gmm_file}")
                raise
            except Exception as e:
                print(f"An unexpected error occurred while loading GMM: {e}")
                raise

            # Initialize self.sigma by sampling from the GMM
            self.initialize_sigma_from_gmm()
        else:
            # Normal PairSampler initialization logic
            self.initialize_sigma()  # Initialize sigma values from BaseMCSampler
            print("Using default sigma initialization (no GMM).")

    def initialize_sigma_from_gmm(self):
        """
        For each required pair type, sample an initial sigma value from the GMM.
        """
        print("Initializing sigma from GMM...")

        def sample_from_gmm_1d(means, covariances, mixture_weights):
            # Normalize the mixture_weights so they sum to 1
            mixture_weights = np.array(mixture_weights)
            mixture_weights = mixture_weights / mixture_weights.sum()
            print(f"  Mixture weights (normalized): {mixture_weights}")

            comp_idx = np.random.choice(range(len(means)), p=mixture_weights)
            print(f"  Selected component index: {comp_idx}")
            mean = means[comp_idx]
            var = covariances[comp_idx]
            print(f"  Selected component mean: {mean}, variance: {var}")
             #check variance
            if var <= 0:
                print("Warning: Variance is zero or negative. Setting candidate to 1.0")
                return 1.0
            candidate = np.random.normal(mean, np.sqrt(var))
            print(f"  Sampled candidate: {candidate}")
            return candidate if candidate > 0 else 1.0

        for pair_type in ["AA", "AB", "BC", "CC"]:
            if pair_type not in self.gmm_data:
                # Fallback if no GMM is present for this type
                print(f"  No GMM data found for {pair_type}. Setting sigma to default 1.0")
                self.sigma[pair_type] = 1.0
                continue

            gmm_params = self.gmm_data[pair_type]
            means = gmm_params["means"]
            covariances = gmm_params["covariances"]
            mixture_weights = gmm_params["mixture_weights"]
            print(f"Sampling sigma for {pair_type} from GMM:")

            # Sample a value from the GMM
            self.sigma[pair_type] = sample_from_gmm_1d(means, covariances, mixture_weights)
            print(f"  Initialized sigma[{pair_type}] = {self.sigma[pair_type]}")

        #print("Sigma initialization complete.")

    def compute_prior_penalty(self, proposed_sigma: dict) -> float:
        """
        Computes a penalty (negative log-likelihood) based on the GMM PDF.
        """
        # Return 0 if not using GMM
        if not self.use_sigma_distribution:
        #    print("compute_prior_penalty: GMM prior not used. Returning 0.0")
            return 0.0

        def pdf_1d_gmm(x, means, covariances, mixture_weights):
            pdf_val = 0.0
            for (m, c, w) in zip(means, covariances, mixture_weights):
                if c <= 0:
                    continue  # Skip components with non-positive variance
                norm_const = 1.0 / (math.sqrt(2.0 * math.pi * c))
                exponent = -((x - m)**2) / (2.0 * c)
                pdf_val += w * norm_const * math.exp(exponent)
            return pdf_val

        penalty = 0.0
        #print("Calculating prior penalty...")
        for param, val in proposed_sigma.items():
        #    print(f"  Checking prior for {param} = {val}")
            # If we didn't fit a GMM for this parameter, skip
            if param not in self.gmm_data:
                print(f"    No GMM data for {param}. Skipping.")
                continue

            gmm_info = self.gmm_data[param]
            prob = pdf_1d_gmm(
                val,
                gmm_info["means"],
                gmm_info["covariances"],
                gmm_info["mixture_weights"]
            )
        #    print(f"    GMM PDF value at {val}: {prob}")

            if prob < 1e-12:
            #    print(f"    Probability too small. Applying large penalty.")
                penalty += 9999  # Or some other large value
            else:
                penalty -= math.log(prob)
            #    print(f"    -log(prob): {-math.log(prob)}, Accumulated penalty: {penalty}")

        #print(f"Total prior penalty: {penalty}")
        return penalty
    
    def propose_tetramer_sigma_move(self):
        """Proposes a new sigma value in log-space, with clipping."""
        proposed_sigma = self.sigma.copy()  # Start with a copy of the current sigmas

        if self.use_sigma_distribution:
            pair_types = list(self.gmm_data.keys())
        else:
            pair_types = list(self.sigma.keys())

        total_combs = 488.0
        # Select pair type based on population
        pair_type = np.random.choice(pair_types, p=[56.0/total_combs, 64.0/total_combs, 128.0/total_combs, 240.0/total_combs])
        current_sigma = self.sigma[pair_type]

        # --- Key Changes: Log-space proposal and clipping ---
        log_current_sigma = np.log(current_sigma)  # Work in log-space
        log_proposal = log_current_sigma + np.random.normal(0, 0.1)

        # Clipping:  Define reasonable bounds for sigma (adjust as needed)
        sigma_min = 0.01 # Example minimum
        sigma_max = 5.0 # Example Maximum
        proposal = np.clip(np.exp(log_proposal), sigma_min, sigma_max)  # Clip after exponentiating

        proposed_sigma[pair_type] = proposal
        # --- End of Key Changes ---

        #print(f"Proposing sigma move: {pair_type} from {current_sigma:.4f} to {proposed_sigma[pair_type]:.4f} "
        #      f"(log: {log_current_sigma:.4f} to {np.log(proposal):.4f})")
        return proposed_sigma, pair_type
                
    def _validate_positions(self, positions):
        """Ensure required components exist"""
        for comp in ['A', 'B', 'C']:
            if comp not in positions:
                raise ValueError(f"Missing {comp} particles in initial positions")
        return positions
    
    def read_positions(self):
        pkl_files = [f for f in os.listdir("output_analysis") if f.startswith("trajectory_PairSampler_chain_") and f.endswith(".pkl")]
        pair_pkl_file = os.path.join("output_analysis", random.choice(pkl_files))
        
        with open(pair_pkl_file, "rb") as f:
            trajectory = pickle.load(f)
        # Skip the first 80% of the steps
        skip_steps = int(0.8 * len(trajectory))
        trajectory = trajectory[skip_steps:]
        # Choose one step at random
        step_data = random.choice(trajectory)
        visualize_3d_configuration(step_data["positions"], self.params.radii, f"Initial tetramer configuration")
        
        return step_data["positions"]
        
    def run_mc(self, n_steps: int = 50000, save_freq: int = 100, pos_read: bool = False) -> Tuple:
        """Monte Carlo sampling with tetramer moves."""
        if pos_read:
            positions = self.read_positions()
        else:
            positions = self.initialize_positions()
        current_score = self.calculate_score(positions, exclusion_weight=1.0, 
                                             pair_weight=1.0, tetramer_weight=1.0)
        best_positions = {k: v.copy() for k, v in positions.items()}
        best_score = current_score
        self.trajectory = []
        accepted = 0
        initial_temp = 5.0
        final_temp = 1.0
        
        tet_moves = 0
        tet_accepted = 0  # Single counter for tetramer moves (combined trans+rot)
        
        for step in range(n_steps):
            temp = initial_temp * (final_temp / initial_temp)**(step / n_steps)
            move_type = np.random.choice(['position', 'sigma', 'tetramer'], 
                                       p=[0.4, 0.1, 0.5])
            if move_type == 'position':
                proposed_positions = self.propose_position_move(positions)
                proposed_sigma = self.sigma.copy()
            elif move_type == 'sigma':
                proposed_positions = positions.copy()
                if self.use_sigma_distribution:
                    #propose a sigma move based on the GMM
                    proposed_sigma, pair_type = self.propose_tetramer_sigma_move()
                else:
                    proposed_sigma, pair_type = self.propose_sigma_move()
            else:  # tetramer move
                tet_moves += 1
                self.tet_trans_acc_rate = tet_accepted / max(1, tet_moves)
                self.tet_rot_acc_rate = tet_accepted / max(1, tet_moves)  # Same rate for simplicity
                proposed_positions = self.propose_tetramer_move(positions)
                proposed_sigma = self.sigma.copy()
            old_sigma = self.sigma.copy()
            self.sigma = proposed_sigma
            prior_penalty_from_distribution = 0.0
            if self.use_sigma_distribution:
                # calculate the prior penalty
                prior_penalty_from_distribution = self.compute_prior_penalty(proposed_sigma)
            proposed_score = self.calculate_score(proposed_positions, 
                                                  prior_penalty_from_distribution)
            #print("proposed score: ", proposed_score)
            delta_e = proposed_score - current_score
            acceptance = 0
            if delta_e < 0 or np.random.random() < np.exp(-delta_e/temp):
                positions = proposed_positions
                current_score = proposed_score
                accepted += 1
                acceptance = 1
                if current_score < best_score:
                    best_score = current_score
                    best_positions = {k: v.copy() for k, v in positions.items()}
            else:
                self.sigma = old_sigma
            
            # Ensure the output directory exists
            os.makedirs("output_analysis", exist_ok=True)
            csv_log_file = "output_analysis/all_info_mcmc_tetramer.csv"  # Different filename
            accept_rate = accepted / (step + 1)

            if step % save_freq == 0:    
                # After accept/reject, record the current sigma in the dictionary
                self.sigma_history_tetramer["AA"].append(self.sigma["AA"])
                self.sigma_history_tetramer["AB"].append(self.sigma["AB"])
                self.sigma_history_tetramer["BC"].append(self.sigma["BC"])
                self.sigma_history_tetramer["CC"].append(self.sigma["CC"])

                # --- Saving Logic (Consolidated and Improved) ---
                if step == 0:
                    with open(csv_log_file, "w") as log_file:
                        header = "Step,T,Score,Accepted\n"  # Simpler header for tetramer
                        log_file.write(header)

                with open(csv_log_file, "a") as log_file:
                    log_file.write(f"{step},{temp:.4f},{current_score:.4f},{acceptance:.4f}\n")

                self.trajectory.append(self.save_state(step, positions, current_score))

                with open("total_score_log_tetramer.csv", "a") as log_file: #separate file for score
                    log_file.write(f"{step},{current_score}\n")
                
                print(f"Step {step}, Score: {current_score:.2f}, "
                      f"Temp: {temp:.4f}, Accept: {accept_rate:.2f}")

        # Save sigma history to CSV (Corrected: use sigma_history_tetramer)
        os.makedirs("output_analysis", exist_ok=True)
        df = pd.DataFrame(self.sigma_history_tetramer)
        df.to_csv("output_analysis/sigma_history_tetramer.csv", index=False)

        return best_positions, self.trajectory, self.save_trajectory()
    
    def get_tetramers(self, positions: Dict[str, np.ndarray], temp=1.0) -> List[Tuple[int, ...]]:
        """Safer probability calculation with validation"""
        tetramers = []
        n_b = len(positions['B'])
        n_c = len(positions['C'])
        
        for a_idx in range(len(positions['A'])):
            # B selection with fallback
            dist_AB = cdist([positions['A'][a_idx]], positions['B'])[0]
            probs_B = np.exp(-dist_AB / temp)
            probs_B = np.nan_to_num(probs_B, nan=1/n_b)  # Handle NaNs
            if probs_B.sum() == 0:
                probs_B = np.ones(n_b)/n_b
            probs_B /= probs_B.sum()
            
            b_idx = np.random.choice(n_b, p=probs_B)

            # C selection with validation
            dist_BC = cdist([positions['B'][b_idx]], positions['C'])[0]
            probs_C = np.exp(-dist_BC / temp)
            probs_C = np.nan_to_num(probs_C, nan=1/n_c)
            if probs_C.sum() == 0 or len(probs_C) < 2:
                c_indices = np.random.choice(n_c, size=2, replace=(n_c < 2))
            else:
                probs_C /= probs_C.sum()
                c_indices = np.random.choice(n_c, size=2, replace=False, p=probs_C)
                
            tetramers.append((a_idx, b_idx, c_indices[0], c_indices[1]))
        
        return tetramers
    
    def propose_tetramer_move(self, positions: Dict[str, np.ndarray]) -> Dict[str, np.ndarray]:
        """
        Propose a coordinated move for an entire ABCC tetramer with adaptive translation and rotation.
        
        Strategy:
        1. Select one tetramer probabilistically using get_tetramers
        2. Apply coordinated translation with adaptive step size
        3. Apply rotation about tetramer's geometric center with adaptive step size
        
        Args:
            positions: Current positions dictionary with particle coordinates by type
            
        Returns:
            New positions dictionary with updated tetramer coordinates
        """
        # Create a deep copy of the current positions to avoid modifying the original
        new_pos = {k: v.copy() for k, v in positions.items()}
        
        # Get all possible tetramers using probabilistic selection based on proximity
        tetramers = self.get_tetramers(positions, temp=1.0)  # temp controls selection softness
        
        # If no tetramers are found, return the unchanged positions
        if not tetramers:
            return new_pos
        
        # Randomly select one tetramer from the list for this move
        a_idx, b_idx, c_idx1, c_idx2 = tetramers[np.random.randint(len(tetramers))]
        c_indices = [c_idx1, c_idx2]  # List of C particle indices for convenience
        
        # Calculate adaptive step size for translation based on acceptance rate
        # Adjustment factor scales between 0.5x and 2x, targeting self.target_acceptance
        trans_adjust = np.clip(1.0 + 2.0 * (self.tet_trans_acc_rate - self.target_acceptance), 0.5, 2.0)
        trans_step = self.tetramer_trans_step * trans_adjust  # Scale base step size
        
        # Calculate adaptive step size for rotation using the same logic
        rot_adjust = np.clip(1.0 + 2.0 * (self.tet_rot_acc_rate - self.target_acceptance), 0.5, 2.0)
        rot_step = self.tetramer_rot_step * rot_adjust  # Scale base rotation step (in radians)
        
        # Generate a random 3D displacement vector for translation
        displacement = np.random.normal(0, trans_step, 3)  # Mean 0, std dev trans_step
        
        # Apply the same displacement to all tetramer components to maintain structure
        for part, idx in [('A', a_idx), ('B', b_idx)] + [('C', c_idx) for c_idx in c_indices]:
            new_pos[part][idx] += displacement  # Translate each particle
        
        # Collect coordinates of the tetramer after translation for rotation
        tetramer_coords = np.array([new_pos['A'][a_idx], new_pos['B'][b_idx]] + 
                                [new_pos['C'][c_idx] for c_idx in c_indices])
        
        # Calculate the geometric center (centroid) of the tetramer
        centroid = np.mean(tetramer_coords, axis=0)
        
        # Generate a random rotation axis (unit vector)
        rotation_axis = np.random.randn(3)  # Random direction in 3D space
        rotation_axis /= np.linalg.norm(rotation_axis)  # Normalize to ensure unit length
        
        # Generate a random rotation angle based on adaptive step size
        rotation_angle = np.random.normal(0, rot_step)  # Angle in radians, mean 0, std dev rot_step
        
        # Compute the rotation matrix using a helper function (e.g., Rodrigues' formula)
        rot_matrix = self._rotation_matrix(rotation_axis, rotation_angle)
        
        # Apply rotation about the centroid to all tetramer components
        for part, idx in [('A', a_idx), ('B', b_idx)] + [('C', c_idx) for c_idx in c_indices]:
            vec = new_pos[part][idx] - centroid  # Vector from centroid to particle
            new_pos[part][idx] = centroid + rot_matrix @ vec  # Rotate and reposition
            # Enforce boundary conditions to keep particles within the simulation box
            new_pos[part][idx] = np.clip(new_pos[part][idx], 0, self.params.box_size)
        
        # Return the updated positions dictionary with the moved tetramer
        return new_pos

    def _rotation_matrix(self, axis: np.ndarray, theta: float) -> np.ndarray:
        """
        Create 3D rotation matrix about arbitrary axis using Rodrigues' formula.
        
        Args:
            axis: Unit vector specifying rotation axis (3D np.array)
            theta: Rotation angle in radians
            
        Returns:
            3x3 rotation matrix
        """
        # Identity matrix
        I = np.eye(3)
        # Cross product matrix for axis
        K = np.array([
            [0, -axis[2], axis[1]],
            [axis[2], 0, -axis[0]],
            [-axis[1], axis[0], 0]
        ])
        return I + np.sin(theta)*K + (1-np.cos(theta))*(K @ K)
        
    def calculate_tetramer_score(self, positions: Dict[str, np.ndarray], 
                                tetramer: Tuple[int, ...]) -> float:
        """Calculate score for a single tetramer."""
        a_idx, b_idx, c1_idx, c2_idx = tetramer
        score = 0.0
        score += self.pair_score_nll(
            positions['A'][a_idx],
            positions['B'][b_idx],
            self.params.pair_distances['AB'],
            self.sigma['AB']
        )
        for c_idx in [c1_idx, c2_idx]:
            score += self.pair_score_nll(
                positions['B'][b_idx],
                positions['C'][c_idx],
                self.params.pair_distances['BC'],
                self.sigma['BC']
            )
        score += self.pair_score_nll(
            positions['C'][c1_idx],
            positions['C'][c2_idx],
            self.params.pair_distances['CC'],
            self.sigma['CC']
        )
        return score
    
    def calculate_score(
        self,
        positions: Dict[str, np.ndarray],
        prior_penalty_from_distribution: float = 0.0,
        exclusion_weight: float = 1.0,
        pair_weight: float = 1.0,
        tetramer_weight: float = 1.0,
        excluded_pairs=None
    ) -> float:
        """
        Calculate the total score for a tetramer system:
        1. Identify the A-B-C-C tetramers and collect the pairs inside each tetramer.
        2. Exclude those pairs from the parent PairSampler calculation.
        3. Add an additional term for tetramer-specific scoring.
        """
        # 1) Identify pairs within tetramers
        tetramers = self.get_tetramers(positions)  # Method that finds tuples like (a_idx, b_idx, c1_idx, c2_idx)
        tetramer_pairs = set()
        for a_idx, b_idx, c1_idx, c2_idx in tetramers:
            tetramer_pairs.add(('A', a_idx, 'B', b_idx))
            tetramer_pairs.add(('B', b_idx, 'C', c1_idx))
            tetramer_pairs.add(('B', b_idx, 'C', c2_idx))
            tetramer_pairs.add(('C', c1_idx, 'C', c2_idx))

        # 2) Calculate score using the parent PairSampler, excluding tetramer-internal pairs
        #    Note the named argument `sigma=self.sigma` to ensure the parent sees a dictionary for sigma.
        score, _, _, _ = super().calculate_score(
            positions=positions,
            sigma=self.sigma,  # Pass the dictionary so parent can do sigma[pair_key]
            exclusion_weight=exclusion_weight,
            pair_weight=pair_weight,
            excluded_pairs=tetramer_pairs,
            use_sigma_distribution=self.use_sigma_distribution,
            prior_penalty_from_distribution=prior_penalty_from_distribution
        )

        # 3) Add tetramer-specific score
        for tetramer in tetramers:
            score += tetramer_weight * self.calculate_tetramer_score(positions, tetramer)

        return score