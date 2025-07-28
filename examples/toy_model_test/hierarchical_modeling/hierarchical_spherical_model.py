import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation as R
from tqdm import tqdm

# --- 1. System Definition and Ground Truth ---

def get_ground_truth_config():
    """Defines the ground truth tetrahedral configuration."""
    # Particle radii
    radii = {'A': 1.5, 'B': 1.0}
    
    # Ground truth hyperparameters
    gt_hyperparams = {
        'mu_d': 5.0, 'sigma_d': 0.2,
        'mu_theta': 109.5, 'sigma_theta': 2.0,
        'sigma_noise_d': 0.1, 'sigma_noise_theta': 1.5
    }

    # Generate ground truth latent parameters from hyperparameters
    np.random.seed(0)
    gt_latent_d = np.random.normal(gt_hyperparams['mu_d'], gt_hyperparams['sigma_d'], 4)
    gt_latent_theta = np.random.normal(gt_hyperparams['mu_theta'], gt_hyperparams['sigma_theta'], 6)

    # Create ground truth positions for a perfect tetrahedron
    # This is a slight simplification; in reality, the positions would be consistent
    # with the latent parameters. We use this as a stable reference.
    positions = {'A': np.array([[0., 0., 0.]]), 'B': np.zeros((4, 3))}
    avg_dist = np.mean(gt_latent_d)
    positions['B'][0] = np.array([1, 1, 1]) * avg_dist / np.sqrt(3)
    positions['B'][1] = np.array([1, -1, -1]) * avg_dist / np.sqrt(3)
    positions['B'][2] = np.array([-1, 1, -1]) * avg_dist / np.sqrt(3)
    positions['B'][3] = np.array([-1, -1, 1]) * avg_dist / np.sqrt(3)
    
    return positions, radii, gt_hyperparams, gt_latent_d, gt_latent_theta

def generate_synthetic_data(latent_d, latent_theta, noise_d_sigma, noise_theta_sigma):
    """Generates noisy observations from latent parameters."""
    np.random.seed(42)
    d_obs = np.random.normal(latent_d, noise_d_sigma)
    theta_obs = np.random.normal(latent_theta, noise_theta_sigma)
    return d_obs, theta_obs

# --- 2. Utility functions for MCMC ---

def get_current_geometry(positions):
    """Calculates distances and angles from current particle positions."""
    pos_A = positions['A'][0]
    pos_B = positions['B']
    
    # Distances (A-B)
    distances = np.linalg.norm(pos_B - pos_A, axis=1)
    
    # Angles (B-A-B)
    angles = []
    combs = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    for i, j in combs:
        v1 = pos_B[i] - pos_A
        v2 = pos_B[j] - pos_A
        cosine_angle = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
        angle = np.arccos(np.clip(cosine_angle, -1.0, 1.0))
        angles.append(np.degrees(angle))
        
    return distances, np.array(angles)

def excluded_volume_penalty(positions, radii):
    """Calculates a penalty for overlapping particles."""
    penalty = 0.0
    # Check B-B overlaps
    for i in range(4):
        for j in range(i + 1, 4):
            dist = np.linalg.norm(positions['B'][i] - positions['B'][j])
            min_dist = radii['B'] + radii['B']
            if dist < min_dist:
                penalty += (min_dist - dist)**2
    return penalty * 100.0 # High penalty weight

# --- 3. MCMC Sampler ---

class HierarchicalSampler:
    def __init__(self, init_positions, radii, observed_d, observed_theta):
        self.positions = init_positions
        self.radii = radii
        self.observed_d = observed_d
        self.observed_theta = observed_theta

        # Initialize parameters to be sampled
        self.params = {
            'mu_d': 4.5, 'sigma_d': 0.5,
            'mu_theta': 105.0, 'sigma_theta': 5.0,
            'latent_d': self.observed_d.copy(),
            'latent_theta': self.observed_theta.copy(),
            'sigma_noise_d': 0.2, 'sigma_noise_theta': 2.0
        }
        
        # Priors (weakly informative)
        self.priors = {
            'mu_d': (4.0, 2.0), 'sigma_d': (0.01, 1.0),
            'mu_theta': (100.0, 10.0), 'sigma_theta': (0.01, 10.0),
            'sigma_noise_d': (0.01, 1.0), 'sigma_noise_theta': (0.01, 5.0)
        }

    def log_likelihood(self, current_d, current_theta):
        # Likelihood of observed data given latent parameters
        ll_d = np.sum(-0.5 * ((self.observed_d - self.params['latent_d']) / self.params['sigma_noise_d'])**2)
        ll_theta = np.sum(-0.5 * ((self.observed_theta - self.params['latent_theta']) / self.params['sigma_noise_theta'])**2)
        
        # Likelihood of latent parameters given current geometry
        ll_latent_d = np.sum(-0.5 * ((self.params['latent_d'] - current_d) / 0.1)**2) # Tight spring
        ll_latent_theta = np.sum(-0.5 * ((self.params['latent_theta'] - current_theta) / 1.0)**2) # Tight spring

        return ll_d + ll_theta + ll_latent_d + ll_latent_theta

    def log_prior(self):
        # Hierarchical priors for latent variables
        lp_d = np.sum(-0.5 * ((self.params['latent_d'] - self.params['mu_d']) / self.params['sigma_d'])**2)
        lp_theta = np.sum(-0.5 * ((self.params['latent_theta'] - self.params['mu_theta']) / self.params['sigma_theta'])**2)

        # Hyperpriors
        lp_hyper = 0
        for p, (mu, sigma) in self.priors.items():
            val = self.params.get(p)
            if val is None: continue

            # Using normal for means, log-normal for std devs (implicitly by sampling in log space)
            if 'mu' in p:
                lp_hyper += -0.5 * ((val - mu) / sigma)**2
            elif 'sigma' in p and val > 0:
                lp_hyper += -np.log(val) # Jeffreys prior proxy
        
        return lp_d + lp_theta + lp_hyper

    def run(self, n_steps=20000, burn_in=5000):
        history = {p: [] for p in self.params}
        history['positions'] = []
        
        for i in tqdm(range(n_steps)):
            # --- Propose moves for each parameter ---
            # (A more sophisticated sampler would tune proposal widths)
            
            # 1. Move particles
            old_pos = {k: v.copy() for k, v in self.positions.items()}
            prop_pos = {k: v + np.random.randn(*v.shape) * 0.05 for k, v in self.positions.items()}
            
            old_d, old_theta = get_current_geometry(old_pos)
            prop_d, prop_theta = get_current_geometry(prop_pos)

            old_log_p = self.log_likelihood(old_d, old_theta) - excluded_volume_penalty(old_pos, self.radii)
            prop_log_p = self.log_likelihood(prop_d, prop_theta) - excluded_volume_penalty(prop_pos, self.radii)

            if prop_log_p > old_log_p or np.random.rand() < np.exp(prop_log_p - old_log_p):
                self.positions = prop_pos

            # 2. Sample other parameters (Gibbs-like steps)
            current_d, current_theta = get_current_geometry(self.positions)
            
            for p_name in self.params:
                old_val = self.params[p_name]
                
                # Propose new value
                if isinstance(old_val, float):
                    prop_val = old_val + np.random.randn() * 0.1
                    if 'sigma' in p_name and prop_val <= 0: continue # sigma must be positive
                else: # array
                    prop_val = old_val + np.random.randn(*old_val.shape) * 0.1

                # Calculate posterior ratio
                old_log_post = self.log_likelihood(current_d, current_theta) + self.log_prior()
                self.params[p_name] = prop_val
                prop_log_post = self.log_likelihood(current_d, current_theta) + self.log_prior()
                
                if prop_log_post > old_log_post or np.random.rand() < np.exp(prop_log_post - old_log_post):
                    pass # Keep proposal
                else:
                    self.params[p_name] = old_val # Reject

            # --- Store history ---
            if i > burn_in:
                for p_name in self.params:
                    history[p_name].append(self.params[p_name])
                history['positions'].append(self.positions['B'].copy())
        
        return history

# --- 4. Analysis and Visualization ---

def analyze_results(history, ground_truth):
    n_params = len(ground_truth)
    fig, axes = plt.subplots(n_params, 2, figsize=(12, 2 * n_params))
    
    param_keys = list(ground_truth.keys())

    for i, key in enumerate(param_keys):
        # Trace plot
        ax = axes[i, 0]
        samples = np.array(history[key])
        ax.plot(samples)
        ax.set_ylabel(key)
        ax.set_title(f'Trace plot for {key}')
        if key in ground_truth:
            ax.axhline(ground_truth[key], color='r', linestyle='--', label='Ground Truth')
        ax.legend()

        # Posterior distribution
        ax = axes[i, 1]
        ax.hist(samples, bins=50, density=True, alpha=0.7, label='Posterior')
        ax.set_title(f'Posterior for {key}')
        if key in ground_truth:
            ax.axvline(ground_truth[key], color='r', linestyle='--', label='Ground Truth')
        ax.legend()

    plt.tight_layout()
    plt.savefig("hierarchical_model_posteriors.png")
    plt.show()

# --- 5. Main Execution ---

if __name__ == '__main__':
    # 1. Setup ground truth and generate data
    gt_pos, radii, gt_hyper, gt_d, gt_theta = get_ground_truth_config()
    d_obs, theta_obs = generate_synthetic_data(gt_d, gt_theta, 
                                               gt_hyper['sigma_noise_d'], 
                                               gt_hyper['sigma_noise_theta'])

    # 2. Initialize sampler with random positions
    init_pos = {'A': np.array([[0., 0., 0.]]), 'B': np.random.rand(4, 3) * 10 - 5}
    sampler = HierarchicalSampler(init_pos, radii, d_obs, theta_obs)

    # 3. Run MCMC
    print("Starting MCMC sampling...")
    history = sampler.run(n_steps=30000, burn_in=10000)
    print("MCMC sampling finished.")

    # 4. Analyze and plot results
    print("Analyzing results...")
    analyze_results(history, gt_hyper)
    print("Analysis complete. Check 'hierarchical_model_posteriors.png'")
