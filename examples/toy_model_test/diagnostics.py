import numpy as np
import pandas as pd
import arviz as az
from scipy import stats
from typing import Dict, List, Tuple, Optional

class AdvancedConvergenceDiagnostics:
    """Convergence diagnostics + automatic burn-in estimation."""
    def __init__(self, chains: Dict[str, List[np.ndarray]]):
        # chains: {param_name: [chain_array, ...]} (each chain_array 1D)
        self.chains = chains

    # -------- Public API --------
    def auto_burnin_fraction(self,
                             min_frac: float = 0.05,
                             max_frac: float = 0.6,
                             step: float = 0.05,
                             z_thresh: float = 2.0) -> float:
        """
        Heuristic burn-in detection:
        For candidate start fractions f, compute Geweke z comparing
        first 10% of remaining chain vs last 50% of remaining chain.
        Pick smallest f where |z| < z_thresh for all chains & params.
        """
        candidates = np.arange(min_frac, max_frac + 1e-9, step)
        # Precompute per param per chain arrays
        per_param = self.chains
        if not per_param:
            return 0.1
        for f in candidates:
            ok_all = True
            for param, chain_list in per_param.items():
                for arr in chain_list:
                    n = len(arr)
                    start = int(f * n)
                    if n - start < 50:  # need enough samples
                        ok_all = False
                        break
                    sub = arr[start:]
                    z = self._geweke_z(sub)
                    if np.isnan(z) or abs(z) > z_thresh:
                        ok_all = False
                        break
                if not ok_all:
                    break
            if ok_all:
                return float(f)
        return float(max_frac)

    def compute_summary_dataframe(self, burnin_fraction: float = 0.0) -> pd.DataFrame:
        """Return simple table: param, chains, samples, Rhat, ESS."""
        rows = []
        for param, chain_list in self.chains.items():
            proc = []
            for c in chain_list:
                n = len(c)
                cut = int(burnin_fraction * n)
                if n - cut < 10:
                    continue
                proc.append(c[cut:])
            if len(proc) < 2:
                continue
            # Equal length truncate
            m = min(len(c) for c in proc)
            arr = np.array([c[:m] for c in proc])  # shape (n_chains, m)
            try:
                rhat = float(az.rhat(arr))
            except Exception:
                rhat = np.nan
            ess = float(az.ess(arr)) if m > 20 else np.nan
            rows.append(dict(Parameter=param,
                             N_chains=len(proc),
                             N_samples=m,
                             Rhat=rhat,
                             ESS=ess,
                             ESS_per_iter=ess / m if m and ess and not np.isnan(ess) else np.nan))
        return pd.DataFrame(rows)

    # -------- Internals --------
    def _geweke_z(self, samples: np.ndarray,
                  first_frac: float = 0.1,
                  last_frac: float = 0.5) -> float:
        n = len(samples)
        if n < 40:
            return np.nan
        a_n = int(first_frac * n)
        b_n = int(last_frac * n)
        if a_n < 5 or b_n < 5:
            return np.nan
        first = samples[:a_n]
        last = samples[-b_n:]
        def spec0(x):
            ac = az.autocorr(x)
            neg = np.where(ac < 0)[0]
            cut = neg[0] if len(neg) else len(ac)
            return np.var(x) * (1 + 2 * np.sum(ac[1:cut]))
        mean_diff = np.mean(first) - np.mean(last)
        var = spec0(first)/a_n + spec0(last)/b_n
        if var <= 0:
            return np.nan
        return mean_diff / np.sqrt(var)
        
    def compute_all_diagnostics(self) -> pd.DataFrame:
        """Compute comprehensive convergence diagnostics."""
        results = []
        
        for param_name, chain_list in self.chains.items():
            # Convert to ArviZ InferenceData for advanced diagnostics
            chain_array = self._prepare_chains(chain_list)
            if chain_array is None:
                continue
                
            idata = az.convert_to_inference_data(chain_array)
            
            diagnostics = {
                'Parameter': param_name,
                'N_chains': len(chain_list),
                'N_samples': chain_array.shape[1],
            }
            
            # 1. Split R-hat (more robust than traditional R-hat)
            diagnostics['Split_Rhat'] = self._compute_split_rhat(chain_array)
            
            # 2. Rank-normalized R-hat (better for heavy-tailed distributions)
            diagnostics['Rank_Rhat'] = self._compute_rank_rhat(chain_array)
            
            # 3. ESS Bulk and Tail (separate measures for center and tails)
            ess_bulk, ess_tail = self._compute_ess_bulk_tail(idata, param_name)
            diagnostics['ESS_Bulk'] = ess_bulk
            diagnostics['ESS_Tail'] = ess_tail
            diagnostics['ESS_per_sec'] = self._compute_ess_per_second(ess_bulk)
            
            # 4. Monte Carlo Standard Error
            diagnostics['MCSE_mean'] = self._compute_mcse_mean(chain_array)
            diagnostics['MCSE_quantiles'] = self._compute_mcse_quantiles(chain_array)
            
            # 5. Geweke diagnostic (test for stationarity)
            diagnostics['Geweke_zscore'] = self._compute_geweke(chain_array)
            diagnostics['Geweke_pvalue'] = 2 * (1 - stats.norm.cdf(abs(diagnostics['Geweke_zscore'])))
            
            # 6. Heidelberger-Welch test (stationarity and halfwidth)
            hw_stat, hw_pval = self._heidelberger_welch_test(chain_array)
            diagnostics['HW_stationarity'] = hw_stat
            diagnostics['HW_pvalue'] = hw_pval
            
            # 7. Raftery-Lewis diagnostic (required sample size)
            diagnostics['RL_required_samples'] = self._raftery_lewis(chain_array)
            
            # 8. Effective sample size per iteration
            diagnostics['Efficiency'] = ess_bulk / chain_array.shape[1]
            
            # 9. Between and within chain variance
            between_var, within_var = self._compute_variance_components(chain_array)
            diagnostics['Between_chain_var'] = between_var
            diagnostics['Within_chain_var'] = within_var
            diagnostics['Var_ratio'] = between_var / within_var if within_var > 0 else np.inf
            
            # 10. Folded rank statistics (for detecting pathological behavior)
            diagnostics['Folded_rank_stat'] = self._compute_folded_rank_stat(chain_array)
            
            results.append(diagnostics)
            
        return pd.DataFrame(results)
    
    def _prepare_chains(self, chain_list: List[np.ndarray]) -> Optional[np.ndarray]:
        """Prepare chains for analysis, handling different lengths."""
        if not chain_list:
            return None
            
        # Find minimum length and truncate all chains
        min_length = min(len(chain) for chain in chain_list)
        if min_length < 100:  # Minimum required samples
            warnings.warn(f"Chain too short ({min_length} samples)")
            return None
            
        # Stack chains: shape (n_chains, n_samples)
        return np.array([chain[:min_length] for chain in chain_list])
    
    def _compute_split_rhat(self, chains: np.ndarray) -> float:
        """Compute split R-hat (more robust than traditional R-hat)."""
        n_chains, n_samples = chains.shape
        
        # Split each chain in half
        split_chains = []
        half = n_samples // 2
        for chain in chains:
            split_chains.extend([chain[:half], chain[half:2*half]])
        
        # Compute R-hat on split chains
        split_array = np.array(split_chains)
        return float(az.rhat(split_array))
    
    def _compute_rank_rhat(self, chains: np.ndarray) -> float:
        """Compute rank-normalized R-hat."""
        # Rank transform all samples
        all_samples = chains.flatten()
        ranks = stats.rankdata(all_samples).reshape(chains.shape)
        
        # Normalize ranks to (0, 1)
        normalized_ranks = (ranks - 0.5) / len(all_samples)
        
        # Transform to normal quantiles
        z_scores = stats.norm.ppf(normalized_ranks)
        z_scores = np.clip(z_scores, -10, 10)  # Avoid infinities
        
        return float(az.rhat(z_scores))
    
    def _compute_ess_bulk_tail(self, idata, param_name: str) -> Tuple[float, float]:
        """Compute ESS for bulk and tail of distribution."""
        try:
            ess_bulk = float(az.ess(idata, method='bulk'))
            ess_tail = float(az.ess(idata, method='tail'))
        except:
            # Fallback to basic ESS
            ess_bulk = float(az.ess(idata))
            ess_tail = ess_bulk
        return ess_bulk, ess_tail
    
    def _compute_ess_per_second(self, ess: float, runtime_seconds: float = 1.0) -> float:
        """Compute computational efficiency (ESS per second)."""
        return ess / runtime_seconds
    
    def _compute_mcse_mean(self, chains: np.ndarray) -> float:
        """Compute Monte Carlo Standard Error for the mean."""
        all_samples = chains.flatten()
        return float(az.mcse(all_samples, method='mean'))
    
    def _compute_mcse_quantiles(self, chains: np.ndarray, quantiles=[0.025, 0.975]) -> Dict[float, float]:
        """Compute MCSE for quantiles."""
        all_samples = chains.flatten()
        mcse_dict = {}
        for q in quantiles:
            mcse_dict[q] = float(az.mcse(all_samples, method='quantile', prob=q))
        return mcse_dict
    
    def _compute_geweke(self, chains: np.ndarray, first_frac: float = 0.1, 
                       last_frac: float = 0.5) -> float:
        """Compute Geweke diagnostic comparing means of first and last parts."""
        all_samples = chains.flatten()
        n = len(all_samples)
        
        first_n = int(n * first_frac)
        last_n = int(n * last_frac)
        
        first_samples = all_samples[:first_n]
        last_samples = all_samples[-last_n:]
        
        # Compute means and spectral density estimates
        mean_diff = np.mean(first_samples) - np.mean(last_samples)
        
        # Use spectral density at zero for variance estimates
        var_first = self._spectral_density_at_zero(first_samples)
        var_last = self._spectral_density_at_zero(last_samples)
        
        se = np.sqrt(var_first/len(first_samples) + var_last/len(last_samples))
        
        return mean_diff / se if se > 0 else 0.0
    
    def _spectral_density_at_zero(self, x: np.ndarray) -> float:
        """Estimate spectral density at frequency zero."""
        # Simple implementation using initial positive sequence estimator
        acf = az.autocorr(x)
        
        # Find first negative autocorrelation
        first_negative = np.where(acf < 0)[0]
        if len(first_negative) > 0:
            cutoff = first_negative[0]
        else:
            cutoff = len(acf)
        
        # Initial positive sequence estimator
        return np.var(x) * (1 + 2 * np.sum(acf[1:cutoff]))
    
    def _heidelberger_welch_test(self, chains: np.ndarray, alpha: float = 0.05) -> Tuple[bool, float]:
        """
        Heidelberger-Welch test for stationarity.
        Returns (is_stationary, p_value)
        """
        all_samples = chains.flatten()
        n = len(all_samples)
        
        # Use Cramér-von Mises test on cumulative means
        cumsum = np.cumsum(all_samples)
        cum_mean = cumsum / np.arange(1, n + 1)
        
        # Brownian bridge transformation
        bridge = cum_mean - cum_mean[-1] * np.arange(1, n + 1) / n
        
        # Cramér-von Mises statistic
        cvm_stat = np.sum(bridge**2) / (n * np.var(all_samples))
        
        # Approximate p-value (simplified)
        p_value = 1 - stats.chi2.cdf(cvm_stat, df=1)
        
        return p_value > alpha, p_value
    
    def _raftery_lewis(self, chains: np.ndarray, quantile: float = 0.025, 
                       accuracy: float = 0.005, probability: float = 0.95) -> int:
        """
        Simplified Raftery-Lewis diagnostic for required sample size.
        """
        all_samples = chains.flatten()
        
        # Binary chain: 1 if sample < quantile, 0 otherwise
        q_val = np.quantile(all_samples, quantile)
        binary = (all_samples < q_val).astype(int)
        
        # Estimate transition probabilities
        transitions = np.column_stack([binary[:-1], binary[1:]])
        p00 = np.mean((transitions[:, 0] == 0) & (transitions[:, 1] == 0))
        p11 = np.mean((transitions[:, 0] == 1) & (transitions[:, 1] == 1))
        
        # Avoid division by zero
        if p00 == 1 or p11 == 1:
            return len(all_samples) * 10  # Need much more samples
        
        # Required sample size (simplified formula)
        z = stats.norm.ppf((1 + probability) / 2)
        required = int((z**2 * quantile * (1 - quantile)) / (accuracy**2))
        
        # Adjust for autocorrelation
        autocorr_factor = (1 + p00 + p11) / (1 - p00 - p11 + 1e-10)
        required = int(required * max(1, autocorr_factor))
        
        return required
    
    def _compute_variance_components(self, chains: np.ndarray) -> Tuple[float, float]:
        """Compute between-chain and within-chain variance."""
        n_chains, n_samples = chains.shape
        
        chain_means = np.mean(chains, axis=1)
        grand_mean = np.mean(chain_means)
        
        # Between-chain variance
        B = n_samples * np.var(chain_means, ddof=1)
        
        # Within-chain variance
        W = np.mean([np.var(chain, ddof=1) for chain in chains])
        
        return B, W
    
    def _compute_folded_rank_stat(self, chains: np.ndarray) -> float:
        """
        Compute folded rank statistic to detect pathological behavior.
        Values far from 0.5 indicate problems.
        """
        all_samples = chains.flatten()
        ranks = stats.rankdata(all_samples)
        
        # Fold ranks around median
        median_rank = len(ranks) / 2
        folded_ranks = np.minimum(ranks, len(ranks) + 1 - ranks)
        
        # Normalize to [0, 1]
        normalized = folded_ranks / median_rank
        
        # Should be uniformly distributed if chains mixed well
        # Use Kolmogorov-Smirnov test against uniform
        ks_stat, _ = stats.kstest(normalized, 'uniform')
        
        return ks_stat
    
    def generate_diagnostic_report(self, df: pd.DataFrame) -> str:
        """Generate a text report with interpretation of diagnostics."""
        report = []
        report.append("=" * 80)
        report.append("MCMC CONVERGENCE DIAGNOSTIC REPORT")
        report.append("=" * 80)
        
        for _, row in df.iterrows():
            report.append(f"\nParameter: {row['Parameter']}")
            report.append("-" * 40)
            
            # Check convergence criteria
            issues = []
            warnings = []
            
            # R-hat checks
            if row['Split_Rhat'] > 1.01:
                issues.append(f"Split R-hat = {row['Split_Rhat']:.3f} (should be < 1.01)")
            elif row['Split_Rhat'] > 1.005:
                warnings.append(f"Split R-hat = {row['Split_Rhat']:.3f} (marginal)")
            
            # ESS checks
            min_ess = 400  # Minimum recommended
            if row['ESS_Bulk'] < min_ess:
                issues.append(f"ESS Bulk = {row['ESS_Bulk']:.0f} (should be > {min_ess})")
            if row['ESS_Tail'] < min_ess:
                issues.append(f"ESS Tail = {row['ESS_Tail']:.0f} (should be > {min_ess})")
            
            # Efficiency check
            if row['Efficiency'] < 0.01:
                warnings.append(f"Low efficiency = {row['Efficiency']:.3f}")
            
            # Geweke test
            if row['Geweke_pvalue'] < 0.05:
                warnings.append(f"Geweke test suggests non-stationarity (p={row['Geweke_pvalue']:.3f})")
            
            # Heidelberger-Welch test
            if not row['HW_stationarity']:
                issues.append(f"Failed Heidelberger-Welch stationarity test")
            
            # Required samples check
            if row['RL_required_samples'] > row['N_samples'] * row['N_chains']:
                warnings.append(f"May need {row['RL_required_samples']} total samples")
            
            # Report findings
            if issues:
                report.append("CONVERGENCE ISSUES:")
                for issue in issues:
                    report.append(f"   - {issue}")
            elif warnings:
                report.append("WARNINGS:")
                for warning in warnings:
                    report.append(f"   - {warning}")
            else:
                report.append("CONVERGENCE LOOKS GOOD")
            
            # Summary statistics
            report.append(f"\nSummary:")
            report.append(f"  Chains: {row['N_chains']}, Samples per chain: {row['N_samples']}")
            report.append(f"  ESS Bulk: {row['ESS_Bulk']:.0f}, ESS Tail: {row['ESS_Tail']:.0f}")
            report.append(f"  Efficiency: {row['Efficiency']:.1%}")
            
        report.append("\n" + "=" * 80)
        return "\n".join(report)