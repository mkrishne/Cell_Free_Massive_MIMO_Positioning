This repository contains the code associated with the paper "Distributed Machine Learning Approach for Low-Latency Localization in Cell-Free Massive MIMO Systems."

Following MATLAB toolboxes are needed to run the code :
	1. Signal Processing Toolbox
	2. Phased Array Toolbox
	3. Statistics and Machine Learning Toolbox
	4. Deep Learning Toolbox
	5. Parallel Computing Toolbox

The simulations were executed in MATLAB, and the results were saved using the diary log. Python was then used to generate the plots for greater control over the plotting environment. An exception is Fig. 9 (CDF plot), which was generated directly in MATLAB. Fig. 10 results were obtained from the CRLB outputs printed by Fig_3_to_6_localization_across_N.m. To ensure reproducibility, set the random number generator (rng) seed to the same value used in the repository. Since the simulations are time-intensive and can take several days to run, the final results have been included in the repository for convenience. 

To generate Fig. 9 (CDF of localization errors):
	1.First, run Fig9_CDF_Localization_error.m to compute the required localization error data.
	2. Then, run Fig9_CDF_Localization_error_plot.m to generate the plot. Be sure to update the variable name and plot title to match the specific localization method being plotted.
	3. Alternatively, to plot directly without re-running the computations, you may import the precomputed file Fig9_CDF_Localization_error_plot.mat, which contains results for all positioning methods.

Although not discussed in the paper, the simulation code for the distributed-mean-z-score method is also included in this repository. This approach applies z-score filtering to eliminate outliers before averaging the position estimates. While it achieves lower localization error compared to the simple mean, the error ellipse area is slightly higher than that of the mean—but still significantly lower than that of the median.

For the FCNN baseline, the referenced paper specifies only the number of layers, without detailing other architectural or training parameters. Therefore, hyperparameters such as the number of nodes per layer, activation functions, optimizer, and number of training epochs were selected based on performance observed across multiple experimental runs for K=225 training points.

In our experimental setup, shadowing effects are modeled using a spatially correlated Gaussian process to reflect realistic signal behavior. The correlation between shadowing coefficients decays exponentially with distance, as specified by the model 
ρ(d)=2^(−d/decorr), where decorr is a scenario-specific decorrelation distance (per 3GPP). A correlation matrix is incrementally built during the offline phase by computing pairwise distances between reference points (RPs), enabling the generation of correlated shadowing terms via conditional Gaussian statistics. For each new RP, the conditional mean and standard deviation are computed based on the existing shadowing values and their spatial correlations, following Theorem 10.2 from Steven Kay’s Estimation Theory (p. 325, a snippet is included in the repo). This process ensures that nearby points are more likely to exhibit similar shadowing, mimicking real-world wireless environments. During the online phase, the same procedure is applied to test points, maintaining consistency in spatial correlation. Implementation details follow the approach outlined in arXiv:2108.02541, particularly page 66, and are available in the accompanying GitHub repository. To better understand and experiment with this model, the MATLAB file shadowing_realisation_example.m illustrates how shadowing is realized in simulations and can be used for hands-on exploration.

