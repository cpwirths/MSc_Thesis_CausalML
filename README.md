ReadMe for Master Thesis “Evaluating Causal Machine Learning methods”

I: Atlantic Causal Inference competition simulation study on ATE, GATE and ITE

	1.	Produce the simulation knobs with the script “create_sim_datasets.R” located in the directory “...acic_simulations/data_simulations/”.
	2.	Execute all Causal ML methods in the corresponding subdirectories. Note that, since the computations are implemented on a high-performance computing cluster, each subdirectory contains an example script for one knob (parameter k) which contains 100 simulations replications. To compute more knobs reproduce scripts and change parameter k or run in a loop with batch scripts. Global ML input parameters are stored in “...acic_simulations/globals_parameters.R”. In case you would like to tune ML input parameters, first run “...acic_simulations/k_tune_sim.R” and then one of the general ML methods with "..._tune", e.g. "...dml/k_dml_acic_tune.R".
	3.	After all methods are computed, merge the results for the ATE, GATE and ITE with the “...merge.R” scripts located in the directory “...acic_simulations/output/”. Final output is retrieved by running “...final_output.R”
II: GSS application  

	1.	Download “GSS7218_R1.DTA” from https://gss.norc.org/ and store raw data in directory “.../gss_application/”
	2.	Run “.../gss_application/data_prep.R” to prepare data as input for Causal ML methods
	3.	Run Causal ML method of interest.
References:

	•	Athey, S., Tibshirani, J., & Wager, S. (2019). Generalized random forests. The Annals of Statistics, 47 (2), 1148–1178
	•	Chernozhukov, Victor, et al. "Double/debiased/neyman machine learning of treatment effects." American Economic Review 107.5 (2017): 261-65
	•	Chernozhukov, V., Demirer, M., Duflo, E., & Fernandez-Val, I. (2017). Generic machine learning inference on heterogeneous treatment effects in randomized experiments (Tech. Rep.). Cemmap
	•	Chipman, H. A., George, E. I., McCulloch, R. E., et al. (2010). BART: Bayesian additive regression trees. The Annals of Applied Statistics, 4 (1), 266–298.
	•	Dorie, V., Hill, J., Shalit, U., Scott, M., Cervone, D., et al. (2019). Automated versus do-it-yourself methods for causal inference: Lessons learned from a data analysis competition. Statistical Science, 34 (1), 43–68.
	•	Hill, J., & Su, Y.-S. (2013). Assessing lack of common support in causal inference using Bayesian nonparametrics: Implications for evaluating the effect of breastfeeding on children’s cognitive outcomes. The Annals of Applied Statistics, 1386–1420.
	•	Knaus, Michael, Michael Lechner, and Anthony Strittmatter. "Machine learning estimation of heterogeneous causal effects: Empirical monte carlo evidence." (2018)
