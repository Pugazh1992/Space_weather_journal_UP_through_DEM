# Space_weather_journal_UP_through_DEM — Project Overview

This repository contains MATLAB and Python scripts for the Space Weather journal uncertainty-propagation workflow. The project studies how uncertainty in selected space-weather inputs, primarily F10.7 and Ap, propagates through a global thermosphere-density prediction neural network for Swarm-C density nowcasting.

## Author information

**Author / script owner:** Dr. Pugazhenthi Sivasankar  
**Affiliation:** XBai Research Group, Department of Mechanical and Aerospace Engineering, Rutgers University  
**Research context:** Space-weather input uncertainty modeling, AETHER-P3/global thermosphere-density prediction, Monte-Carlo uncertainty propagation, and Unscented Transform validation.

## Repository purpose

The repository supports a three-stage workflow:

1. **Input generation:** Create MC and UT input `.mat` files by perturbing F10.7 or Ap with a base-10 lognormal uncertainty model.
2. **Model execution:** Run the trained evidential neural network on the MC/UT input files and save density-prediction results.
3. **Post-processing and diagnostics:** Compare MC, UT, and baseline density statistics; validate input distributions; check theoretical-vs-UT moment behavior.

For full file-by-file details, see [`README_detailed_script_documentation.md`](README_detailed_script_documentation.md).

## Main workflow scripts

These files are the main pipeline scripts. They either generate input files or run the neural-network density model.

| File | Role | Main output |
|---|---|---|
| `sd_swarmc_f107_lognormal_mc_ut_v1.m` | MATLAB input-generation script for F10.7 uncertainty. | `Noisy_f107_MC_1000000.mat`, `sigma_points_f107_ut_3.mat`, and histogram figures. |
| `sd_swarmc_ap_lognormal_mc_ut_v1.m` | MATLAB input-generation script for Ap uncertainty. | `Noisy_ap_MC_1000000.mat`, `sigma_points_ap_ut_3.mat`, and histogram figures. |
| `mt_swarmc_f107_lognormal_mc_ut_v1.py` | Python model-execution script for F10.7 MC/UT inputs. | `mc_results.mat` and `ut_results.mat` in each F10.7 case/noise-ratio folder. |
| `mt_swarmc_ap_lognormal_mc_ut_v1.py` | Python model-execution script for Ap MC/UT inputs. | `mc_results.mat` and `ut_results.mat` in each Ap case/noise-ratio folder. |

## Post-processing and diagnostic scripts

These files analyze existing inputs/results. They generally create figures and command-window summaries rather than new `.mat` result files.

| File | Role | Main output |
|---|---|---|
| `pr_swarmc_f107_lognormal_mc_ut_v1.m` | Post-processes F10.7 NN density outputs and compares MC, UT, and baseline physical-density statistics. | 2-by-2 histogram figures and printed MC/UT/baseline summaries; no new `.mat` file by default. |
| `pr_swarmc_ap_lognormal_mc_ut_v1.m` | Post-processes Ap NN density outputs and compares MC, UT, and baseline physical-density statistics. | 2-by-2 histogram figures and printed MC/UT/baseline summaries; no new `.mat` file by default. |
| `train_vs_input_f107_dist.m` | Diagnostic script comparing training F10.7 data, lognormal input samples, UT sigma points, and nonlinear-system outputs. | Three MATLAB figures and command-window statistics tables. |
| `train_vs_input_ap_dist.m` | Diagnostic script comparing training Ap data, lognormal input samples, UT sigma points, and nonlinear-system outputs. | Three MATLAB figures and command-window statistics tables. |
| `th_vs_ut_moments.m` | Self-contained theoretical-vs-UT moment-checking script for linear and nonlinear scalar systems. | Two command-window tables; no figures or saved files by default. |

## Recommended execution order

```text
1. Generate MC/UT input files
   sd_swarmc_f107_lognormal_mc_ut_v1.m
   sd_swarmc_ap_lognormal_mc_ut_v1.m

2. Run the trained density neural network
   mt_swarmc_f107_lognormal_mc_ut_v1.py
   mt_swarmc_ap_lognormal_mc_ut_v1.py

3. Post-process density results
   pr_swarmc_f107_lognormal_mc_ut_v1.m
   pr_swarmc_ap_lognormal_mc_ut_v1.m

4. Run optional diagnostics
   train_vs_input_f107_dist.m
   train_vs_input_ap_dist.m
   th_vs_ut_moments.m
```

## Main data and model folders

```text
Output/Swarm/                          # Swarm-C case folders and generated MC/UT/result files
Output/Model/Model1/                   # GPinput.mat and GPoutput.mat for normalization
Model/                                 # Trained TensorFlow/evidential NN weights
figures/                               # Optional exported figures
```

## Important notes

- Large `.mat` files are usually better kept outside Git or tracked with Git LFS.
- The MATLAB input-generation scripts save MC and UT `.mat` files but do not automatically save their figures.
- The Python model-execution scripts save `mc_results.mat` and `ut_results.mat`.
- The post-processing/diagnostic scripts mainly display figures and print tables/summaries.
- Keep the uncertain input type, noise ratio, MC sample count, UT sigma-point count, and case name visible in output filenames or figure captions.
