# Space_weather_journal_UP_through_DEM — Detailed Script Documentation

This file contains detailed documentation for the MATLAB and Python scripts used in the Swarm-C space-weather input uncertainty propagation workflow. For a shorter project-level overview and script classification, see [`README_project_overview.md`](README_project_overview.md).

## Author information

**Author / script owner:** Dr. Pugazhenthi Sivasankar  
**Affiliation:** XBai Research Group, Department of Mechanical and Aerospace Engineering, Rutgers University  
**Research context:** Space-weather input uncertainty modeling, AETHER-P3/global thermosphere-density prediction, Monte-Carlo uncertainty propagation, and Unscented Transform validation.

## Workflow context

The scripts are organized around three tasks:

1. MATLAB input-generation scripts create MC and UT `.mat` files for uncertain F10.7 or Ap inputs.
2. Python model-execution scripts run the trained evidential neural network on those MC/UT files and save density-prediction result structures.
3. MATLAB post-processing and diagnostic scripts compare MC/UT/baseline density results and validate the lognormal/UT uncertainty formulation.

## Common MATLAB dependencies

The scripts use standard MATLAB functionality such as:

- `load`
- `save`
- `mkdir`
- `table`
- `mean`
- `std`
- `min`
- `max`
- `randn`
- `repmat`
- `histogram`
- `xline`
- `tiledlayout`
- `figure` and other basic plotting utilities

The distribution-comparison scripts also use Statistics and Machine Learning Toolbox functions:

- `skewness`
- `ksdensity`

No Orekit data, Java class, ONNX Runtime call, SP3 file, or orbit-propagation execution is required by these MATLAB scripts.

## Common Python dependencies

The Python model-execution scripts use:

- Python 3
- `numpy`
- `scipy`
- `h5py`
- `pandas`
- `tensorflow`
- `evidential_deep_learning`
- `pathlib`, `os`, `random`, and `time` from the Python standard library

The Python scripts require the trained Keras/TensorFlow evidential model weights and model-normalization files:

- `Model/model_Seed_14.weights.h5`
- `Output/Model/Model1/GPinput.mat`
- `Output/Model/Model1/GPoutput.mat`

They also require the MC/UT `.mat` files produced by the MATLAB input-generation scripts. The input MC/UT files are MATLAB v7.3 files and are loaded using `h5py`; the result files are saved using `scipy.io.savemat`.

## Script summary

| File | Main purpose | Required input files | Output produced |
|---|---|---|---|
| `sd_swarmc_f107_lognormal_mc_ut_v1.m` | Creates MC and UT input datasets by perturbing the F10.7 column with a base-10 lognormal uncertainty model for Swarm-C `High_f107`, `Low_f107`, and `Medium_f107` cases. | `Output/Swarm/High_f107/GPinputS.mat`; `Output/Swarm/High_f107/GPoutputS.mat`; `Output/Swarm/Low_f107/GPinputS.mat`; `Output/Swarm/Low_f107/GPoutputS.mat`; `Output/Swarm/Medium_f107/GPinputS.mat`; `Output/Swarm/Medium_f107/GPoutputS.mat`. | Histogram figures; `Noisy_f107_MC_1000000.mat`; `sigma_points_f107_ut_3.mat`; command-window `tblSpaceWeather`. |
| `sd_swarmc_ap_lognormal_mc_ut_v1.m` | Creates MC and UT input datasets by perturbing the Ap column with a base-10 lognormal uncertainty model for Swarm-C `High_Ap`, `Low_Ap`, and `Medium_Ap` cases. | `Output/Swarm/High_Ap/GPinputS.mat`; `Output/Swarm/High_Ap/GPoutputS.mat`; `Output/Swarm/Low_Ap/GPinputS.mat`; `Output/Swarm/Low_Ap/GPoutputS.mat`; `Output/Swarm/Medium_Ap/GPinputS.mat`; `Output/Swarm/Medium_Ap/GPoutputS.mat`. | Histogram figures; `Noisy_ap_MC_1000000.mat`; `sigma_points_ap_ut_3.mat`; command-window `tblSpaceWeather`. |
| `mt_swarmc_f107_lognormal_mc_ut_v1.py` | Runs the trained evidential global thermosphere-density neural network for Swarm-C F10.7 lognormal MC and UT input cases. | `Model/model_Seed_14.weights.h5`; `Output/Model/Model1/GPinput.mat`; `Output/Model/Model1/GPoutput.mat`; `Noisy_f107_MC_1000000.mat`; `sigma_points_f107_ut_3.mat` in each F10.7 case/noise-ratio folder. | `mc_results.mat`; `ut_results.mat`; printed case table, file paths, and timing diagnostics. No figures are produced. |
| `mt_swarmc_ap_lognormal_mc_ut_v1.py` | Runs the trained evidential global thermosphere-density neural network for Swarm-C Ap lognormal MC and UT input cases. | `Model/model_Seed_14.weights.h5`; `Output/Model/Model1/GPinput.mat`; `Output/Model/Model1/GPoutput.mat`; `Noisy_ap_MC_1000000.mat`; `sigma_points_ap_ut_3.mat` in each Ap case/noise-ratio folder. | `mc_results.mat`; `ut_results.mat`; printed case table, file paths, and timing diagnostics. No figures are produced. |
| `pr_swarmc_f107_lognormal_mc_ut_v1.m` | Post-processes F10.7 lognormal MC/UT neural-network density outputs, converts log-density predictions to physical density, and compares MC, UT, and baseline density statistics. | `Output/Model/Model1/GPoutput.mat`; `Noisy_f107_MC_1000000.mat`; `sigma_points_f107_ut_3.mat`; `mc_results.mat`; `ut_results.mat` in each F10.7 case/noise-ratio folder. | 2-by-2 histogram figures and command-window MC/UT/baseline density summaries. No new `.mat` file is saved by default. |
| `pr_swarmc_ap_lognormal_mc_ut_v1.m` | Post-processes Ap lognormal MC/UT neural-network density outputs, converts log-density predictions to physical density, and compares MC, UT, and baseline density statistics. | `Output/Model/Model1/GPoutput.mat`; `Noisy_ap_MC_1000000.mat`; `sigma_points_ap_ut_3.mat`; `mc_results.mat`; `ut_results.mat` in each Ap case/noise-ratio folder. | 2-by-2 histogram figures and command-window MC/UT/baseline density summaries. No new `.mat` file is saved by default. |
| `train_vs_input_f107_dist.m` | Compares NN/AETHER-P3 training-data F10.7 statistics against synthetic lognormal F10.7 uncertain inputs; compares MC samples with UT sigma points; propagates the uncertain input through a nonlinear scalar test system. | `GPinput.mat`; `High_f107/GPinputS.mat`; `High_f107/GPoutputS.mat`; `Low_f107/GPinputS.mat`; `Low_f107/GPoutputS.mat`; `Medium_f107/GPinputS.mat`; `Medium_f107/GPoutputS.mat`. | Three MATLAB figures and several command-window tables. No `.mat`, `.fig`, `.png`, or `.csv` file is saved by default. |
| `train_vs_input_ap_dist.m` | Performs the analogous distribution, MC-vs-UT, and nonlinear-output diagnostics for the Ap geomagnetic index. | `GPinput.mat`; `High_Ap/GPinputS.mat`; `High_Ap/GPoutputS.mat`; `Low_Ap/GPinputS.mat`; `Low_Ap/GPoutputS.mat`; `Medium_Ap/GPinputS.mat`; `Medium_Ap/GPoutputS.mat`. | Three MATLAB figures and several command-window tables. No `.mat`, `.fig`, `.png`, or `.csv` file is saved by default. |
| `th_vs_ut_moments.m` | Compares theoretical mean/variance/skewness propagation with UT-reconstructed moments for a base-10 lognormal variable passed through linear and nonlinear scalar systems. | No external input file is required. Inputs are defined inside the script. | Two command-window tables. No figures or saved files are produced by default. |

## Detailed script documentation

### 1. `sd_swarmc_f107_lognormal_mc_ut_v1.m`

#### Purpose

This MATLAB script prepares uncertain input data for MC and UT thermosphere-density prediction studies using a Swarm-C test case. It treats F10.7 as the single uncertain space-weather input while keeping the other model features deterministic.

The script considers:

- `High_f107`
- `Low_f107`
- `Medium_f107`

For each case and active noise ratio, the script:

1. Loads the deterministic Swarm-C input vector and output density value.
2. Creates an MC ensemble by repeating the deterministic input vector.
3. Replaces only column 12, the F10.7 column, with base-10 lognormal MC samples.
4. Creates three UT sigma points for the same base-10 lognormal F10.7 uncertainty model.
5. Replaces column 12 in the UT input matrix with the three UT sigma points.
6. Displays a histogram comparing the MC F10.7 distribution and UT sigma-point locations.
7. Saves the MC and UT input `.mat` files for downstream global thermosphere-density prediction runs.

#### External file dependencies

MATLAB dependencies:

- MATLAB base functionality: `load`, `table`, `repmat`, `randn`, `figure`, `tiledlayout`, `histogram`, `xline`, `mkdir`, `save`, and basic array operations.

Project/data dependencies:

- Parent directory: `Output/Swarm/`
- Swarm-C case folders:
  - `High_f107/`
  - `Low_f107/`
  - `Medium_f107/`

#### Input files required

Required files:

- `Output/Swarm/High_f107/GPinputS.mat`
- `Output/Swarm/High_f107/GPoutputS.mat`
- `Output/Swarm/Low_f107/GPinputS.mat`
- `Output/Swarm/Low_f107/GPoutputS.mat`
- `Output/Swarm/Medium_f107/GPinputS.mat`
- `Output/Swarm/Medium_f107/GPoutputS.mat`

Required variables:

- `xdata` in each `GPinputS.mat` file
- `ydata` in each `GPoutputS.mat` file

The script assumes this `xdata` column convention:

| Column | Quantity |
|---:|---|
| 7 | Satellite latitude |
| 8 | Satellite longitude |
| 9 | Satellite altitude |
| 12 | F10.7 |
| 15 | Dst |
| 16 | Ap |

#### Output produced: figures, tables, and `.mat` files

Figures:

- One histogram figure per case and noise ratio.
- With the active setup, the script generates 6 figures: 3 F10.7 cases × 2 noise ratios.
- Figures are displayed but not automatically saved.

Command-window output:

- `tblSpaceWeather`, containing case number, satellite name, selected space-weather values, satellite location values, and noise-ratio information.

Saved output directory:

```text
Output/Swarm/<case>/Noise_ratio_<noise_ratio>/lognormal_mc_ut/
```

Saved files:

```text
Noisy_f107_MC_1000000.mat
sigma_points_f107_ut_3.mat
```

The MC file contains:

- `xdata1_mc_lognormal`
- `ydata1_mc`
- `f107_selected`
- `Dst_selected`
- `Ap_selected`
- `noise_ratio_current`

The UT file contains:

- `xdata1_ut`
- `ydata1_ut`
- `f107_selected`
- `Dst_selected`
- `Ap_selected`
- `noise_ratio_current`
- `Wm`
- `Wc`

#### Author information

Author / script owner: Dr. Pugazhenthi Sivasankar  
Research context: Space-weather input uncertainty propagation through a global thermosphere-density prediction framework using MC and UT.

---

### 2. `sd_swarmc_ap_lognormal_mc_ut_v1.m`

#### Purpose

This MATLAB script prepares uncertain input data for MC and UT thermosphere-density prediction studies using a Swarm-C test case. It treats Ap as the single uncertain space-weather input while keeping the other model features deterministic.

The script considers:

- `High_Ap`
- `Low_Ap`
- `Medium_Ap`

For each case and active noise ratio, the script:

1. Loads the deterministic Swarm-C input vector and output density value.
2. Creates an MC ensemble by repeating the deterministic input vector.
3. Replaces only column 16, the Ap column, with base-10 lognormal MC samples.
4. Creates three UT sigma points for the same base-10 lognormal Ap uncertainty model.
5. Replaces column 16 in the UT input matrix with the three UT sigma points.
6. Displays a histogram comparing the MC Ap distribution and UT sigma-point locations.
7. Saves the MC and UT input `.mat` files for downstream global thermosphere-density prediction runs.

#### External file dependencies

MATLAB dependencies:

- MATLAB base functionality: `load`, `table`, `repmat`, `randn`, `figure`, `tiledlayout`, `histogram`, `xline`, `mkdir`, `save`, and basic array operations.

Project/data dependencies:

- Parent directory: `Output/Swarm/`
- Swarm-C case folders:
  - `High_Ap/`
  - `Low_Ap/`
  - `Medium_Ap/`

#### Input files required

Required files:

- `Output/Swarm/High_Ap/GPinputS.mat`
- `Output/Swarm/High_Ap/GPoutputS.mat`
- `Output/Swarm/Low_Ap/GPinputS.mat`
- `Output/Swarm/Low_Ap/GPoutputS.mat`
- `Output/Swarm/Medium_Ap/GPinputS.mat`
- `Output/Swarm/Medium_Ap/GPoutputS.mat`

Required variables:

- `xdata` in each `GPinputS.mat` file
- `ydata` in each `GPoutputS.mat` file

The script assumes this `xdata` column convention:

| Column | Quantity |
|---:|---|
| 7 | Satellite latitude |
| 8 | Satellite longitude |
| 9 | Satellite altitude |
| 12 | F10.7 |
| 15 | Dst |
| 16 | Ap |

#### Output produced: figures, tables, and `.mat` files

Figures:

- One histogram figure per case and noise ratio.
- With the active setup, the script generates 6 figures: 3 Ap cases × 2 noise ratios.
- Figures are displayed but not automatically saved.

Command-window output:

- `tblSpaceWeather`, containing case number, satellite name, selected space-weather values, satellite location values, and noise-ratio information.

Saved output directory:

```text
Output/Swarm/<case>/Noise_ratio_<noise_ratio>/lognormal_mc_ut/
```

Saved files:

```text
Noisy_ap_MC_1000000.mat
sigma_points_ap_ut_3.mat
```

The MC file contains:

- `xdata1_mc_lognormal`
- `ydata1_mc`
- `f107_selected`
- `Dst_selected`
- `Ap_selected`
- `noise_ratio_current`

The UT file contains:

- `xdata1_ut`
- `ydata1_ut`
- `f107_selected`
- `Dst_selected`
- `Ap_selected`
- `noise_ratio_current`
- `Wm`
- `Wc`

#### Author information

Author / script owner: Dr. Pugazhenthi Sivasankar  
Research context: Space-weather input uncertainty propagation through a global thermosphere-density prediction framework using MC and UT.

---

### 3. `mt_swarmc_f107_lognormal_mc_ut_v1.py`

#### Purpose

This Python script performs the density-nowcasting/model-execution stage for Swarm-C cases in which F10.7 is the uncertain space-weather input. It loads the lognormal MC input samples and UT sigma-point input files generated by `sd_swarmc_f107_lognormal_mc_ut_v1.m`, normalizes the 65-dimensional inputs using the training-data statistics, runs the trained evidential neural network, and saves density-prediction and uncertainty outputs as MATLAB-readable `.mat` files.

The active cases are:

- `High_f107`
- `Low_f107`
- `Medium_f107`

The active noise ratios are:

- `0.1`
- `0.25`

#### External file dependencies

Python/package dependencies:

- Python 3
- `numpy`
- `scipy`
- `h5py`
- `pandas`
- `tensorflow`
- `evidential_deep_learning`

Model and normalization dependencies:

- `Model/model_Seed_14.weights.h5`
- `Output/Model/Model1/GPinput.mat`
- `Output/Model/Model1/GPoutput.mat`

Upstream generated input directories:

- `Output/Swarm/High_f107/Noise_ratio_<noise_ratio>/lognormal_mc_ut/`
- `Output/Swarm/Low_f107/Noise_ratio_<noise_ratio>/lognormal_mc_ut/`
- `Output/Swarm/Medium_f107/Noise_ratio_<noise_ratio>/lognormal_mc_ut/`

#### Input files required

For each F10.7 case and noise ratio, the script expects:

- `Noisy_f107_MC_1000000.mat`
- `sigma_points_f107_ut_3.mat`

The MC input file should contain `xdata1_mc_lognormal` and `ydata1_mc`. The UT input file should contain `xdata1_ut`, `ydata1_ut`, `Wm`, and `Wc`.

#### Output produced: figures, tables, and `.mat` files

Figures produced:

- None. This script does not create figures.

Command-window output:

- Swarm-C F10.7 case table.
- Input and output folder paths for each case.
- Completion messages for MC and UT model evaluations.
- Average execution times for MC and UT inference.

Saved `.mat` files:

For each F10.7 case and noise ratio, the script saves:

- `mc_results.mat`
- `ut_results.mat`

inside:

```text
Output/Swarm/<High_f107 | Low_f107 | Medium_f107>/Noise_ratio_<0.1 | 0.25>/lognormal_mc_ut/
```

The saved structures contain:

- `mu_1`
- `v_1`
- `alpha_1`
- `beta_1`
- `var_1`
- `sigma_1`

Here, `var_1` is the clipped epistemic standard-deviation estimate and `sigma_1` is the clipped aleatoric standard-deviation estimate from the evidential model output.

#### Author information

Author / script owner: Dr. Pugazhenthi Sivasankar  
Research context: F10.7 input uncertainty propagation through an evidential global thermosphere-density neural network for Swarm-C density nowcasting.

### 4. `mt_swarmc_ap_lognormal_mc_ut_v1.py`

#### Purpose

This Python script performs the density-nowcasting/model-execution stage for Swarm-C cases in which Ap is the uncertain space-weather input. It loads the lognormal MC input samples and UT sigma-point input files generated by `sd_swarmc_ap_lognormal_mc_ut_v1.m`, normalizes the 65-dimensional inputs using the training-data statistics, runs the trained evidential neural network, and saves density-prediction and uncertainty outputs as MATLAB-readable `.mat` files.

The active cases are:

- `High_Ap`
- `Low_Ap`
- `Medium_Ap`

The active noise ratios are:

- `0.1`
- `0.25`

#### External file dependencies

Python/package dependencies:

- Python 3
- `numpy`
- `scipy`
- `h5py`
- `pandas`
- `tensorflow`
- `evidential_deep_learning`

Model and normalization dependencies:

- `Model/model_Seed_14.weights.h5`
- `Output/Model/Model1/GPinput.mat`
- `Output/Model/Model1/GPoutput.mat`

Upstream generated input directories:

- `Output/Swarm/High_Ap/Noise_ratio_<noise_ratio>/lognormal_mc_ut/`
- `Output/Swarm/Low_Ap/Noise_ratio_<noise_ratio>/lognormal_mc_ut/`
- `Output/Swarm/Medium_Ap/Noise_ratio_<noise_ratio>/lognormal_mc_ut/`

#### Input files required

For each Ap case and noise ratio, the script expects:

- `Noisy_ap_MC_1000000.mat`
- `sigma_points_ap_ut_3.mat`

The MC input file should contain `xdata1_mc_lognormal` and `ydata1_mc`. The UT input file should contain `xdata1_ut`, `ydata1_ut`, `Wm`, and `Wc`.

#### Output produced: figures, tables, and `.mat` files

Figures produced:

- None. This script does not create figures.

Command-window output:

- Swarm-C Ap case table.
- Input and output folder paths for each case.
- Completion messages for MC and UT model evaluations.
- Average execution times for MC and UT inference.

Saved `.mat` files:

For each Ap case and noise ratio, the script saves:

- `mc_results.mat`
- `ut_results.mat`

inside:

```text
Output/Swarm/<High_Ap | Low_Ap | Medium_Ap>/Noise_ratio_<0.1 | 0.25>/lognormal_mc_ut/
```

The saved structures contain:

- `mu_1`
- `v_1`
- `alpha_1`
- `beta_1`
- `var_1`
- `sigma_1`

Here, `var_1` is the clipped epistemic standard-deviation estimate and `sigma_1` is the clipped aleatoric standard-deviation estimate from the evidential model output.

#### Author information

Author / script owner: Dr. Pugazhenthi Sivasankar  
Research context: Ap input uncertainty propagation through an evidential global thermosphere-density neural network for Swarm-C density nowcasting.

### 5. `train_vs_input_f107_dist.m`

#### Purpose

This MATLAB diagnostic script compares the distribution of the F10.7 solar-radio flux index in the NN/AETHER-P3 training input data against a base-10 lognormal uncertain-input distribution generated for Monte-Carlo-style space-weather uncertainty studies.

The script considers three Swarm-C representative cases:

1. `High_f107`
2. `Low_f107`
3. `Medium_f107`

For each case, the script:

- loads the deterministic Swarm-C input/output sample,
- perturbs only the F10.7 input using a base-10 lognormal distribution,
- generates a large Monte-Carlo sample set,
- constructs Unscented Transform sigma points for the same lognormal input,
- compares MC and UT input statistics,
- propagates the uncertain F10.7 input through the nonlinear test system `y = 5*x^2 + 3`,
- compares output statistics from MC samples and UT recombination.

#### External file dependencies

MATLAB dependencies:

- MATLAB base functionality
- Statistics and Machine Learning Toolbox:
  - `skewness`
  - `ksdensity`

Project data dependencies:

- `GPinput.mat`
- Swarm-C case files under `parent_dir_str`:
  - `High_f107/GPinputS.mat`
  - `High_f107/GPoutputS.mat`
  - `Low_f107/GPinputS.mat`
  - `Low_f107/GPoutputS.mat`
  - `Medium_f107/GPinputS.mat`
  - `Medium_f107/GPoutputS.mat`

#### Input files required

The following files must be available before running the script:

1. `GPinput.mat` in the active MATLAB path or current working directory.
2. The case-specific Swarm-C files listed above inside the directory defined by `parent_dir_str`.

Required variables:

- `xdata` in `GPinput.mat` and each `GPinputS.mat` file
- `ydata` in each `GPoutputS.mat` file

#### Output produced: figures, tables, and files

Figures produced:

1. Comparative histograms of NN training F10.7 data and uncertain lognormal F10.7 input samples for the three Swarm-C cases.
2. Comparison of Monte-Carlo F10.7 input histograms against the UT sigma points.
3. Comparison of nonlinear-system outputs generated from MC F10.7 samples and UT sigma-point outputs.

Command-window tables/diagnostics:

- `train_SW_stats`
- `tblSpaceWeather`
- `input_SW_stats`
- `ut_SW_stats`
- `output_mc_stats`
- `output_ut_stats`
- `ovl_coeff`

Saved output files:

- None by default. The script displays figures and prints tables, but it does not save `.mat`, `.fig`, `.png`, `.csv`, or other output files unless explicit save commands are added later.

#### Author information

Author / script owner: Dr. Pugazhenthi Sivasankar  
Research context: Space-weather input uncertainty modeling and uncertainty propagation for AETHER-P3 thermospheric-density/orbit-prediction studies.

---

### 6. `train_vs_input_ap_dist.m`

#### Purpose

This MATLAB diagnostic script compares the distribution of the Ap geomagnetic index in the NN/AETHER-P3 training input data against a base-10 lognormal uncertain-input distribution generated for Monte-Carlo-style space-weather uncertainty studies.

The script considers three Swarm-C representative cases:

1. `High_Ap`
2. `Low_Ap`
3. `Medium_Ap`

For each case, the script:

- loads the deterministic Swarm-C input/output sample,
- perturbs only the Ap input using a base-10 lognormal distribution,
- generates a large Monte-Carlo sample set,
- constructs Unscented Transform sigma points for the same lognormal input,
- compares MC and UT input statistics,
- propagates the uncertain Ap input through the nonlinear test system `y = 5*x^2 + 3`,
- compares output statistics from MC samples and UT recombination.

#### External file dependencies

MATLAB dependencies:

- MATLAB base functionality
- Statistics and Machine Learning Toolbox:
  - `skewness`
  - `ksdensity`

Project data dependencies:

- `GPinput.mat`
- Swarm-C case files under `parent_dir_str`:
  - `High_Ap/GPinputS.mat`
  - `High_Ap/GPoutputS.mat`
  - `Low_Ap/GPinputS.mat`
  - `Low_Ap/GPoutputS.mat`
  - `Medium_Ap/GPinputS.mat`
  - `Medium_Ap/GPoutputS.mat`

#### Input files required

The following files must be available before running the script:

1. `GPinput.mat` in the active MATLAB path or current working directory.
2. The case-specific Swarm-C files listed above inside the directory defined by `parent_dir_str`.

Required variables:

- `xdata` in `GPinput.mat` and each `GPinputS.mat` file
- `ydata` in each `GPoutputS.mat` file

#### Output produced: figures, tables, and files

Figures produced:

1. Comparative histograms of NN training Ap data and uncertain lognormal Ap input samples for the three Swarm-C cases.
2. Comparison of Monte-Carlo Ap input histograms against the UT sigma points.
3. Comparison of nonlinear-system outputs generated from MC Ap samples and UT sigma-point outputs.

Command-window tables/diagnostics:

- `train_SW_stats`
- `tblSpaceWeather`
- `input_SW_stats`
- `ut_SW_stats`
- `output_mc_stats`
- `output_ut_stats`
- `ovl_coeff`

Saved output files:

- None by default. The script displays figures and prints tables, but it does not save `.mat`, `.fig`, `.png`, `.csv`, or other output files unless explicit save commands are added later.

#### Author information

Author / script owner: Dr. Pugazhenthi Sivasankar  
Research context: Space-weather input uncertainty modeling and uncertainty propagation for AETHER-P3 thermospheric-density/orbit-prediction studies.

---

### 7. `th_vs_ut_moments.m`

#### Purpose

This MATLAB verification script compares theoretical moment propagation with moments reconstructed using Unscented Transform recombination equations.

The input random variable is represented as a base-10 lognormal variable. The script computes the theoretical mean, variance, and skewness of the input and compares them with the corresponding UT-reconstructed moments obtained from log-domain sigma points mapped back to the linear domain.

The same comparison is repeated after propagating the input through two scalar systems:

1. Linear system: `y = 5*x + 3`
2. Nonlinear system: `y = 5*x^2 + 3`

This script is intended as a compact mathematical consistency check for the lognormal UT construction used in the larger space-weather and thermospheric-density uncertainty-propagation workflow.

#### External file dependencies

MATLAB dependencies:

- MATLAB base functionality
- No external MATLAB toolboxes are explicitly required by the current script.

Project data/model dependencies:

- None. This script is self-contained and does not load `GPinput.mat`, `GPoutput.mat`, Orekit `.mat` files, ONNX model files, SP3 files, or CHAMP/Swarm density files.

#### Input files required

No external input file is required.

The input distribution and system coefficients are defined directly in the script:

| Variable | Value |
|---|---:|
| `target_mean` | 100 |
| `target_variance` | 4 |
| `sys_coeff_a` | 5 |
| `sys_coeff_b` | 3 |

The local function `ut_sigma_points_log10_lognormal` constructs the UT sigma points and weights from the specified linear-domain mean and variance.

#### Output produced: figures, tables, and files

Figures produced:

- None. This script does not generate figures.

Command-window tables/diagnostics:

- `input_table`: theoretical and UT-reconstructed mean, variance, and skewness of the base-10 lognormal input random variable.
- `output_table`: theoretical and UT-reconstructed mean, variance, and skewness after propagation through the linear and nonlinear scalar systems.

Saved output files:

- None by default. The script prints MATLAB tables to the Command Window, but it does not save `.mat`, `.fig`, `.png`, `.csv`, or other output files unless explicit save commands are added later.

#### Author information

Author / script owner: Dr. Pugazhenthi Sivasankar  
Research context: Verification of Unscented Transform moment propagation for base-10 lognormal uncertainty variables used in space-weather and orbit-uncertainty studies.

### 8. `pr_swarmc_f107_lognormal_mc_ut_v1.m`

#### Purpose

This MATLAB post-processing script analyzes the neural-network density-prediction results for the Swarm-C F10.7 lognormal uncertainty cases. It treats F10.7 as the only uncertain input among the 65 neural-network features and compares the resulting MC, UT, and baseline density predictions.

The script uses the MC/UT input files generated by `sd_swarmc_f107_lognormal_mc_ut_v1.m` and the neural-network output files generated by `mt_swarmc_f107_lognormal_mc_ut_v1.py`. It converts the normalized/logarithmic density outputs back to physical density space, combines aleatoric and epistemic uncertainty, forms consolidated MC density statistics, recombines the UT density outputs, and prints MC-vs-UT relative differences.

The active cases are:

- `High_f107`
- `Low_f107`
- `Medium_f107`

The active noise-ratio setting in the uploaded script is `noise_ratio = 0.25`.

#### External file dependencies

MATLAB/toolbox dependencies:

- MATLAB base functionality for loading `.mat` files, array operations, tables, command-window output, and plotting.
- Statistics and Machine Learning Toolbox functions, including `makedist` and `lognrnd`.
- Parallel Computing Toolbox is recommended because the conversion from lognormal parameters to physical density uses `parfor`.

Project/data dependencies:

- `Output/Model/Model1/GPoutput.mat`
- `Output/Swarm/<F10.7 case>/Noise_ratio_<ratio>/lognormal_mc_ut/`
- MC/UT input `.mat` files from the F10.7 input-generation script.
- MC/UT neural-network result `.mat` files from the F10.7 Python model-execution script.

#### Input files required

Normalization/reference file:

- `Output/Model/Model1/GPoutput.mat`

For each F10.7 case and active noise ratio, the script expects the following files inside `Output/Swarm/<case>/Noise_ratio_<noise_ratio>/lognormal_mc_ut/`:

- `Noisy_f107_MC_1000000.mat`
- `sigma_points_f107_ut_3.mat`
- `mc_results.mat`
- `ut_results.mat`

The supported case folders are `High_f107`, `Low_f107`, and `Medium_f107`.

#### Output produced: figures, tables, and `.mat` files

Figures produced:

- One MATLAB figure per active case/noise-ratio combination, named `MC Results Histograms (2x3)`, using a 2-by-2 layout:
  1. Input lognormal F10.7 MC distribution with UT sigma-point locations.
  2. Histogram of predicted normalized/logarithmic density mean, `mu_hat_z`.
  3. Histogram of predicted normalized/logarithmic density standard deviation, `sigma_hat_z`.
  4. Consolidated physical-density distribution comparing MC, UT, and baseline mean/one-sigma intervals.

Tables/printed outputs:

- Command-window display of `tblSpaceWeather`.
- For each active case/noise-ratio combination, printed MC, baseline, and UT mean/std density values in physical space.
- Printed relative errors between UT and MC density mean and standard deviation.

`.mat` files produced:

- No new `.mat` file is saved by default. The script loads existing input/result `.mat` files and produces figures plus command-window diagnostics.

#### Author information

**Author / script owner:** Dr. Pugazhenthi Sivasankar  
**Affiliation:** XBai Research Group, Department of Mechanical and Aerospace Engineering, Rutgers University  
**Research context:** Propagation of lognormal F10.7 input uncertainty through a global thermosphere-density neural network using MC and UT methods.

### 9. `pr_swarmc_ap_lognormal_mc_ut_v1.m`

#### Purpose

This MATLAB post-processing script analyzes the neural-network density-prediction results for the Swarm-C Ap lognormal uncertainty cases. It treats Ap as the only uncertain input among the 65 neural-network features and compares the resulting MC, UT, and baseline density predictions.

The script uses the MC/UT input files generated by `sd_swarmc_ap_lognormal_mc_ut_v1.m` and the neural-network output files generated by `mt_swarmc_ap_lognormal_mc_ut_v1.py`. It converts the normalized/logarithmic density outputs back to physical density space, combines aleatoric and epistemic uncertainty, forms consolidated MC density statistics, recombines the UT density outputs, and prints MC-vs-UT relative differences.

The active cases are:

- `High_Ap`
- `Low_Ap`
- `Medium_Ap`

The active noise-ratio setting in the uploaded script is `noise_ratio = 0.25`.

#### External file dependencies

MATLAB/toolbox dependencies:

- MATLAB base functionality for loading `.mat` files, array operations, tables, command-window output, and plotting.
- Statistics and Machine Learning Toolbox functions, including `makedist` and `lognrnd`.
- Parallel Computing Toolbox is recommended because the conversion from lognormal parameters to physical density uses `parfor`.

Project/data dependencies:

- `Output/Model/Model1/GPoutput.mat`
- `Output/Swarm/<Ap case>/Noise_ratio_<ratio>/lognormal_mc_ut/`
- MC/UT input `.mat` files from the Ap input-generation script.
- MC/UT neural-network result `.mat` files from the Ap Python model-execution script.

#### Input files required

Normalization/reference file:

- `Output/Model/Model1/GPoutput.mat`

For each Ap case and active noise ratio, the script expects the following files inside `Output/Swarm/<case>/Noise_ratio_<noise_ratio>/lognormal_mc_ut/`:

- `Noisy_ap_MC_1000000.mat`
- `sigma_points_ap_ut_3.mat`
- `mc_results.mat`
- `ut_results.mat`

The supported case folders are `High_Ap`, `Low_Ap`, and `Medium_Ap`.

#### Output produced: figures, tables, and `.mat` files

Figures produced:

- One MATLAB figure per active case/noise-ratio combination, named `MC Results Histograms (2x3)`, using a 2-by-2 layout:
  1. Input lognormal Ap MC distribution with UT sigma-point locations.
  2. Histogram of predicted normalized/logarithmic density mean, `mu_hat_z`.
  3. Histogram of predicted normalized/logarithmic density standard deviation, `sigma_hat_z`.
  4. Consolidated physical-density distribution comparing MC, UT, and baseline mean/one-sigma intervals.

Tables/printed outputs:

- Command-window display of `tblSpaceWeather`.
- For each active case/noise-ratio combination, printed MC, baseline, and UT mean/std density values in physical space.
- Printed relative errors between UT and MC density mean and standard deviation.

`.mat` files produced:

- No new `.mat` file is saved by default. The script loads existing input/result `.mat` files and produces figures plus command-window diagnostics.

#### Author information

**Author / script owner:** Dr. Pugazhenthi Sivasankar  
**Affiliation:** XBai Research Group, Department of Mechanical and Aerospace Engineering, Rutgers University  
**Research context:** Propagation of lognormal Ap input uncertainty through a global thermosphere-density neural network using MC and UT methods.

## Recommended repository organization

A clean repository structure is:

```text
Space_weather_journal_UP_through_DEM/
├── README.md
├── sd_swarmc_f107_lognormal_mc_ut_v1.m
├── sd_swarmc_ap_lognormal_mc_ut_v1.m
├── mt_swarmc_f107_lognormal_mc_ut_v1.py
├── mt_swarmc_ap_lognormal_mc_ut_v1.py
├── train_vs_input_f107_dist.m
├── train_vs_input_ap_dist.m
├── th_vs_ut_moments.m
├── GPinput.mat                         # local/project data; usually not tracked if large
├── Model/
│   └── model_Seed_14.weights.h5
├── Output/
│   ├── Model/
│   │   └── Model1/
│   │       ├── GPinput.mat
│   │       └── GPoutput.mat
│   └── Swarm/
│       ├── High_f107/
│       │   ├── GPinputS.mat
│       │   ├── GPoutputS.mat
│       │   └── Noise_ratio_<value>/
│       │       └── lognormal_mc_ut/
│       │           ├── Noisy_f107_MC_1000000.mat
│       │           ├── sigma_points_f107_ut_3.mat
│       │           ├── mc_results.mat
│       │           └── ut_results.mat
│       ├── Low_f107/
│       │   ├── GPinputS.mat
│       │   └── GPoutputS.mat
│       ├── Medium_f107/
│       │   ├── GPinputS.mat
│       │   └── GPoutputS.mat
│       ├── High_Ap/
│       │   ├── GPinputS.mat
│       │   ├── GPoutputS.mat
│       │   └── Noise_ratio_<value>/
│       │       └── lognormal_mc_ut/
│       │           ├── Noisy_ap_MC_1000000.mat
│       │           ├── sigma_points_ap_ut_3.mat
│       │           ├── mc_results.mat
│       │           └── ut_results.mat
│       ├── Low_Ap/
│       │   ├── GPinputS.mat
│       │   └── GPoutputS.mat
│       └── Medium_Ap/
│           ├── GPinputS.mat
│           └── GPoutputS.mat
└── figures/                            # optional exported figures
```

Large `.mat` files are often better kept outside Git or tracked with Git LFS.

## Example MATLAB usage

```matlab
% From the repository or script folder:
run('sd_swarmc_f107_lognormal_mc_ut_v1.m')
run('sd_swarmc_ap_lognormal_mc_ut_v1.m')
run('train_vs_input_f107_dist.m')
run('train_vs_input_ap_dist.m')
run('th_vs_ut_moments.m')
```

Example Python workflow after the MATLAB input-generation files have been created:

```bash
python mt_swarmc_f107_lognormal_mc_ut_v1.py
python mt_swarmc_ap_lognormal_mc_ut_v1.py
```

Before running scripts that use Swarm-C case files, check that `parent_dir_str` points to the directory containing the required folders.

## Reproducibility notes

- Keep the noise ratio, MC sample count, UT sigma-point count, and case names visible in the script or figure captions.
- Record whether the uncertain input is F10.7 or Ap.
- Clearly distinguish MC-sampled input files from UT sigma-point input files.
- Clearly distinguish MC-sampled moments from UT-reconstructed moments.
- When exporting figures for the manuscript, include the script name, input index, noise ratio, and case name in the file name.
- The input-generation scripts save MC and UT `.mat` files, but their histogram figures are displayed only unless explicit export commands are added.
- The Python model-execution scripts save `mc_results.mat` and `ut_results.mat` for each case/noise-ratio folder.
- The diagnostic scripts display figures and/or tables but do not save them automatically.
