
'''

This script takes the space weather input and runs it through global thermosphere density prediction neural network to
nowcast the thermosphere density. The space weather input with variations in Ap (MC samples and UT sigma points)
are considered with three different cases for Swarm-C satellite. For the MC & UT simulation, Ap samples are obtained from the 
lognormal distribution

Version 1:
- The outputs are collected as MATLAB structures and saved.
- a function is created to execute the NN model to predict the density.
- The key 'xdata1_mc' is changed to 'xdata1_mc_lognormal' in the function compute_density_prediction_results(),
- It is run for only one noise_ratio = [0.1 0.25],
- '/lognormal_mc_ut' is included in child_dir_mc_ap,
- '/lognormal_mc_ut' is included in save_folder_mc,
- '/lognormal_mc_ut' is included in child_dir_ut_ap,
- '/lognormal_mc_ut' is included in save_folder_ut,

'''


# Import the packages:
import numpy as np
import evidential_deep_learning as edl
import tensorflow as tf
import scipy.io as io
from numpy.random import seed
import random
import os
from pathlib import Path
import pandas as pd
import time
import h5py

def loadmat_v73(filename):
    """
    Load MATLAB v7.3 .mat files using h5py
    """
    data = {}
    with h5py.File(filename, 'r') as f:
        for key in f.keys():
            # Get the data
            dataset = f[key]
            if isinstance(dataset, h5py.Dataset):
                data[key] = np.array(dataset)
                # Transpose if it's a matrix (MATLAB uses column-major)
                if data[key].ndim >= 2:
                    data[key] = data[key].T
    return data

def compute_density_prediction_results(filepath_ip, x_mean, x_std, y_mean, y_std, model, uq_method):
    """
    Load an input (.mat v7.3), normalize, run prediction, split into
    (mu, v, alpha, beta), compute epistemic & aleatoric stds, and return
    a dict suitable for saving as a MATLAB struct.

    Parameters
    ----------
    filepath_ip : str or Path
        Path to the MATLAB v7.3 file.
        For MC, expects keys 'xdata1_mc_lognormal' and 'ydata1_mc'.
        For UT, expects keys 'xdata1_ut' and 'ydata1_ut'.
    x_mean, x_std, y_mean, y_std : np.ndarray
        Normalization statistics (vectors) for inputs/targets.
    model : tf.keras.Model
        Trained Keras model that outputs [mu, v, alpha, beta] concatenated.
    uq_method : str
        'MC' or 'UT' — selects which dataset keys to read.

    Returns
    -------
    density_prediction_results : dict
        Dict with fields:
        mu_1, v_1, alpha_1, beta_1, var_1, sigma_1.
    """

    # --- Load data (choose keys based on UQ method) ---
    xt = loadmat_v73(str(filepath_ip))
    if uq_method == 'MC':
        xkey, ykey = 'xdata1_mc_lognormal', 'ydata1_mc'
    else:
        xkey, ykey = 'xdata1_ut', 'ydata1_ut'

    xt_1 = np.real(xt[xkey])
    yt_1 = np.real(xt[ykey])

    # --- Normalization (inputs used for prediction; targets kept for parity) ---
    xt1_nor = (xt_1 - x_mean) / x_std
    yt1_nor = (yt_1 - y_mean) / y_std  # not used later, preserved for completeness

    # --- Prepare result containers (lists -> can be saved as cell arrays in MATLAB if desired) ---
    mu_1, v_1, alpha_1, beta_1 = [], [], [], []
    sigma_1, var_1 = [], []

    # --- Reproducibility & TF session handling (same behavior as your script) ---
    seedvalue = 17
    random.seed(seedvalue)
    tf.random.set_seed(seedvalue)
    print(f'Seed Value = {seedvalue}\n')
    tf.keras.backend.clear_session()

    print(f'Processing input file: {filepath_ip}  |  UQ method: {uq_method}\n')

    # --- Model prediction ---
    cur_pred_1 = model.predict(xt1_nor, verbose=0)

    # Split into (mu, v, alpha, beta)
    cur_mu_1, cur_v_1, cur_alpha_1, cur_beta_1 = tf.split(cur_pred_1, 4, axis=-1)

    # --- Epistemic Std: sqrt(beta / (v * (alpha - 1))) ---
    cur_var_1 = np.sqrt(cur_beta_1 / (cur_v_1 * (cur_alpha_1 - 1)))
    cur_var_1 = np.minimum(cur_var_1, 1e3)  # clip as in your code
    var_1 += [cur_var_1]

    # --- Aleatoric Std: sqrt(beta / (alpha - 1)) ---
    cur_sigma_1 = np.sqrt(cur_beta_1 / (cur_alpha_1 - 1))
    cur_sigma_1 = np.minimum(cur_sigma_1, 1e3)  # clip as in your code
    sigma_1 += [cur_sigma_1]

    # Collect other parameters (keep array nature)
    mu_1    += [np.array(cur_mu_1)]
    v_1     += [cur_v_1]
    alpha_1 += [cur_alpha_1]
    beta_1  += [cur_beta_1]

    # --- Pack into a MATLAB-struct-compatible dict ---
    density_prediction_results = {
        "mu_1":    np.array(mu_1),
        "v_1":     np.array(v_1),
        "alpha_1": np.array(alpha_1),
        "beta_1":  np.array(beta_1),
        "var_1":   np.array(var_1),
        "sigma_1": np.array(sigma_1),
    }

    return density_prediction_results


# set the seed for randomly split the training dataset and validation dataset
randomseed = 9
seed(randomseed)

# create model structure and load saved model
model = tf.keras.models.Sequential()
model.add(tf.keras.Input(shape=(65,)))  
model.add(tf.keras.layers.Dense(512, activation='relu6', kernel_regularizer=tf.keras.regularizers.L2(2e-4)))
model.add(tf.keras.layers.Dropout(0.3))
model.add(tf.keras.layers.Dense(256, activation='relu6', kernel_regularizer=tf.keras.regularizers.L2(2e-4)))
model.add(tf.keras.layers.Dropout(0.3))
model.add(tf.keras.layers.Dense(128, activation='relu6', kernel_regularizer=tf.keras.regularizers.L2(2e-4)))
model.add(tf.keras.layers.Dropout(0.2))
model.add(edl.layers.DenseNormalGamma(1))

model.load_weights('Model/model_Seed_14.weights.h5')


# Build the data as dictionary using information provided by Ruochen.
data = {
    "Report_case":      [1, 2, 3],
    "Satellite":        ["Swarm-C", "Swarm-C", "Swarm-C"],
    "Ap_indication":  ["High_Ap", "Low_Ap", "Medium_Ap"]
}

# Create DataFrame with a suitable name
df_space_weather = pd.DataFrame(data)

# Show the table
print(df_space_weather)


# load dataset from mat file
x = io.loadmat('Output/Model/Model1/GPinput.mat')
y = io.loadmat('Output/Model/Model1/GPoutput.mat')

# training data (For normalization)
x = np.real(x['xdata'])
y = np.real(y['ydata'])

# Obtain the mean and the standard deviation for normalization:
x_mean = np.mean(x, axis=0)
y_mean = np.mean(y, axis=0)
x_std = np.std(x, axis=0)
y_std = np.std(y, axis=0)


# Parent directory string: (constant)
parent_dir_str = 'Output/Swarm/'

#  Set the Gaussian noise for the MC and UT inputs:
# noise_ratio = [0.1, 0.25, 0.5, 0.75, 1]
noise_ratio = [0.1, 0.25]

# Number of Monte-Carlo samples:
n_samples_MC = int(1e6)

# Set the number of sigma points for UT:
n_samples_ut = 3

# String for the noisy data for the current case (MC):
data_str_mc = '/Noisy_ap_MC_' + str(n_samples_MC) + '.mat'

# String for the noisy data for the current case (UT):
data_str_ut = '/sigma_points_ap_ut_' + str(n_samples_ut) + '.mat'

# Start timer
start_time_mc = time.time()

# Loop over the different noise ratios:
for noise_ratio_current in noise_ratio:

    # --- Option B: explicit per-row logic (more like the MATLAB loop) ---
    for _, row in df_space_weather.iterrows():

        # Obtain the Ap test case name:
        rc = row['Ap_indication']

        # Combine the strings to obtain the file path:
        child_dir_mc_ap = rc + '/Noise_ratio_' + str(noise_ratio_current) + '/lognormal_mc_ut' + data_str_mc
        filepath_mc = os.path.join(parent_dir_str, child_dir_mc_ap)

        # Print the filepath_mc for confirmation:
        print(filepath_mc)

        # Save the output into respective files:
        save_folder_mc = parent_dir_str + rc + '/Noise_ratio_' + str(noise_ratio_current) + '/lognormal_mc_ut' + '/'       

        print(save_folder_mc)

        # 2) Create the folder inside the present working directory
        out_dir = Path.cwd() / save_folder_mc
        out_dir.mkdir(parents=True, exist_ok=True)

        # Call the function that uses NN to predict the density:
        mc_results = compute_density_prediction_results(
            filepath_ip = filepath_mc,   # your constructed path
            x_mean = x_mean, x_std = x_std,
            y_mean = y_mean, y_std = y_std,
            model = model, uq_method = 'MC'
        )        

        # Write a single file named 'mc_results.mat' in the case folder
        io.savemat(out_dir / "mc_results.mat", {"mc_results": mc_results})        

        print('Completed the MC execution of test case ' + str(rc) + '\n')
        


# End timer
end_time_mc = time.time()

print(len(df_space_weather))
print(end_time_mc - start_time_mc)

# Divide the total time by the total number of cases:
execution_time_mc = (end_time_mc - start_time_mc)/(len(df_space_weather)*len(noise_ratio))
print(f"Average execution time for {n_samples_MC} MC samples after loading Ruochen's NN: {execution_time_mc:.4f} seconds")


## For processing the UT output files:
# Start timer
start_time_ut = time.time()

# Loop over the different noise ratios:
for noise_ratio_current in noise_ratio:

    # --- Option B: explicit per-row logic (more like the MATLAB loop) ---
    for _, row in df_space_weather.iterrows():

        # Obtain the F107 test case name:
        rc = row['Ap_indication']

        # Combine the strings to obtain the file path:
        child_dir_ut_ap = rc + '/Noise_ratio_' + str(noise_ratio_current) + '/lognormal_mc_ut' + data_str_ut
        filepath_ut = os.path.join(parent_dir_str, child_dir_ut_ap)

        # Print the filepath_mc for confirmation:
        print(filepath_ut)

        # Save the output into respective files:
        save_folder_ut = parent_dir_str + rc + '/Noise_ratio_' + str(noise_ratio_current) + '/lognormal_mc_ut' + '/'       

        print(save_folder_ut)

        # 2) Create the folder inside the present working directory
        out_dir = Path.cwd() / save_folder_ut
        out_dir.mkdir(parents=True, exist_ok=True)

        # Call the function that uses NN to predict the density:
        ut_results = compute_density_prediction_results(
            filepath_ip = filepath_ut,   # your constructed path
            x_mean = x_mean, x_std = x_std,
            y_mean = y_mean, y_std = y_std,
            model = model, uq_method = 'UT'
        )        

        # Write a single file named 'mc_results.mat' in the case folder
        io.savemat(out_dir / "ut_results.mat", {"ut_results": ut_results})        

        print('Completed the UT execution of test case ' + str(rc) + '\n')
        


# End timer
end_time_ut = time.time()

print(len(df_space_weather))
print(end_time_ut - start_time_ut)

# Divide the total time by the total number of cases:
execution_time_ut = (end_time_ut - start_time_ut)/(len(df_space_weather)*len(noise_ratio))
print(f"Average execution time for {n_samples_ut} UT sigma points after loading Ruochen's NN: {execution_time_ut:.4f} seconds")

