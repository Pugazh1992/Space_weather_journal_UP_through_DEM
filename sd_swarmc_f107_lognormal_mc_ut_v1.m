%{
===============================================================================
File name
===============================================================================
sd_swarmc_f107_lognormal_mc_ut_v1.m

===============================================================================
Purpose
===============================================================================
This MATLAB script prepares uncertain input data for Monte-Carlo (MC) and
Unscented Transform (UT) thermosphere-density prediction studies using a
Swarm-C test case.  The script treats the F10.7 solar radio flux index as the
single uncertain space-weather input and keeps the remaining model input
features deterministic.

For each selected Swarm-C case and noise ratio, the script:
  1. Loads the deterministic input vector and output density value.
  2. Repeats the deterministic input vector to form an MC input ensemble.
  3. Replaces the F10.7 column with base-10 lognormal MC samples.
  4. Constructs three UT sigma points for the same base-10 lognormal F10.7
     uncertainty model.
  5. Replaces the F10.7 column in the UT input matrix with these sigma points.
  6. Generates a histogram showing the MC F10.7 distribution and the UT sigma
     point locations.
  7. Saves the MC and UT input .mat files for downstream global thermosphere
     density-prediction runs.

The active cases are:
  - High_f107
  - Low_f107
  - Medium_f107

The active noise ratios are:
  - 0.10
  - 0.25

The active MC sample count is 1,000,000, and the UT uses 3 sigma points because
only one scalar random variable, F10.7, is being perturbed.

===============================================================================
External file dependencies
===============================================================================
MATLAB/toolbox dependencies:
  - MATLAB base functionality: load, table, repmat, randn, figure,
    tiledlayout, histogram, xline, mkdir, save, and basic array operations.
  - No Orekit, Java, ONNX Runtime, or SP3 reader is called by this script.

Project/data dependencies:
  - The parent data directory is set by:
        parent_dir_str = 'Output/Swarm/';
  - The following Swarm-C case folders must exist under parent_dir_str:
        High_f107/
        Low_f107/
        Medium_f107/

===============================================================================
Input files required
===============================================================================
For each case folder, the script expects:
  - GPinputS.mat
  - GPoutputS.mat

Therefore, the required files are:
  - Output/Swarm/High_f107/GPinputS.mat
  - Output/Swarm/High_f107/GPoutputS.mat
  - Output/Swarm/Low_f107/GPinputS.mat
  - Output/Swarm/Low_f107/GPoutputS.mat
  - Output/Swarm/Medium_f107/GPinputS.mat
  - Output/Swarm/Medium_f107/GPoutputS.mat

Required variables inside the .mat files:
  - xdata in each GPinputS.mat file
  - ydata in each GPoutputS.mat file

The script assumes the following xdata column convention:
  - Column 7  : satellite latitude
  - Column 8  : satellite longitude
  - Column 9  : satellite altitude
  - Column 12 : F10.7 solar radio flux index
  - Column 15 : Dst geomagnetic index
  - Column 16 : Ap geomagnetic index

===============================================================================
Output produced: figures, tables, and .mat files
===============================================================================
Figures:
  - One MATLAB histogram figure is generated for each case and noise ratio.
  - Each figure shows the base-10 lognormal MC F10.7 input distribution and
    overlays the three UT sigma-point locations.
  - With the active setup, the script generates 6 figures:
        3 F10.7 cases x 2 noise ratios.
  - The figures are displayed but are not automatically saved as .fig, .png,
    or .pdf files.

Command-window output:
  - The table tblSpaceWeather is displayed after the loop.  It contains the
    case number, satellite name, selected F10.7/Dst/Ap values, satellite
    location information, and the active noise-ratio value stored during the
    loop.

Saved .mat files:
  For each case and noise ratio, output files are saved in:
        Output/Swarm/<case>/Noise_ratio_<noise_ratio>/lognormal_mc_ut/

  MC input file:
        Noisy_f107_MC_1000000.mat
    containing:
        xdata1_mc_lognormal, ydata1_mc, f107_selected, Dst_selected,
        Ap_selected, noise_ratio_current

  UT input file:
        sigma_points_f107_ut_3.mat
    containing:
        xdata1_ut, ydata1_ut, f107_selected, Dst_selected, Ap_selected,
        noise_ratio_current, Wm, Wc

These .mat files are intended to be used as downstream inputs to the global
thermosphere-density prediction framework.

===============================================================================
Author information
===============================================================================
Author / script owner: Dr. Pugazhenthi Sivasankar
Affiliation: XBai Research Group, Department of Mechanical and Aerospace
             Engineering, Rutgers University
Research context: Space-weather input uncertainty propagation through a global
                  thermosphere-density prediction framework using MC and UT.
Documentation note: Header documentation prepared with ChatGPT assistance.
===============================================================================
%}

% Begin with a clean workspace, figures and a command window:
clear variables; close all; clc;

% Address of the saved variables: (parent directory string)
parent_dir_str = 'Output/Swarm/';

% Set the Gaussian noise for the MC and UT inputs:
% noise_ratio = [0.1 0.25 0.5 0.75 1];
% noise_ratio = [0.5 0.75 1];
% noise_ratio = [0.75];
noise_ratio = [0.1 0.25];


%% Create a table of cases based on the data given by Ruochen:

% Build the table with the desired content
Report_case     = [1; 2; 3];
Satellite       = ["Swarm-C"; "Swarm-C"; "Swarm-C"];    % string type
sw_indicator    = ["High_f107"; "Low_f107"; "Medium_f107"];

% Suitable table name
tblSpaceWeather = table(Report_case, Satellite, sw_indicator);

%% Create and save the MC samples and UT sigma points:

% Set the number of samples for MC:
n_samples_MC = 1e6;

% Set the number of sigma points for UT:
n_samples_ut = 3;

% String for the noisy data for the current case (MC):
data_str_mc = strcat('/Noisy_f107_MC_', num2str(n_samples_MC), '.mat');

% String for the noisy data for the current case (UT):
data_str_ut = strcat('/sigma_points_f107_ut_', num2str(n_samples_ut), '.mat');

% Index of F10.7 according to the description given by Ruochen:
ind_interest = 12;

% Bin count suitable for histogram (adjust if desired)
nbins = 100;

% Loop over the noise ratios:
for jj = 1:length(noise_ratio)

    % Get the current noise ratio:
    noise_ratio_current = noise_ratio(jj);

    % Loop over each row of the space weather table:
    for ii = 1:length(tblSpaceWeather.Report_case)
    % for ii = 1:1
    
        child_dir_str_input = strcat(char(tblSpaceWeather.sw_indicator(ii)),...
        '/GPinputS.mat');
        child_dir_str_output = strcat(char(tblSpaceWeather.sw_indicator(ii)),...
        '/GPoutputS.mat');
    
        % Load the data:
        load(strcat(parent_dir_str, child_dir_str_input));
        load(strcat(parent_dir_str, child_dir_str_output));
    
        % Extract the very first row:
        x1 = xdata(1, :);
        y1 = ydata(1);
    
        % Extract the SW indices:
        f107_selected = x1(:, 12);
        Dst_selected = x1(:, 15);  
        Ap_selected = x1(:, 16);
        
        % Include the SW indices as part of the table:
        tblSpaceWeather.F107_value(ii) = f107_selected;
        tblSpaceWeather.Dst_value(ii) = Dst_selected;
        tblSpaceWeather.Ap_value(ii) = Ap_selected;
    
        % Include the satellite's latitude, longitude and altitude as part of
        % the table:
        tblSpaceWeather.sat_lat(ii) = x1(:, 7);
        tblSpaceWeather.sat_lon(ii) = x1(:, 8);
        tblSpaceWeather.sat_alt(ii) = x1(:, 9);  
    
        % Include the noise ratio as part of the table:
        tblSpaceWeather.noise_ratio(ii) = noise_ratio_current;
    
        % Repeat the matrix for easy insertion of MC noise:
        xdata1_mc = repmat(x1, n_samples_MC, 1);
        ydata1_mc = repmat(y1, n_samples_MC, 1);
    
        % Draw F107 from lognormal distribution:
        xdata1_mc_lognormal = xdata1_mc;
        mean_linear = tblSpaceWeather.F107_value(ii);
        variance_linear = (noise_ratio_current*tblSpaceWeather.F107_value(ii))^2;
        xdata1_mc_lognormal(:, ind_interest) = generate_log10_lognormal_samples(mean_linear, variance_linear, ...
            n_samples_MC);

        % Repeat the matrix for easy insertion of UT sigma points:
        xdata1_ut = repmat(x1, n_samples_ut, 1);
        ydata1_ut = repmat(y1, n_samples_ut, 1);

        % Creating sigma points for F10.7:
        [ut_xdata1_f107, Wm, Wc, details] = ut_sigma_points_log10_lognormal( ...
        mean_linear, variance_linear, 'Verbose', false);

        % Substituting with sigma points:
        xdata1_ut(:, ind_interest) = ut_xdata1_f107;         

        figure('Color','w','Name','Input Histograms (1x1)');
        tiledlayout(1,1,'TileSpacing','compact','Padding','compact');
       
        
        %%% 1️⃣ Histogram of F10.7 
        nexttile;
        data1 = xdata1_mc_lognormal(:, ind_interest);
        mu1 = mean(data1, 'omitnan');
        
        histogram(data1, nbins, 'Normalization','percentage', ...
            'FaceColor',[0.2 0.45 0.85], 'EdgeColor','none', 'FaceAlpha',0.85);
        grid on;
        hold on;
        xline(mu1, '--r', 'LineWidth',1.5);
        xline(ut_xdata1_f107(1), '--b', 'LineWidth',1);
        xline(ut_xdata1_f107(2), '-b', 'LineWidth',1);
        xline(ut_xdata1_f107(3), '-b', 'LineWidth',1);   
        legend( 'Mean MC data', 'UT sigma point 1', ...
            'UT sigma point 2', 'UT sigma point 3')        
        xlabel('$F_{10.7}$ index','Interpreter','latex');
        ylabel('% of the total samples');
        title(sprintf('Histogram of input lognormal $F_{10.7}$ and its UT sigma points(case %d %s with a noise ratio of %0.2f)', ...
            tblSpaceWeather.Report_case(ii), ...
            tblSpaceWeather.sw_indicator{ii}, noise_ratio_current ),'Interpreter','latex');
        ax = gca;
        ax.FontName = 'Times New Roman';
        ax.FontSize = 15;        
        
       
        % Obtain the full file path for output directory and create one if it
        % does not already exist:
        output_dir = strcat(parent_dir_str, char(tblSpaceWeather.sw_indicator(ii)), ...
            '/Noise_ratio_', num2str(noise_ratio_current),'/lognormal_mc_ut');
        if ~exist(output_dir, 'dir')
            mkdir(output_dir);
        end
    
        % Save the MC samples in the corresponding folder:
        save(strcat(output_dir, data_str_mc), 'xdata1_mc_lognormal', 'ydata1_mc', ...
            'f107_selected', 'Dst_selected', 'Ap_selected', 'noise_ratio_current', "-v7.3");    

        % Save the UT sigma points in the corresponding folder:
        save(strcat(output_dir, data_str_ut), 'xdata1_ut', 'ydata1_ut', ...
            'f107_selected', 'Dst_selected', 'Ap_selected', 'noise_ratio_current',  ...
            'Wm', 'Wc', "-v7.3");     
    
    end

end

disp(tblSpaceWeather)

%% Local functions:

function samples = generate_log10_lognormal_samples(mean_linear, variance_linear, ...
    num_samples, varargin)
% GENERATE_LOG10_LOGNORMAL_SAMPLES  Draw positive samples using a base-10 lognormal model.
%
%   samples = generate_log10_lognormal_samples(mean_linear, variance_linear, num_samples)
%   samples = generate_log10_lognormal_samples(..., 'Verbose', true/false, 'IndexName', 'F10.7')
%
% Inputs (linear domain):
%   mean_linear     - Target arithmetic mean of the index (scalar, > 0)
%   variance_linear - Target variance of the index (scalar, >= 0)
%   num_samples     - Number of samples to generate (positive integer)
%
% Name-Value Pairs (optional):
%   'Verbose'   : logical, default false. If true, prints parameter details and realized stats.
%   'IndexName' : char/string, default 'Index'. Used only in printed messages.
%
% Output:
%   samples - [num_samples x 1] positive samples drawn from X = 10^Y,
%             with Y ~ N(mu10, sigma10^2) chosen to match the given
%             (mean, variance) in the linear domain.
%
% Model:
%   Let Y = log10(X). Assume Y ~ N(mu10, sigma10^2).
%   Given target arithmetic mean m and variance s^2 in the *linear* domain:
%       sigma_e^2 = ln(1 + s^2/m^2),          (natural-log variance)
%       mu_e      = ln(m) - 0.5*sigma_e^2,    (natural-log mean)
%   Convert to base-10 parameters:
%       mu10      = mu_e / ln(10)
%       sigma10   = sqrt(sigma_e^2) / ln(10)
%   Sampling:
%       X = 10^( mu10 + sigma10 * Z ), with Z ~ N(0,1).
%
% Notes:
%   - This is convenient when your workflows use log10 scaling.
%   - Ensure any downstream statistics are consistent with base-10 usage.

    % ---------- Parse optional args ----------
    p = inputParser;
    p.addParameter('Verbose',   false, @(x)islogical(x) && isscalar(x));
    p.addParameter('IndexName', 'Index', @(x)ischar(x) || isstring(x));
    p.parse(varargin{:});
    verbose   = p.Results.Verbose;
    indexName = char(p.Results.IndexName);

    % ---------- Input validation ----------
    if ~isscalar(mean_linear) || ~(mean_linear > 0)
        error('mean_linear must be a positive scalar.');
    end
    if ~isscalar(variance_linear) || ~(variance_linear >= 0)
        error('variance_linear must be a nonnegative scalar.');
    end
    if ~isscalar(num_samples) || num_samples <= 0 || floor(num_samples) ~= num_samples
        error('num_samples must be a positive integer scalar.');
    end

    % ---------- Natural-log parameters from linear (mean, var) ----------
    m  = mean_linear;
    s2 = variance_linear;
    cv2       = s2 / (m^2);           % coefficient of variation squared
    sigma_e2  = log(1 + cv2);         % Var[ln X]
    mu_e      = log(m) - 0.5*sigma_e2;% E[ln X]

    % ---------- Convert to base-10 log parameters ----------
    ln10    = log(10);
    mu10    = mu_e   / ln10;
    sigma10 = sqrt(sigma_e2) / ln10;

    % ---------- Generate samples ----------
    Z = randn(num_samples, 1);            % standard normal
    samples = 10.^(mu10 + sigma10 * Z);   % X = 10^Y

    % ---------- Optional reporting ----------
    if verbose
        fprintf('%s Sample Generation (base-10 lognormal)\n', indexName);
        fprintf('  Target mean (linear):       %.6f\n', m);
        fprintf('  Target variance (linear):   %.6f\n', s2);
        fprintf('  Natural-log  mu_e:          %.6f\n', mu_e);
        fprintf('  Natural-log  sigma_e:       %.6f\n', sqrt(sigma_e2));
        fprintf('  Base-10      mu10:          %.6f\n', mu10);
        fprintf('  Base-10      sigma10:       %.6f\n', sigma10);
        fprintf('  Sample mean (linear):       %.6f\n', mean(samples));
        fprintf('  Sample variance (linear):   %.6f\n', var(samples));
    end
end


function [xSigma, Wm, Wc, details] = ut_sigma_points_log10_lognormal( ...
        mean_linear, variance_linear, varargin)
% UT_SIGMA_POINTS_LOG10_LOGNORMAL
%   Build Unscented Transform sigma points for a base-10 lognormal RV X.
%   The UT is applied in log10 space (Gaussian), then mapped to X-space.
%   The function chooses lambda (i.e., alpha) to MATCH the target mean
%   exactly, and chooses beta to MATCH the target variance exactly.
%
%   [xSigma, Wm, Wc, details] = ut_sigma_points_log10_lognormal(m, s2)
%   [ ... ] = ut_sigma_points_log10_lognormal(m, s2, 'Verbose', true)
%
% Inputs:
%   mean_linear     m  > 0  : target arithmetic mean of X
%   variance_linear s2 >= 0 : target variance of X
%
% Name-Value pairs (optional):
%   'Kappa'   : default 0 (standard choice). For n=1, lambda = alpha^2*(1+kappa)-1.
%   'Verbose' : default false. If true, prints diagnostics.
%
% Outputs:
%   xSigma : [3x1] sigma points in X-space (strictly positive)
%   Wm     : [3x1] UT weights for the mean
%   Wc     : [3x1] UT weights for the covariance
%   details: struct with fields
%              .mu10, .sigma10, .alpha, .beta, .kappa, .lambda
%              .Xmean_check, .Xvar_check   (reconstructed mean/var)
%
% Theory/References:
%   - UT/Scaled-UT: Julier & Uhlmann (2004), Proc. IEEE 92(3); 
%                   Wan & van der Merwe (2000).
%   - Base-10 lognormal parameterization and moments:
%       If log10 X ~ N(mu10, sigma10^2), then ln X ~ N(mu, sigma^2)
%       with mu = mu10*ln(10), sigma = sigma10*ln(10).
%       E[X] = exp(mu + 0.5*sigma^2);
%       Var[X] = (exp(sigma^2)-1) * exp(2*mu + sigma^2).
%
% Notes:
%   - n=1 (one random variable) -> three sigma points.
%   - We solve for lambda to match the MEAN exactly; then choose beta to
%     match the VARIANCE exactly (because W0^c differs linearly via beta).
%   - If s2=0, all sigma points collapse to the mean and Wm/Wc are [1;0;0].

    % ---------- Parse inputs ----------
    p = inputParser;
    p.addParameter('Kappa',   0, @(x) isscalar(x) && isreal(x));
    p.addParameter('Verbose', false, @(x) isscalar(x) && islogical(x));
    p.parse(varargin{:});
    kappa   = p.Results.Kappa;
    verbose = p.Results.Verbose;

    % ---------- Validate ----------
    m  = mean_linear;
    s2 = variance_linear;
    if ~(isscalar(m)  && m  > 0), error('mean_linear must be positive scalar.'); end
    if ~(isscalar(s2) && s2 >= 0), error('variance_linear must be nonnegative scalar.'); end

    % ---------- Natural-log parameters from linear (m, s2) ----------
    % See e.g. Aitchison & Brown (1957); Johnson, Kotz, Balakrishnan (1994).
    cv2      = s2 / (m^2);
    sigma_e2 = log(1 + cv2);           % Var[ln X]
    mu_e     = log(m) - 0.5*sigma_e2;  % E[ln X]

    % ---------- Convert to base-10 log parameters ----------
    ln10   = log(10);
    mu10   = mu_e   / ln10;
    sigma10 = sqrt(sigma_e2) / ln10;

    % Handle the degenerate case separately
    if s2 == 0
        xSigma = [m; m; m];
        Wm     = [1; 0; 0];
        Wc     = [1; 0; 0];  % any set summing to 1 works; this is fine
        details = struct('mu10',mu10,'sigma10',sigma10, ...
                         'alpha',0,'beta',0,'kappa',kappa,'lambda',-1, ...
                         'Xmean_check',m,'Xvar_check',0);
        if verbose
            fprintf('[UT base-10 lognormal] Degenerate s2=0 -> all points at mean.\n');
        end
        return;
    end

    % ---------- n=1, choose lambda (equiv. alpha) to match MEAN exactly ----------
    n = 1;
    % For given lambda (> -n), the sigma points in log10-space are:
    %   z0 = mu10, z± = mu10 ± sqrt(n+lambda)*sigma10
    % and mean-weights:
    %   W0m = lambda/(n+lambda), W±m = 1/(2*(n+lambda))
    % The X-space sigma points are x_i = 10^(z_i).
    % We solve f(lambda) = sum_i Wm(lambda)*x_i(lambda) - m = 0.

    % Ensure n+lambda > 0, so lambda > -1. We'll search on (-0.999, Lmax).
    Llo = -0.999; Lhi = 50; % wide range; adjust if needed
    mean_fun = @(lam) ut_mean_X_given_lambda(lam, mu10, sigma10, n) - m;
    % Try to bracket a root; expand upper bound if necessary.
    fLlo = mean_fun(Llo); fLhi = mean_fun(Lhi);
    tries = 0;
    while sign(fLlo) == sign(fLhi) && tries < 6
        Lhi = Lhi * 2;
        fLhi = mean_fun(Lhi);
        tries = tries + 1;
    end
    if sign(fLlo) == sign(fLhi)
        error('Failed to bracket a root for lambda when matching the mean.');
    end
    lambda = fzero(mean_fun, [Llo, Lhi]);

    % From lambda, with given kappa, get alpha (n=1):
    % lambda = alpha^2*(1+kappa) - n  =>  alpha = sqrt((lambda + n)/(1+kappa))
    alpha = sqrt((lambda + n)/(1 + kappa));

    % ---------- Build sigma points in log10-space ----------
    c   = sqrt(n + lambda) * sigma10;   % spread in log10-space
    z0  = mu10;
    zP  = mu10 + c;
    zM  = mu10 - c;

    % Map to X-space (strictly positive)
    x0 = 10.^z0;
    xP = 10.^zP;
    xM = 10.^zM;
    xSigma = [x0; xP; xM];

    % ---------- Mean weights ----------
    W0m = lambda/(n + lambda);
    Wpm = 1/(2*(n + lambda));
    Wm  = [W0m; Wpm; Wpm];

    % Reconstructed mean (should match m within solver tolerance)
    Xmean = Wm.' * xSigma;

    % ---------- Covariance weights: choose beta to match VARIANCE exactly ----------
    % Wc = Wm except W0c = W0m + (1 - alpha^2 + beta).
    % For n=1, alpha^2 = (lambda + 1)/(1 + kappa).
    % Let V0 = sum(Wm .* (x_i - m)^2); then
    %   s2_target = V0 + (1 - alpha^2 + beta) * (x0 - m)^2
    % Hence solve for beta (linear, closed form):
    V0 = sum(Wm .* (xSigma - m).^2);
    denom = (x0 - m)^2;
    if denom < eps
        % If x0 happens to equal the mean, W0c adjustment has no effect;
        % in this rare case the UT with mean-weights already matches var.
        beta = 0;
    else
        beta = (s2 - V0)/denom + alpha^2 - 1;
    end
    W0c = W0m + (1 - alpha^2 + beta);
    Wc  = [W0c; Wpm; Wpm];

    % Reconstructed variance check
    Xvar = sum(Wc .* (xSigma - m).^2);

    % ---------- Reporting ----------
    details = struct('mu10',mu10,'sigma10',sigma10, ...
                     'alpha',alpha,'beta',beta,'kappa',kappa,'lambda',lambda, ...
                     'Xmean_check',Xmean,'Xvar_check',Xvar);

    if verbose
        fprintf('[UT base-10 lognormal]\n');
        fprintf('  Target  mean: %.9g | Reconstructed mean: %.9g\n', m, Xmean);
        fprintf('  Target  var : %.9g | Reconstructed var : %.9g\n', s2, Xvar);
        fprintf('  mu10=%.6g, sigma10=%.6g, alpha=%.6g, beta=%.6g, kappa=%.6g, lambda=%.6g\n', ...
                mu10, sigma10, alpha, beta, kappa, lambda);
    end
end

% ---------- helper ----------
function mX = ut_mean_X_given_lambda(lambda, mu10, sigma10, n)
    % Build mean weights and points in log10 space for given lambda (n=1).
    c    = sqrt(n + lambda) * sigma10;
    z0   = mu10;
    zP   = mu10 + c;
    zM   = mu10 - c;
    x0   = 10.^z0;
    xP   = 10.^zP;
    xM   = 10.^zM;
    W0m  = lambda/(n + lambda);
    Wpm  = 1/(2*(n + lambda));
    mX   = W0m*x0 + Wpm*(xP + xM);
end

