
%{

This script compares the statistics of the Ap from the training
data of the NN with that of the statistics of the lognormal input SW
indices distribution created for the MC simulation.
 
version 1:
 - Three different cases of Ap are considered with the SWARM-C satellite.
 - Monte-Carlo samples of Ap are created for each case and appended to
 the rest of the deterministic outputs. The number of Monte-Carlo samples
 is user-defined.


Authors:
1) Ruochen Wang, Ph.D. student, XBai Research Group, Rutgers University.
2) Dr. Pugazhenthi Sivasankar, Post-doctoral researcher, XBai Research
Group, Rutgers University.

%}

% Begin with a clean workspace, figures and a command window:
clear variables; close all; clc;

% load the training dataset for inverse normalization
load('GPinput.mat')

% Extract the F10.7 values from the training data:
train_f107 = xdata(:, 12);

% Extract the Ap values from the training data:
train_ap = xdata(:, 16);

% Clear the 'xdata' to save some RAM:
clear('xdata');

% Compute the first three statistical moments of F10.7:
f107_stats.mean_f107_train = mean(train_f107);
f107_stats.std_f107_train = std(train_f107);
f107_stats.skew_f107_train = skewness(train_f107);

% Compute the extreme values of F10.7:
f107_stats.min_f107_train = min(train_f107);
f107_stats.max_f107_train = max(train_f107);
f107_stats.range_f107_train = f107_stats.max_f107_train - f107_stats.min_f107_train;

% Compute the first three statistical moments of Ap:
ap_stats.mean_ap_train = mean(train_ap);
ap_stats.std_ap_train = std(train_ap);
ap_stats.skew_ap_train = skewness(train_ap);

% Compute the extreme values of Ap:
ap_stats.min_ap_train = min(train_ap);
ap_stats.max_ap_train = max(train_ap);
ap_stats.range_ap_train = ap_stats.max_ap_train - ap_stats.min_ap_train;

% Build the table with the desired content
SW_index     = ["F107"; "Ap";];
Mean         = [f107_stats.mean_f107_train; ap_stats.mean_ap_train];
Std          = [f107_stats.mean_f107_train; ap_stats.std_ap_train];
Skewness     = [f107_stats.mean_f107_train; ap_stats.skew_ap_train];
Minimum      = [f107_stats.mean_f107_train; ap_stats.min_ap_train];
Maximum      = [f107_stats.mean_f107_train; ap_stats.max_ap_train];

% Suitable table name
train_SW_stats = table(SW_index, Mean, Std, Skewness, Minimum, Maximum);

% Display the table:
disp('Statistics of F10.7 and Ap in the NN training data:')
disp(train_SW_stats)

% Address of the saved variables: (parent directory string)
parent_dir_str = 'C:\Users\Pugazhenthi\Box\Xiaoli_Bai_Rutgers\Thermosphere_density_predictions\Pugal_new_v1\Output\Swarm\';

noise_ratio = [0.25];
% noise_ratio = [0.1];

% Set the linear system coefficients:
lin_sys_coeff_a = 5;
lin_sys_coeff_b = 3;

% Set the nonlinear system coefficients:
nonlin_sys_coeff_a = 5;
nonlin_sys_coeff_b = 3;

%% Create a table of cases based on the data given by Ruochen:

% Build the table with the desired content
Report_case     = [1; 2; 3];
Satellite       = ["Swarm-C"; "Swarm-C"; "Swarm-C"];    % string type
sw_indicator    = ["High_Ap"; "Low_Ap"; "Medium_Ap"];

% Suitable table name
tblSpaceWeather = table(Report_case, Satellite, sw_indicator);
disp(tblSpaceWeather)

%% Create and save the MC samples:

% Set the number of samples for MC:
n_samples_MC = 1e6;

% Bin count suitable for histogram (adjust if desired)
nbins = 100;

% Create a figure for the comparitive histograms:
fig1 = figure('Color','w','Name','Comparitive Histograms (1x3)');
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
sgtitle(sprintf(['Comparitive histograms of the NN training $A_p$ data and' ...
    ' uncertain input $A_p$ (lognormal) with a noise ratio of %0.2f)'], ...
     noise_ratio), 'Interpreter','latex', 'FontSize',15);

fig2 = figure('Color','w','Name','MC - UT comparison (1x3)');
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
sgtitle(sprintf(['Comparision of the uncertain lognormal input $A_p$ histogram (MC)' ...
    ' and the sigma points (UT) with a noise ratio of %0.2f)'], ...
     noise_ratio), 'Interpreter','latex', 'FontSize',15);

fig3 = figure('Color','w','Name','MC - UT comparison (1x3)');
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
% sgtitle(sprintf(['Comparision of the linear system output using uncertain lognormal input $A_p$ histogram (MC)' ...
%     ' and the sigma points (UT) with a noise ratio of %0.2f)'], ...
%      noise_ratio), 'Interpreter','latex', 'FontSize',15);
sgtitle(sprintf(['Comparision of the nonlinear system output using uncertain lognormal input $A_p$ histogram (MC)' ...
    ' and the sigma points (UT) with a noise ratio of %0.2f)'], ...
     noise_ratio), 'Interpreter','latex', 'FontSize',15);

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
        mean_linear = tblSpaceWeather.Ap_value(ii);
        variance_linear = (noise_ratio_current*tblSpaceWeather.Ap_value(ii))^2;
        xdata1_mc_lognormal(:, 16) = generate_log10_lognormal_samples(mean_linear, variance_linear, ...
            n_samples_MC);  

        % Extract the F10.7 values from the training data:
        input_ap = xdata1_mc_lognormal(:, 16);     

        % Compute the first three statistical moments of F10.7:
        input_ap_stats{ii}.mean_ap_train = mean(input_ap);
        input_ap_stats{ii}.std_ap_train = std(input_ap);
        input_ap_stats{ii}.skew_ap_train = skewness(input_ap);
        
        % Compute the extreme values of F10.7:
        input_ap_stats{ii}.min_ap_train = min(input_ap);
        input_ap_stats{ii}.max_ap_train = max(input_ap);
        input_ap_stats{ii}.range_ap_train = input_ap_stats{ii}.max_ap_train - ...
        input_ap_stats{ii}.min_ap_train;

        %%% 1️⃣ Histogram of histogram of training and case(ii) Ap
        figure(fig1);
        nexttile;
        histogram(train_ap, nbins, 'Normalization','percentage', ...
            'FaceColor','blue', 'EdgeColor','none', 'FaceAlpha',0.5);
        grid on;
        hold on;
        xline(ap_stats.mean_ap_train, '--b', 'LineWidth',1);
        xline(ap_stats.min_ap_train, '-b', 'LineWidth',1);
        xline(ap_stats.max_ap_train, '-b', 'LineWidth',1);
        histogram(input_ap, nbins, 'Normalization','percentage', ...
            'FaceColor','red', 'EdgeColor','none', 'FaceAlpha',0.5);
        xline(input_ap_stats{ii}.mean_ap_train, '--r', 'LineWidth',1);
        xline(input_ap_stats{ii}.min_ap_train, '-r', 'LineWidth',1);
        xline(input_ap_stats{ii}.max_ap_train, '-r', 'LineWidth',1);        
        legend('Training data', 'Mean train data', 'Min. train data', ...
            'Max. train data', 'uncertain i/p data', 'Mean i/p data', ...
            'Min. i/p data', 'Max. i/p data')
        xlabel('$A_p$ index','Interpreter','latex');
        ylabel('% of the total samples');
        title(sprintf('Training data and uncertain input (case %d %s)', ...
            tblSpaceWeather.Report_case(ii), tblSpaceWeather.sw_indicator{ii}),'Interpreter','latex');
        ax = gca;
        ax.FontName = 'Times New Roman';
        ax.FontSize = 15;  

        ovl_coeff(ii) = kde_overlap(train_ap, input_ap)

        [xSigma, Wm, Wc, details] = ut_sigma_points_log10_lognormal( ...
        mean_linear, variance_linear, 'Verbose', false);

        % Moments from sigma points:
        ut_mean = Wm.' * xSigma;
        ut_ap_stats{ii}.mean_ap_train = ut_mean;                       % mean uses Wm
        ut_var  = sum( Wc .* (xSigma - ut_mean).^2 );  % variance uses Wc
        ut_ap_stats{ii}.std_ap_train = sqrt(ut_var);
        ut_mu3  = sum( Wm .* (xSigma - ut_mean).^3 );  % third central (use Wm)
        ut_ap_stats{ii}.skew_ap_train = ut_mu3 / (ut_var^(3/2));  

        figure(fig2);
        nexttile;
        histogram(input_ap, nbins, 'Normalization','percentage', ...
            'FaceColor','red', 'EdgeColor','none', 'FaceAlpha',0.5);
        grid on;
        hold on;
        xline(input_ap_stats{ii}.mean_ap_train, '--r', 'LineWidth',1);
        xline(input_ap_stats{ii}.min_ap_train, '-r', 'LineWidth',1);
        xline(input_ap_stats{ii}.max_ap_train, '-r', 'LineWidth',1); 
        xline(xSigma(1), '--b', 'LineWidth',1);
        xline(xSigma(2), '-b', 'LineWidth',1);
        xline(xSigma(3), '-b', 'LineWidth',1);        
        legend('uncertain i/p data', 'Mean i/p data', ...
            'Min. i/p data', 'Max. i/p data', 'UT sigma point 1', ...
            'UT sigma point 2', 'UT sigma point 3')
        xlabel('$A_p$ index','Interpreter','latex');
        ylabel('% of the total samples');
        title(sprintf('MC (vs.) UT (case %d %s)', ...
            tblSpaceWeather.Report_case(ii), tblSpaceWeather.sw_indicator{ii}),'Interpreter','latex');
        ax = gca;
        ax.FontName = 'Times New Roman';
        ax.FontSize = 15;   

        % % Pass the MC samples through the linear system:
        % output_mc = lin_system(input_ap, lin_sys_coeff_a, lin_sys_coeff_b);
        % 
        % % Pass the UT sigma points through the linear system:
        % output_ut = lin_system(xSigma, lin_sys_coeff_a, lin_sys_coeff_b);

        % Pass the MC samples through the nonlinear system:
        output_mc = non_lin_system(input_ap, nonlin_sys_coeff_a, nonlin_sys_coeff_b);

        % Pass the UT sigma points through the nonlinear system:
        output_ut = non_lin_system(xSigma, nonlin_sys_coeff_a, nonlin_sys_coeff_b);        

        % Compute the first three statistical moments of the output MC:
        output_mc_ap_stats{ii}.mean_ap = mean(output_mc);
        output_mc_ap_stats{ii}.std_ap = std(output_mc);
        output_mc_ap_stats{ii}.skew_ap = skewness(output_mc);
        output_mc_ap_stats{ii}.min_ap = min(output_mc);
        output_mc_ap_stats{ii}.max_ap = max(output_mc);        

        % Output statistical moments from sigma points:
        ut_mean_output = Wm.' * output_ut;
        output_ut_ap_stats{ii}.mean_ap = ut_mean_output;   % mean uses Wm
        ut_var_output  = sum( Wc .* (output_ut - ut_mean_output).^2 );  % variance uses Wc
        output_ut_ap_stats{ii}.std_ap = sqrt(ut_var_output);
        ut_mu3_output  = sum( Wm .* (output_ut - ut_mean_output).^3 );  % third central (use Wm)
        output_ut_ap_stats{ii}.skew_ap = ut_mu3_output / (ut_var_output^(3/2));         

      
        figure(fig3);
        nexttile;
        histogram(output_mc, nbins, 'Normalization','percentage', ...
            'FaceColor','red', 'EdgeColor','none', 'FaceAlpha',0.5);
        grid on;
        hold on;
        xline(output_mc_ap_stats{ii}.mean_ap, '--r', 'LineWidth',1);
        xline(output_mc_ap_stats{ii}.min_ap, '-r', 'LineWidth',1);
        xline(output_mc_ap_stats{ii}.max_ap, '-r', 'LineWidth',1); 
        xline(output_ut(1), '--b', 'LineWidth',1);
        xline(output_ut(2), '-b', 'LineWidth',1);
        xline(output_ut(3), '-b', 'LineWidth',1);        
        legend('MC output', 'Mean MC output', ...
            'Min. MC output', 'Max. MC output', 'UT sigma point 1', ...
            'UT sigma point 2', 'UT sigma point 3')
        % xlabel('Linear system output','Interpreter','latex');
        xlabel('Nonlinear system output','Interpreter','latex');
        ylabel('% of the total samples');
        title(sprintf('MC (vs.) UT (case %d %s)', ...
            tblSpaceWeather.Report_case(ii), tblSpaceWeather.sw_indicator{ii}),'Interpreter','latex');
        ax = gca;
        ax.FontName = 'Times New Roman';
        ax.FontSize = 15;

    end
end



% Build the table with the desired content
SW_index     = ["Ap"; "Ap"; "Ap"];
Mean         = [input_ap_stats{1}.mean_ap_train; input_ap_stats{2}.mean_ap_train;...
    input_ap_stats{3}.mean_ap_train;];
Std          = [input_ap_stats{1}.std_ap_train; input_ap_stats{2}.std_ap_train;...
    input_ap_stats{3}.std_ap_train;];
Skewness     = [input_ap_stats{1}.skew_ap_train; input_ap_stats{2}.skew_ap_train;...
    input_ap_stats{3}.skew_ap_train;];
Minimum      = [input_ap_stats{1}.min_ap_train; input_ap_stats{2}.min_ap_train;...
    input_ap_stats{3}.min_ap_train;];
Maximum      = [input_ap_stats{1}.max_ap_train; input_ap_stats{2}.max_ap_train;...
    input_ap_stats{3}.max_ap_train;];


Mean_ut         = [ut_ap_stats{1}.mean_ap_train; ut_ap_stats{2}.mean_ap_train;...
    ut_ap_stats{3}.mean_ap_train;];
Std_ut          = [ut_ap_stats{1}.std_ap_train; ut_ap_stats{2}.std_ap_train;...
    ut_ap_stats{3}.std_ap_train;];
Skewness_ut     = [ut_ap_stats{1}.skew_ap_train; ut_ap_stats{2}.skew_ap_train;...
    ut_ap_stats{3}.skew_ap_train;];

% Suitable table name
input_SW_stats = table(Report_case, SW_index, Mean, Std, Skewness, Minimum, Maximum);
ut_SW_stats = table(Report_case, SW_index, Mean_ut, Std_ut, Skewness_ut);

% Display the table:
fprintf('Statistics of F10.7 and Ap in the input data (MC) with noise ratio %0.2f\n:', noise_ratio_current)
disp(input_SW_stats) 

fprintf('Statistics of F10.7 and Ap in the input data (UT) with noise ratio %0.2f\n:', noise_ratio_current)
disp(ut_SW_stats) 

% Build the output table for the linear/nonlinear system:
% For the MC samples:
Mean_mc         = [output_mc_ap_stats{1}.mean_ap; output_mc_ap_stats{2}.mean_ap;...
    output_mc_ap_stats{3}.mean_ap];
Std_mc          = [output_mc_ap_stats{1}.std_ap; output_mc_ap_stats{2}.std_ap;...
    output_mc_ap_stats{3}.std_ap];
Skewness_mc     = [output_mc_ap_stats{1}.skew_ap; output_mc_ap_stats{2}.skew_ap;...
    output_mc_ap_stats{3}.skew_ap];
output_mc_stats = table(Report_case, SW_index, Mean_mc, Std_mc, Skewness_mc);

% For the UT sigma points:
Mean_ut         = [output_ut_ap_stats{1}.mean_ap; output_ut_ap_stats{2}.mean_ap;...
    output_ut_ap_stats{3}.mean_ap];
Std_ut          = [output_ut_ap_stats{1}.std_ap; output_ut_ap_stats{2}.std_ap;...
    output_ut_ap_stats{3}.std_ap];
Skewness_ut     = [output_ut_ap_stats{1}.skew_ap; output_ut_ap_stats{2}.skew_ap;...
    output_ut_ap_stats{3}.skew_ap];
output_ut_stats = table(Report_case, SW_index, Mean_ut, Std_ut, Skewness_ut);

% Display the table:
% fprintf('Statistics of the linear system output for the MC samples with noise ratio %0.2f:\n', noise_ratio_current)
fprintf('Statistics of the nonlinear system output for the MC samples with noise ratio %0.2f:\n', noise_ratio_current)
disp(output_mc_stats) 

% fprintf('Statistics of the linear system output for the UT sigma points with noise ratio %0.2f:\n', noise_ratio_current)
fprintf('Statistics of the nonlinear system output for the UT sigma points with noise ratio %0.2f:\n', noise_ratio_current)
disp(output_ut_stats)
        
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

function ovl = kde_overlap(data1, data2)
    % grid spanning both sets
    lo = min([data1(:); data2(:)]);
    hi = max([data1(:); data2(:)]);
    x  = linspace(lo, hi, 2000);
    f1 = ksdensity(data1, x, 'Function','pdf');
    f2 = ksdensity(data2,  x, 'Function','pdf');
    ovl = trapz(x, min(f1,f2));   % in [0,1]
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

function output = lin_system(input, a, b)

%{

This is a simple linear system that follows the equation: y = a*x + b.

Inputs:
input: a (1 x n) or (n x 1) vector

Outputs:
output: a (1 x n) or (n x 1) vector

%}

    output = a.*input + b;
end

function output = non_lin_system(input, a, b)

%{

This is a simple nonlinear system that follows the equation: y = a*x^2 + b.

Inputs:
input: a (1 x n) or (n x 1) vector

Outputs:
output: a (1 x n) or (n x 1) vector

%}

    output = a.*input.^2 + b;
end

