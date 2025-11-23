
%{

This script checks the first three statistical moments propagation through
linear and nonlinear systems. The theoretical moments are compared with
the moments obtained from the UT recombination equations.


Authors:
1) Dr. Pugazhenthi Sivasankar, Post-doctoral researcher, XBai Research
Group, Rutgers University.

%}


% Begin with a clean workspace, figures and a command window:
clear variables; close all; clc;

% Set the system coefficients:
sys_coeff_a = 5;
sys_coeff_b = 3;

% Build the table with the desired content
System_type     = ["linear"; "Nonlinear"];  % string type
% target_mean       = [10; 10];    
% target_variance    = [8; 8]

target_mean       = [100; 100];
% target_variance    = [2; 2];
target_variance    = [4; 4];
% target_variance    = [0.5; 0.5];

% Obtain the mean and the variance in the base-10 exponential domain:
sigma_e2 = log(1 + target_variance(1)/target_mean(1)^2);
mu_e       = log(target_mean(1)) - 0.5*sigma_e2;
mean_10    = mu_e/log(10);
variance_10    = (sqrt(sigma_e2)/log(10))^2;
mean_10 = [mean_10; mean_10];
variance_10 = [variance_10; variance_10];

% Create a table:
input_table = table(System_type, target_mean, target_variance);

% Compute the theoretical mean, variance, and skewness of the input:
mean_x_ip_th = (10^(mean_10(1)))*exp(0.5*(log(10)^2)*variance_10(1));
var_x_ip_th = (exp(log(10)^2*variance_10(1)) - 1)*(10^(2*mean_10(1)))*exp(log(10)^2*variance_10(1));
skew_x_ip_th = (exp(log(10)^2*variance_10(1)) + 2)*sqrt(exp(log(10)^2*variance_10(1)) - 1);

% Append the theoretical mean, variance, and skewness of the input to the
% table:
input_table.mean_x_ip_th = [mean_x_ip_th; mean_x_ip_th];
input_table.var_x_ip_th = [var_x_ip_th; var_x_ip_th];
input_table.skew_x_ip_th = [skew_x_ip_th; skew_x_ip_th];

% Obtain the sigma points using log-domain UT:
[in_data_ut, Wm, Wc, details] = ut_sigma_points_log10_lognormal( ...
target_mean(1), target_variance(1), 'Verbose', false);

% Input statistical moments from sigma points:
mean_x_ip_ut = Wm.' * in_data_ut;
var_x_ip_ut  = sum( Wc .* (in_data_ut - mean_x_ip_ut).^2 );  % variance uses Wc
ut_mu3  = sum( Wm .* (in_data_ut - mean_x_ip_ut).^3 );  % third central (use Wm)
skew_x_ip_ut = ut_mu3 / (var_x_ip_ut^(3/2)); 

% Append the input statistical moments from sigma points to the table:
input_table.mean_x_ip_ut = [mean_x_ip_ut; mean_x_ip_ut];
input_table.var_x_ip_ut = [var_x_ip_ut; var_x_ip_ut];
input_table.skew_x_ip_ut = [skew_x_ip_ut; skew_x_ip_ut];
disp('Input random variable statistics:')
disp(input_table)

% Create the output table:
output_table = table(System_type);

% Compute the theoretical mean, variance and skewness for the linear and
% nonlinear systems:
mean_x_op_th(1) = sys_coeff_a*input_table.mean_x_ip_th(1) + sys_coeff_b;
mean_x_op_th(2) = sys_coeff_a*(10^(2*mean_10(1)))*exp(2*(log(10)^2)*variance_10(1)) + sys_coeff_b;

var_x_op_th(1) = sys_coeff_a^2*var_x_ip_th(1);
var_x_op_th(2) = sys_coeff_a^2*(exp(4*log(10)^2*variance_10(1)) - 1)*(10^(2*2*mean_10(1)))*exp(4*log(10)^2*variance_10(1));

skew_x_op_th(1) = sign(sys_coeff_a)*skew_x_ip_th(1);
skew_x_op_th(2) = sign(sys_coeff_a)*(exp(4*log(10)^2*variance_10(1)) + 2)*sqrt(exp(4*log(10)^2*variance_10(1)) - 1);

% Append to the table:
output_table.mean_x_op_th = mean_x_op_th';
output_table.var_x_op_th = var_x_op_th';
output_table.skew_x_op_th = skew_x_op_th';

% Compute the linear system output:
lin_sys_op_data_ut = lin_system(in_data_ut, sys_coeff_a, sys_coeff_b);

% Output statistical moments from sigma points for the linear system:
mean_x_op_ut(1) = Wm.' * lin_sys_op_data_ut;
var_x_op_ut(1)  = sum( Wc .* (lin_sys_op_data_ut - mean_x_op_ut(1)).^2 );  % variance uses Wc
ut_mu3  = sum( Wm .* (lin_sys_op_data_ut - mean_x_op_ut(1)).^3 );  % third central (use Wm)
skew_x_op_ut(1) = ut_mu3 / (var_x_op_ut(1)^(3/2));

% Compute the nonlinear system output:
nonlin_sys_op_data_ut = non_lin_system(in_data_ut, sys_coeff_a, sys_coeff_b);

% Output statistical moments from sigma points for the nonlinear system:
mean_x_op_ut(2) = Wm.' * nonlin_sys_op_data_ut;
var_x_op_ut(2)  = sum( Wc .* (nonlin_sys_op_data_ut - mean_x_op_ut(2)).^2 );  % variance uses Wc
ut_mu3  = sum( Wm .* (nonlin_sys_op_data_ut - mean_x_op_ut(2)).^3 );  % third central (use Wm)
skew_x_op_ut(2) = ut_mu3 / (var_x_op_ut(2)^(3/2));

% Append to the table:
output_table.mean_x_op_ut = mean_x_op_ut';
output_table.var_x_op_ut = var_x_op_ut';
output_table.skew_x_op_ut = skew_x_op_ut';

% Display the final results:
disp('Output random variable statistics:')
disp(output_table)

%% Local functions:

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