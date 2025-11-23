
%{

This script processes the output from Ruochen's NN for the MC and UT
simulations. Out of the 65 different inputs, only F10.7 is considered as a
random variable with lognormal distribution. Both MC & UT simulation use
the lognormal distributions

version 1:
 - Three different cases of F10.7 are considered with a single satellite:
 Swarm C.
 - The number of Monte-Carlo samples is set to 1 million and the number of
 draws for the consolidated MC is varied for each mean and standard
 deviation combination.
 - Baseline (i.e., without input space weather uncertainty) results
 included in the density post-processing.
 - Data from the respective lognormal folders is post-processed.
 - Modified post-processing with lognormal parametric formulations is used
 for density.
 - input F10.7 data is shown along with the UT sigma points.


Authors:
1) Ruochen Wang, Ph.D. student, XBai Research Group, Rutgers University.
2) Dr. Pugazhenthi Sivasankar, Post-doctoral researcher, XBai Research
Group, Rutgers University.

%}

% Begin with a clean workspace, figures and a command window:
clear variables; close all; clc;

% load the training dataset for inverse normalization
load('Output/Model/Model1/GPoutput.mat')
train_y = ydata';
Ystd = std(train_y);
Ymean = mean(train_y);

% Parent directory string: (constant)
parent_dir_str = 'Output/Swarm/';

% Set the Gaussian noise for the MC and UT inputs:
% noise_ratio = [0.1 0.25 0.5 0.75 1];
% noise_ratio = [0.1 0.5 1];
% noise_ratio = [0.1 0.25];
% noise_ratio = [0.1];
noise_ratio = [0.25];  % to circumvent the memory problem.

% Index of F10.7 in the input data to the NN:
sw_index = 12;

% Set the font size for all figures:
fig_font_size = 20;

%% Create a table of cases based on the data given by Ruochen:

% Build the table with the desired content
Report_case     = [1; 2; 3];
Satellite       = ["Swarm-C"; "Swarm-C"; "Swarm-C"];    % string type
F107_indication = ["High_f107"; "Low_f107"; "Medium_f107"];

% Suitable table name
tblSpaceWeather = table(Report_case, Satellite, F107_indication);

% Display the table:
disp(tblSpaceWeather)

%% load data

% Set the number of samples for MC:
n_samples_MC = 1e6;

% Set the number of sigma points for UT:
n_samples_ut = 3;

% String for the noisy data for the current case (MC):
data_str_mc = strcat('/Noisy_f107_MC_', num2str(n_samples_MC), '.mat');

% String for the noisy data for the current case (UT):
data_str_ut = strcat('/sigma_points_f107_ut_', num2str(n_samples_ut), '.mat');

% Loop over the noise ratios:
for jj = 1:length(noise_ratio)

    % Get the current noise ratio:
    noise_ratio_current = noise_ratio(jj);

    % Switch file pathways based on the different cases:
    for ii = 1:length(tblSpaceWeather.Report_case)
    % for ii = 1:1
    
        % A portion of the file pathway string:
        child_dir_str = strcat(char(tblSpaceWeather.F107_indication(ii)), ...
            '/Noise_ratio_', num2str(noise_ratio_current),'/lognormal_mc_ut');
    
        %%%% Post-processing of MC results %%%%

        % load the true y dataset for MC:
        mc_data = load(strcat(parent_dir_str, child_dir_str, data_str_mc));
    
        % File pathway for MC simulation output from Ruochen's NN:
        save_folder_mc = strcat(parent_dir_str, child_dir_str, '/mc_results.mat');
    
        % Load the output of the MC simulation from the neural network:
        load(save_folder_mc)

        % Convert the struct fields into doubles format:
        mc_results = structfun(@double, mc_results, 'UniformOutput', false);
    
        % Obtain the true accelerometer density:
        mc_results.true_density_after_log = mc_data.ydata1_mc';  
    
        % Obtain the predicted density after logarithmic transformation:
        mc_results.pred_density_after_log = mc_results.mu_1*Ystd + Ymean;
    
        % Obtain the total standard deviation of the predicted density after logarithmic transformation:
        % Combine the aleatoric and epistemic uncertainty appropriately:
        mc_results.sigma_total = sqrt(mc_results.sigma_1.^2 + ...
            mc_results.var_1.^2);
        mc_results.std_pred_density_after_log = mc_results.sigma_total*Ystd;
    
        % Using MATLAB's inbuilt function for base-10 lognormal distribution:
        density_mean_physical_mc = zeros(1,length(mc_results.pred_density_after_log));
        density_std_physical_mc = zeros(1,length(mc_results.pred_density_after_log));
        parfor iii = 1:length(mc_results.pred_density_after_log)
            pd = makedist('Lognormal','mu',mc_results.pred_density_after_log(iii).*log(10), ...
                'sigma',mc_results.std_pred_density_after_log(iii).*log(10));
            density_mean_physical_mc(iii) = mean(pd);
            density_std_physical_mc(iii) = std(pd); 
        end
       mc_results_using_lognormal.density_mean_physical = density_mean_physical_mc;
       mc_results_using_lognormal.density_std_physical = density_std_physical_mc;

        % Histograms for mc_results fields in a single figure (2x2 subplots)
        % --- Safety checks and reshape to column vectors
        vars = {'mu_1','sigma_total','pred_density_after_log',...
            'std_pred_density_after_log'};
        for k = 1:numel(vars)
            assert(isfield(mc_results, vars{k}), ...
                'mc_results.%s is missing.', vars{k});
            v = mc_results.(vars{k});
            assert(isvector(v), 'mc_results.%s must be a vector.', vars{k});
            mc_results.(vars{k}) = v(:); % force column
        end
    
        % Bin count suitable for 1e6 samples (adjust if desired)
        nbins = 100;
    
        % Drawing 100 samples from each lognormal distribution for consolidated statistics:
        n_samples_consolidated = 100;  % default
        % n_samples_consolidated = 1000;
        % n_samples_consolidated = 2000;
        mc_results_using_lognormal.all_samples = sample_positive_lognormal(mc_results_using_lognormal.density_mean_physical, ...
        mc_results_using_lognormal.density_std_physical, n_samples_consolidated);
        mc_results_using_lognormal.density_mean_physical_final = mean(mc_results_using_lognormal.all_samples);
        mc_results_using_lognormal.density_std_physical_final = std(mc_results_using_lognormal.all_samples);
        mc_results_using_lognormal.density_var_physical_final = mc_results_using_lognormal.density_std_physical_final.^2;

        %%%%% Post-processing of UT results %%%%

        % Load the true dataset for the UT:
        ut_data = load(strcat(parent_dir_str, child_dir_str, data_str_ut));    
    
        % File pathway for UT simulation output from Ruochen's NN:  
        save_folder_ut = strcat(parent_dir_str, child_dir_str, '/ut_results.mat');

        % Load the output of the MC simulation from the neural network:
        load(save_folder_ut)

        % Convert the struct fields into doubles format:
        ut_results = structfun(@double, ut_results, 'UniformOutput', false);
    
        % Obtain the true accelerometer density in logarithmic space:
        ut_results.true_density_after_log = ut_data.ydata1_ut';  
   
        % Obtaining the mean of the predicted density in logarithmic space:
        ut_results.pred_density_after_log = ut_results.mu_1*Ystd + Ymean;
    
        % Obtaining the aleatoric standard deviation in logarithmic space:
        ut_results.sigma_std = sqrt(ut_results.beta_1./(ut_results.alpha_1 - 1));
    
        % Obtaining the epistemic standard deviation in logarithmic space:
        ut_results.variance_std = sqrt(ut_results.beta_1./(ut_results.v_1.*(ut_results.alpha_1 - 1)));
    
        % Obtain the total uncertainty in logarithmic space:
        ut_results.sigma_total = sqrt(ut_results.sigma_std.^2 + ut_results.variance_std.^2);
        ut_results.std_pred_density_after_log = ut_results.sigma_total*Ystd;
    
        % Using MATLAB's inbuilt function for base-10 lognormal distribution:
        density_mean_physical_ut = zeros(1,length(ut_results.pred_density_after_log));
        density_std_physical_ut = zeros(1,length(ut_results.pred_density_after_log));
        parfor iii = 1:length(ut_results.pred_density_after_log)
            pd = makedist('Lognormal','mu',ut_results.pred_density_after_log(iii).*log(10), ...
                'sigma',ut_results.std_pred_density_after_log(iii).*log(10));
            density_mean_physical_ut(iii) = mean(pd);
            density_std_physical_ut(iii) = std(pd); 
        end
       ut_results_using_lognormal.density_mean_physical = density_mean_physical_ut;
       ut_results_using_lognormal.density_std_physical = density_std_physical_ut;
       ut_results_using_lognormal.density_var_physical = ut_results_using_lognormal.density_std_physical.^2;

        % Collect and recombine the UT output:
        ut_results_using_lognormal.collected = [ut_results_using_lognormal.density_mean_physical;...
            ut_results_using_lognormal.density_var_physical];
        [ut_results_using_lognormal.mu_recomb, ut_results_using_lognormal.P_recomb] = ut_recombine_outputs(ut_results_using_lognormal.collected, ...
            ut_data.Wm, ut_data.Wc);
        ut_results_using_lognormal.density_mean_physical_final = ut_results_using_lognormal.mu_recomb(1);
        
        % Combine the variances (regular variance + variance of the
        % conditional mean):
        ut_results_using_lognormal.density_var_physical_final = ut_results_using_lognormal.mu_recomb(2) + ...
            ut_results_using_lognormal.P_recomb(1,1);
        ut_results_using_lognormal.density_std_physical_final = sqrt(ut_results_using_lognormal.density_var_physical_final);


        %%%%% Post-processing of baseline results %%%%
        % Collecting the baseline results:
        baseline_results.mu = ut_results.mu_1(1);
        baseline_results.v = ut_results.v_1(1);
        baseline_results.alpha = ut_results.alpha_1(1);
        baseline_results.beta = ut_results.beta_1(1);

        % Obtaining the mean of the predicted density (baseline):
        baseline_results.pred_density_after_log = baseline_results.mu*Ystd + Ymean;        

        % Obtaining the aleatoric standard deviation (baseline):
        baseline_results.sigma_std = sqrt(baseline_results.beta/(baseline_results.alpha - 1));
    
        % Obtaining the epistemic standard deviation (baseline):
        baseline_results.variance_std = sqrt(baseline_results.beta/(baseline_results.v*(baseline_results.alpha - 1)));
    
        % Obtain the total uncertainty (baseline):
        baseline_results.sigma_total = sqrt(baseline_results.sigma_std.^2 + baseline_results.variance_std.^2);
        baseline_results.std_pred_density_after_log = baseline_results.sigma_total*Ystd; 

        % Using MATLAB's inbuilt function for base-10 lognormal distribution:
        baseline_results_using_lognormal.pd = makedist('Lognormal','mu',baseline_results.pred_density_after_log * log(10), ...
            'sigma',baseline_results.std_pred_density_after_log * log(10));
        baseline_results_using_lognormal.density_mean_physical_final = mean(baseline_results_using_lognormal.pd);
        baseline_results_using_lognormal.density_std_physical_final = std(baseline_results_using_lognormal.pd);
        baseline_results_using_lognormal.density_var_physical_final = baseline_results_using_lognormal.density_std_physical_final.^2;
        
        % Final results for the sake of plotting:
        final_results.input_variable = mc_data.xdata1_mc_lognormal(:, sw_index);
        final_results.pred_density_after_log = mc_results.pred_density_after_log;
        final_results.std_pred_density_after_log = mc_results.std_pred_density_after_log;        
        %%% Please complete the rest %%%%

        %%%%% Presentation of results %%%%%    
        display_result_space = 'physical_space';
        % display_result_space = 'logarithmic_space';
    

        %%%  Alternate figure  %%%%%

        % --- Single figure with 2x3 grid
        figure('Color','w','Name','MC Results Histograms (2x3)');
        % tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
        tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

        nexttile;
        data1 = mc_data.xdata1_mc_lognormal(:, sw_index);
        data2 = ut_data.xdata1_ut(:, sw_index);
        mu1 = mean(data1, 'omitnan');
        
        histogram(data1, nbins, 'Normalization','percentage', ...
            'FaceColor',[0.2 0.45 0.85], 'EdgeColor','none', 'FaceAlpha',0.85);
        grid on;
        hold on;
        xline(mu1, '--r', 'LineWidth',1.5);
        xline(data2(1), '--b', 'LineWidth',1);
        xline(data2(2), '-b', 'LineWidth',1);
        xline(data2(3), '-b', 'LineWidth',1);   
        legend('input MC data', 'Mean MC data', 'UT sigma point 1', ...
            'UT sigma point 2', 'UT sigma point 3')        
        xlabel('$F_{10.7}$ index','Interpreter','latex');
        ylabel('% of the total samples');
        title('Histogram of input lognormal $F_{10.7}$ and its UT sigma points','Interpreter','latex');
        ax = gca;
        ax.FontName = 'Times New Roman';
        ax.FontSize = fig_font_size;         
    
        % 2) Distribution of the output mean:
        nexttile;
        histogram(final_results.pred_density_after_log, nbins, 'Normalization','percentage', ...
            'FaceColor',[0.2 0.45 0.85], 'EdgeColor','none', 'FaceAlpha',0.85);
        grid on;
        title('Histogram of $\mu_{\hat{z}}$','Interpreter','latex');
        xlabel('$\mu_{\hat{z}}$','Interpreter','latex');
        ylabel('% of the total samples');
        hold on;
        xline(mean(final_results.pred_density_after_log,'omitnan'), '--r', 'LineWidth',1.5);
        hold off;
        ax = gca;
        ax.FontName = 'Times New Roman';
        ax.FontSize = fig_font_size;
    
        % 3) Distribution of the output standard deviation:
        nexttile;
        histogram(final_results.std_pred_density_after_log, nbins, 'Normalization','percentage', ...
            'FaceColor',[0.2 0.45 0.85], 'EdgeColor','none', 'FaceAlpha',0.85);
        grid on;
        title('Histogram of $\sigma_{\hat{z}}$','Interpreter','latex');
        xlabel('$\sigma_{\hat{z}}$','Interpreter','latex');
        ylabel('% of the total samples');
        hold on;
        xline(mean(final_results.std_pred_density_after_log,'omitnan'), '--r', 'LineWidth',1.5);
        hold off;
        ax = gca;
        ax.FontName = 'Times New Roman';
        ax.FontSize = fig_font_size;
    
        % 4) Distribution of the consolidated MC output:
        nexttile;

        if isequal(display_result_space, 'logarithmic')
    
            % In the logarithmic space:
            histogram(all_samples, nbins, 'Normalization','pdf', ...
                'FaceColor',[0.2 0.45 0.85], 'EdgeColor','none', 'FaceAlpha',0.85);
            grid on; hold on;
    
            % Mean line (choose your preferred color)
            hMean = xline(mu_hat, '--r', 'LineWidth', 1.5, 'DisplayName','Consolidated MC mean');
    
            % Use the exact same color as the mean for ±1σ lines, dash-dot style
            meanColor = get(hMean,'Color');
            hM1 = xline(mu_hat - sigma_hat, '-', 'LineWidth', 1.5, 'Color', meanColor, 'DisplayName','-1\sigma');
            hP1 = xline(mu_hat + sigma_hat, '-', 'LineWidth', 1.5, 'Color', meanColor, 'DisplayName','+1\sigma');
    
            % Mean line for unscented transformation:
            hMean_ut = xline(ut_results.pred_density_after_log, '--k', ...
                'LineWidth', 1.5, 'DisplayName','UT mean');
    
            % Use the exact same color as the mean for ±1σ lines, dash-dot style
            meanColor_ut = get(hMean_ut,'Color');
            hM1_ut = xline(ut_results.pred_density_after_log - ut_results.std_pred_density_after_log, ...
                '-', 'LineWidth', 1.5, 'Color', meanColor_ut, 'DisplayName','-1\sigma UT');
            hP1_ut = xline(ut_results.pred_density_after_log + ut_results.std_pred_density_after_log, ...
                '-', 'LineWidth', 1.5, 'Color', meanColor_ut, 'DisplayName','+1\sigma UT');
    
            % Title / labels
            title(sprintf(['$\rho_{pred}$ from consolidated MC and UT'], ...
                  tblSpaceWeather.Report_case(ii)), 'Interpreter','latex');
            xlabel('$\log_{10}[Normalized\, \rho_{pred}]\ \mathrm{(kg/m^3)}$','Interpreter','latex'); 
            ylabel('PDF');
            legend([hMean hM1 hP1 hMean_ut hM1_ut hP1_ut], {'Consolidated MC mean','-1\sigma MC','+1\sigma MC', ...
               'UT mean','-1\sigma UT','+1\sigma UT' }, 'Location','best');
    
            ax = gca; ax.FontName = 'Times New Roman'; ax.FontSize = fig_font_size;
            hold off;

        else
    
            % In the physical space:
            histogram(mc_results_using_lognormal.all_samples, nbins, 'Normalization','percentage', ...
                'FaceColor',[0.2 0.45 0.85], 'EdgeColor','none', 'FaceAlpha',0.85);
            grid on; hold on;
        
            % Mean line (choose your preferred color)
            mu_hat_physical = mc_results_using_lognormal.density_mean_physical_final;
            hMean = xline(mu_hat_physical, '--r', 'LineWidth', 1.5, 'DisplayName','Consolidated MC mean');
        
            % Use the exact same color as the mean for ±1σ lines, dash-dot style
            meanColor = get(hMean,'Color');
            mu_hat_minus_sigma_phy = mc_results_using_lognormal.density_mean_physical_final - ...
            mc_results_using_lognormal.density_std_physical_final;
            mu_hat_plus_sigma_phy = mc_results_using_lognormal.density_mean_physical_final + ...
            mc_results_using_lognormal.density_std_physical_final;
            hM1 = xline(mu_hat_minus_sigma_phy, '-', 'LineWidth', 1.5, 'Color', meanColor, 'DisplayName','-1\sigma');
            hP1 = xline(mu_hat_plus_sigma_phy, '-', 'LineWidth', 1.5, 'Color', meanColor, 'DisplayName','+1\sigma');
        
            % Mean line for unscented transformation:
            ut_mean_physical = ut_results_using_lognormal.density_mean_physical_final;
            hMean_ut = xline(ut_mean_physical, '--k', ...
                'LineWidth', 1.5, 'DisplayName','UT mean');

            % Use the exact same color as the mean for ±1σ lines, dash-dot style
            meanColor_ut = get(hMean_ut,'Color');
            ut_mean_minus_sigma_phy = ut_results_using_lognormal.density_mean_physical_final - ...
                ut_results_using_lognormal.density_std_physical_final;
            ut_mean_plus_sigma_phy = ut_results_using_lognormal.density_mean_physical_final + ...
                ut_results_using_lognormal.density_std_physical_final;
            hM1_ut = xline(ut_mean_minus_sigma_phy, ...
                '-', 'LineWidth', 1.5, 'Color', meanColor_ut, 'DisplayName','-1\sigma UT');
            hP1_ut = xline(ut_mean_plus_sigma_phy, ...
                '-', 'LineWidth', 1.5, 'Color', meanColor_ut, 'DisplayName','+1\sigma UT');
        
            % Mean line for the baseline (i.e. without input uncertainty):
            hMean_baseline = xline(baseline_results_using_lognormal.density_mean_physical_final, '--g', ...
                'LineWidth', 1.5, 'DisplayName','UT mean');

            % Use the exact same color as the mean for ±1σ lines, dash-dot style
            meanColor_baseline = get(hMean_baseline,'Color');
            baseline_results_using_lognormal.mean_minus_sigma_phy = baseline_results_using_lognormal.density_mean_physical_final - ...
            baseline_results_using_lognormal.density_std_physical_final;
            baseline_results_using_lognormal.mean_plus_sigma_phy  = baseline_results_using_lognormal.density_mean_physical_final + ...
            baseline_results_using_lognormal.density_std_physical_final;
            hM1_baseline = xline(baseline_results_using_lognormal.mean_minus_sigma_phy, ...
                '-', 'LineWidth', 1.5, 'Color', meanColor_baseline, 'DisplayName','-1\sigma UT');
            hP1_baseline = xline(baseline_results_using_lognormal.mean_plus_sigma_phy, ...
                '-', 'LineWidth', 1.5, 'Color', meanColor_baseline, 'DisplayName','+1\sigma UT'); 

            % Title / labels
            title('$\rho_{pred}$ from consolidated MC and UT', 'Interpreter','latex');
            xlabel('$\rho_{pred}\ \mathrm{(kg/m^3)}$','Interpreter','latex'); 
            ylabel('% of the total samples');
            legend([hMean hM1 hP1 hMean_ut hM1_ut hP1_ut hMean_baseline ...
                hM1_baseline hP1_baseline], {'Consolidated MC mean','-1\sigma MC','+1\sigma MC', ...
               'UT mean','-1\sigma UT','+1\sigma UT', 'Baseline mean', ...
               '-1\sigma baseline','+1\sigma baseline' }, 'Location','best');
        
            ax = gca; ax.FontName = 'Times New Roman'; ax.FontSize = fig_font_size;
            hold off;    
        end


        % % 5) Input F10.7 vs mean of predicted density
        % nexttile;
        % plot(final_results.input_variable, final_results.pred_density_after_log, 'b.')
        % grid on
        % hold on
        % xlabel(sprintf('$F_{10.7}$ index'),'Interpreter','latex');
        % ylabel('$\mu_{\hat{z}}$','Interpreter','latex');
        % title('Input $F_{10.7}$ index (vs) $\mu_{\hat{z}}$', ...
        %     'Interpreter','latex')
        % ax = gca; ax.FontName = 'Times New Roman'; ax.FontSize = 15;
        % 
        % % 6) Input F10.7 vs std of predicted density
        % nexttile;
        % plot(final_results.input_variable, final_results.std_pred_density_after_log, 'b.')
        % grid on
        % hold on
        % xlabel(sprintf('$F_{10.7}$ index'),'Interpreter','latex');
        % ylabel('$\sigma_{\hat{z}}$','Interpreter','latex');
        % title('Input $F_{10.7}$ index (vs) $\sigma_{\hat{z}}$', ...
        %     'Interpreter','latex')
        % ax = gca; ax.FontName = 'Times New Roman'; ax.FontSize = 15;
        % 
        % sgtitle(sprintf(['Comparisons b/w MC and UT with NEW post-processing' ...
        %     '(case %d %s): %d samples from %d distributions for consolidated MC for noise ratio %0.2f'],  ...
        %     tblSpaceWeather.Report_case(ii), tblSpaceWeather.F107_indication(ii), ...
        %     n_samples_consolidated, n_samples_MC, noise_ratio_current), ...
        %     'Interpreter','latex','FontSize', 15)        

        if isequal(display_result_space, 'logarithmic')
    
            disp('Consolidated MC mean:')
            disp(mu_hat)
    
            disp('UT mean:') 
            disp(ut_results.pred_density_after_log)
    
            disp('Consolidated MC std.:')
            disp(sigma_hat)
    
            disp('UT std.:')
            disp(ut_results.std_pred_density_after_log)

        else
    
            % Display results in the physical space:
            disp('Noise ratio')
            disp(noise_ratio_current)

            disp('Consolidated MC mean in physical space:')
            disp(mc_results_using_lognormal.density_mean_physical_final)

            disp('Baseline mean in physical space:')
            disp(baseline_results_using_lognormal.density_mean_physical_final)            
        
            disp('UT mean in physical space:') 
            disp(ut_results_using_lognormal.density_mean_physical_final)
        
            disp('Relative error in the mean density in physical space:')
            disp(abs((ut_results_using_lognormal.density_mean_physical_final/mc_results_using_lognormal.density_mean_physical_final) - 1))
        
            disp('Consolidated MC std. in physical space:')
            disp(mc_results_using_lognormal.density_std_physical_final)

            disp('Baseline std. in physical space:')
            disp(baseline_results_using_lognormal.density_std_physical_final)             
        
            disp('UT std. in physical space:')
            disp(ut_results_using_lognormal.density_std_physical_final)   
        
            disp('Relative error in the std. density in physical space:')
            disp(abs((ut_results_using_lognormal.density_std_physical_final/mc_results_using_lognormal.density_std_physical_final) - 1))      

        end

    end

end

% %%%% For checking purposes %%%%%
% [sorted_input, sort_idx] = sort(final_results.input_variable);
% sorted_mu = final_results.pred_density_after_log(sort_idx);
% 
% figure
% plot(sorted_input, sorted_mu, 'b.')  % Use line with points
% grid on
% xlabel(sprintf('$F_{10.7}$ index'),'Interpreter','latex');
% ylabel('$\mu_{\hat{z}}$','Interpreter','latex');
% title('Input $F_{10.7}$ index (vs) $\mu_{\hat{z}}$', 'Interpreter','latex')
% ax = gca; ax.FontName = 'Times New Roman'; ax.FontSize = 15;


%% Local functions:

function [mu_y, P_yy] = ut_recombine_outputs(Y, Wm, Wc)
% UT_RECOMBINE_OUTPUTS  Mean & covariance from propagated UT sigma outputs.
%   [mu_y, P_yy] = ut_recombine_outputs(Y, Wm, Wc)
%
% Inputs
%   Y  : ny x L matrix of outputs, one column per sigma point
%        (here ny = 4, L = 3 for 1-D UT)
%   Wm : 1 x L mean weights
%   Wc : 1 x L covariance weights
%
% Outputs
%   mu_y : ny x 1 UT mean of outputs
%   P_yy : ny x ny UT covariance of outputs
%
% Notes
%   - Works for any output dimension ny and any number of sigma points L
%     as long as Wm/Wc are compatible.

    [ny, L] = size(Y);
    assert(length(Wm) == L && length(Wc) == L, 'Weights must match columns of Y');

    % Mean
    mu_y = Y * Wm(:);       % ny x 1

    % Covariance
    P_yy = zeros(ny, ny);
    for i = 1:L
        d = Y(:,i) - mu_y;
        P_yy = P_yy + Wc(i) * (d * d.');
    end
end


function [all_samples, mu_hat, sigma_hat] = draw_from_many_gaussians_loop_v2(mc_results, m)
%DRAW_FROM_MANY_GAUSSIANS_LOOP
% Sample m times from each of n Gaussians (with their own mean and std),
% stack all samples into one vector, and return the consolidated statistics.
%
% Inputs
%   mc_results.pred_density_after_log    : [n x 1] or [1 x n] vector of means
%   mc_results.std_pred_density_after_log: [n x 1] or [1 x n] vector of std devs
%   m                                    : positive integer, samples per Gaussian
%
% Outputs
%   all_samples : [m*n x 1] vector of stacked samples (m from each of n Gaussians)
%   mu_hat      : scalar, sample mean of all m*n draws (uses 'omitnan')
%   sigma_hat   : scalar, sample std (unbiased, N-1) of all m*n draws (uses 'omitnan')
%
% Notes
% - This function intentionally uses a FOR-LOOP to (1) draw m samples for
%   each Gaussian and (2) stack them in order. This preserves the notion of
%   "m samples per distribution," rather than vectorizing across all n.
% - The stacked distribution is a mixture of Gaussians (generally NOT
%   Gaussian). The returned mu_hat and sigma_hat are the empirical statistics
%   of that mixture sample set.

    % -------------------------------
    % 0) Input checks and preparation
    % -------------------------------
    if nargin < 2 || isempty(m)
        error('Provide m = samples per Gaussian (positive integer).');
    end

    assert(isfield(mc_results,'pred_density_after_log'), ...
        'mc_results.pred_density_after_log is missing.');
    assert(isfield(mc_results,'std_pred_density_after_log'), ...
        'mc_results.std_pred_density_after_log is missing.');

    mu_vec  = mc_results.pred_density_after_log(:);      % [n x 1] means
    sig_vec = mc_results.std_pred_density_after_log(:);  % [n x 1] std devs
    n = numel(mu_vec);

    assert(numel(sig_vec) == n, ...
        'Lengths of means and std devs must match (n means and n stds).');
    assert(all(isfinite(mu_vec)), 'Means must be finite.');
    assert(all(isfinite(sig_vec) & sig_vec >= 0), ...
        'Standard deviations must be finite and nonnegative.');
    assert(isscalar(m) && m == floor(m) && m > 0, ...
        'm must be a positive integer.');

    % ------------------------------------------
    % 1) For each (mu_i, sigma_i), draw m samples
    % 2) Stack the samples via a for-loop
    % ------------------------------------------
    all_samples = zeros(m*n, 1);  % preallocate for speed
    for i = 1:n
        mu_i  = mu_vec(i);
        sig_i = sig_vec(i);

        % Draw m samples from N(mu_i, sig_i^2)
        samples_i = mvnrnd(mu_i, sig_i^2, m);

        % Put these m samples into their block in all_samples
        idx = (i-1)*m + (1:m);
        all_samples(idx) = samples_i;
    end

    % -----------------------------------------------------
    % 3) Compute consolidated mean and standard deviation
    % -----------------------------------------------------
    mu_hat    = mean(all_samples, 'omitnan');
    sigma_hat = std(all_samples, 0, 'omitnan');   % 0 => unbiased (N-1)

    % -----------------------------------------------------
    % 4) Return outputs: all_samples, mu_hat, sigma_hat
    % (done via function outputs)
    % -----------------------------------------------------
end

function x_all = sample_positive_lognormal(means, stds, m_per_pair, rngSeed)
% means, stds: column vectors (n x 1) of physical-domain mean & std
% m_per_pair: number of samples per (mean,std)
% Returns x_all: (m_per_pair*n) x 1 vector of positive samples

    if nargin>3, rng(rngSeed); end
    n = numel(means);
    x_all = zeros(m_per_pair*n,1);
    idx = 1;
    
    for i = 1:n
        m = means(i); s = stds(i);
        sigma2 = log(1 + (s/m)^2);
        mu_e   = log(m) - 0.5*sigma2;
        x_all(idx:idx+m_per_pair-1) = lognrnd(mu_e, sqrt(sigma2), [m_per_pair, 1]);
        idx = idx + m_per_pair;
    end
end