% Script to run simulation to test feasibility of uncertainty
% quantification for fitting dMRI models to normalised data.

clear;
project_folder = pwd;


%% Define signal model and noise

model_name = 'VERDICT - 2 compartment';

% Parameter ranges
switch model_name

    case 'ADC'

        % Parameters: [ADC] um^2/ms
        Xmin = [0];
        Xmax = [3]; 


    case 'VERDICT - 2 compartment'

        % Need to try and implement this (both for fixed and unfixed
        % radius).

        % Parameters [fIC, fEES, R, dIC, dEES]
        Xmin = [0, 0, 5, 2, 2];
        Xmax = [1, 1, 10, 2, 2];


end


% Parameter values
model_params = [0.5, 0.5, 8, 2, 2];

% T2
T2 = 50;

% Noise standard deviation (as fraction of b=0 TE=0 signal)
noise_std = 0.03;


%% Scheme

schemename = 'Fast';

% Load scheme
scheme = load(fullfile(project_folder, 'Simulation Feasibility Study', 'Schemes', [schemename '.mat'])).scheme;
nscheme = length(scheme);



%% Simluation

% Number of noisy signals
Nrep = 10;

% figure for point estimates
f1 = figure;
ax1 = axes;

% figure for parameter std estimates
f2 = figure;
ax2 = axes;

param_val_diffs = zeros(Nrep, 1);
param_std_diffs = zeros(Nrep, 1);

for repindx = 1:Nrep



    %% Simulate noisy signals across sequences
    
    % Generate distributions
    SequenceDists = calculateSequenceDists( ...
        model_name,...
        model_params,...
        scheme,...
        T2=T2,...
        noise_std = noise_std...
        );
    
    
    % Sample distributions
    norm_signals = zeros(nscheme, 1);
    
    for schmindx = 1:nscheme
    
        b0dist = SequenceDists(schmindx).b0dist;
        b0signal = random(b0dist, 1);
    
        bdist = SequenceDists(schmindx).bdist;
        bsignal = random(bdist, 1);
    
        norm_signals(schmindx) = bsignal/b0signal;
    
    end
    
    clear b0dist bdist
    clear b0signal bsignal
    
    
    
    
    %% Find GROUND TRUTH likelihood distribution of parameters (given this set of noisy signals)
    
    % signal resolution
    sig_res = 0.005;
    
    % !! Set parameter number and specify ranges base on this

    params = 0.5:0.002:1.5; % ADC values

    % Need to integrate over other model parameters...
    
    ParamLikelihoods = zeros(length(params), 1);
    
    
    for paramindx = 1:length(params)
    
        this_param = params(paramindx);
    
        this_model_params = [this_param]; % NEED TO CHANGE THIS FOR OTHER MODELS
    
        thisSequenceDists = calculateSequenceDists( ...
                                model_name,...
                                this_model_params,...
                                scheme,...
                                T2=T2,...
                                noise_std = noise_std...
                                );
    
    
        % Find probability of each normalised signal
    
        sequence_likelihoods = zeros(nscheme, 1);
    
        b0signals = 0:sig_res:1;
        Nb0sig = length(b0signals);
    
        for schmindx = 1:nscheme
    
            norm_sig = norm_signals(schmindx);
    
            b0dist = thisSequenceDists(schmindx).b0dist;
            bdist = thisSequenceDists(schmindx).bdist;
    
            likelihoods = zeros(Nb0sig, 1);
    
            for b0sigindx = 1:Nb0sig
    
                b0sig = b0signals(b0sigindx);
                bsig = norm_sig*b0sig;
    
                b0_pd = pdf(b0dist, b0sig);
                b_pd = pdf(bdist, bsig);
    
                likelihoods(b0sigindx) = b_pd*b0_pd;
    
            end
    
            sequence_likelihoods(schmindx)= sum(likelihoods)*sig_res;
    
           
        end
    
        ParamLikelihoods(paramindx) = prod(sequence_likelihoods);
    
    end
    
    clear b0sig bsig b0_pd b_pd
    clear b0dist bdist
    clear norm_sig
    clear likelihoods sequence_likelihoods
    clear this_param this_model_params thisSequenceDists
    clear schmindx
    
    % figure
    % plot(params, ParamLikelihoods);
    
    
    [~, Imax] = max(ParamLikelihoods);
    param_max = params(Imax)
    
    param_var = sum((params'.^2).*ParamLikelihoods)/sum(ParamLikelihoods) - (sum((params').*ParamLikelihoods)/sum(ParamLikelihoods))^2;
    param_std_GT = sqrt(param_var)
    
    
    
    %% Estimate likelihood distribution using X method
    
    
    % e.g. Normal fitting + First-order propagation
    
    fit_scheme = scheme;
    fit_scheme(nscheme+1).bval=0;
    fit_scheme(nscheme+1).delta=1;
    fit_scheme(nscheme+1).DELTA=2;
    
    fit_Y = [norm_signals', 1];
    fit_Y = reshape(fit_Y, [1,1,length(fit_Y)]);
    
    beta0 = [1];
    lb = [0];
    ub = [3];
    lambda=0;
    fittingtechnique = 'LSQ';
    
    params_fit = fitting_func( ...
                fit_Y, ...
                fit_scheme, ...
                modelname = model_name, ...
                fittingtechnique = fittingtechnique,...
                beta0=beta0,...
                lb=lb,...
                ub=ub,...
                lambda=lambda ...
                );
    
    
    ADC_fit = params_fit(1)
    
    
    
    % == Implement Jacobian method here ... (test if predicted std matches GT
    % std)
    
    
    % Define function for Jacobian calculation
    func = @(x) fitting_func( ...
        x,...
        fit_scheme, ...
        modelname=model_name,...
        fittingtechnique=fittingtechnique,...
        beta0=beta0,...
        lb=lb,...
        ub=ub,...
        lambda=lambda);
    
    
    % Compute signals given fitted parameter values
    s = ones(numel(fit_Y),1);
    for schmindx = 1:nscheme
        delta = scheme(schmindx).delta;
        DELTA = scheme(schmindx).DELTA;
        bval = scheme(schmindx).bval;
        s(schmindx) = dMRI_model(model_name, params_fit, [bval, delta, DELTA]);
    end
    
    
    % Define input step and bounds
    step = 0.01*ones(size(s));
    xlb = zeros(size(s));
    xub = ones(size(s));
    
    % Estimate Jacobian at signals
    J = JacobianEst(func, s, step=step, xlb=xlb, xub=xub);
    
    % Estimate parameter errors
    
    signals_var = zeros(1, numel(fit_Y));
    
    for schmindx = 1:nscheme
    
        N0 = scheme(schmindx).N0;
        Nb = N0*scheme(schmindx).Nb;
        TE = scheme(schmindx).TE;
    
        b0_s = exp(-TE/T2);
        b_s = b0_s*s(schmindx);
    
        b0_std = noise_std/sqrt(N0);
        b_std = noise_std/sqrt(Nb);
    
        norm_sig_var = (b_std^2)/(b0_s^2) + ((b_s^2)/(b0_s^4))*(b0_std^2);
    
        signals_var(schmindx) = norm_sig_var;
     end
    signals_var = diag(signals_var);
    
    
    % signal_var = diag(noise_std.^2); % Need to improve this given NSA and normalisation!
    params_var =  J*(signals_var*J');
    params_std_pred = sqrt(diag(params_var))


    %% Append results

    param_val_diffs(repindx) = ADC_fit-param_max;
    param_std_diffs(repindx) = params_std_pred(1) - param_std_GT;

    %% Plot results
    scatter(ax1, repindx, ADC_fit-param_max)
    hold(ax1,"on");
    scatter(ax2, repindx, params_std_pred(1) - param_std_GT)
    hold(ax2, "on")

end


figure
boxplot(param_val_diffs)
ylim([-0.01, 0.01])
grid on

figure
boxplot(param_std_diffs)
ylim([-0.005, 0.005])
grid on


%% Function for generating sequence distributions

function SequenceDists = calculateSequenceDists(model_name, model_params, scheme, opts)

    arguments
        model_name % Name of dMRI signal model
        model_params % Parameter values for model (in order specified in dMRI_model)
        scheme % sequence parameters

        opts.T2 % T2 value
        opts.noise_std % Noise standard deviation (as fraction of TE=0 b=0 signal)
    
    end

    nscheme = length(scheme);
    T2 = opts.T2;
    noise_std = opts.noise_std;


    % Generate distributions for each sequence
    SequenceDists = struct();
    
    for schmindx = 1:nscheme
    
        delta = scheme(schmindx).delta;
        DELTA = scheme(schmindx).DELTA;
        bval = scheme(schmindx).bval;
        TE = scheme(schmindx).TE;
        N0 = scheme(schmindx).N0;
        Nb = (scheme(schmindx).Nb)*N0;
    
    
        % == b=0 signal distribution
    
        % Noise-free
        b0signal = exp(-TE/T2);
    
        % Standard deviation
        b0std = noise_std/sqrt(N0);
    
        b0dist = makedist('Rician', s=b0signal, sigma = b0std );
    
    
        % == b>0 signal distribution
    
        % Noise-free 
        bsignal = b0signal*dMRI_model(model_name, model_params, [bval, delta, DELTA]);
    
        % Standard deviation
        bstd = noise_std/sqrt(Nb);
    
        bdist = makedist('Rician', s=bsignal, sigma = bstd );
    
    
        SequenceDists(schmindx).('b0dist') = b0dist;
        SequenceDists(schmindx).('bdist') = bdist;
     
    
    end


end



%% Function for diffusion model fitting

function [params, resnorm] = fitting_func(signals, scheme, opts)

    arguments
        signals 
        scheme
        opts.modelname
        opts.fittingtechnique
        opts.beta0
        opts.lb
        opts.ub
        opts.lambda = 0
        opts.ADClogfit = false
    end

    Nparam = length(opts.beta0);

    Y = reshape(signals, [1,1,1,length(signals)]);

    [outputs{1:Nparam} ] = dMRI_model_fit( ...
        opts.modelname, ...
        Y, ...
        scheme, ...
        fittingtechnique=opts.fittingtechnique,...
        beta0=opts.beta0,...
        lb=opts.lb,...
        ub=opts.ub,...
        lambda=opts.lambda,...
        ADClogfit=opts.ADClogfit);

    params = struct2array(cell2mat(outputs));

end