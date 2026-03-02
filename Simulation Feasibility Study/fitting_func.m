%% Function for diffusion model fitting

function [params, resnorm, jacobians] = fitting_func(signals, scheme, opts)

    arguments
        signals 
        scheme
        opts.modelname
        opts.fittingtechnique
        opts.beta0
        opts.lb
        opts.ub
        opts.lambda = 0
        opts.return_jacobian = false
        opts.ADClogfit = false
    end

    Nparam = length(opts.beta0);

    Y = reshape(signals, [1,1,1,length(signals)]);

    outputs = dMRI_model_fit( ...
        opts.modelname, ...
        Y, ...
        scheme, ...
        fittingtechnique=opts.fittingtechnique,...
        beta0=opts.beta0,...
        lb=opts.lb,...
        ub=opts.ub,...
        lambda=opts.lambda,...
        return_jacobian = opts.return_jacobian,...
        ADClogfit=opts.ADClogfit);

    f = fields(outputs);
    nf = length(f);
    params = zeros(nf, 1);
    for indx = 1:nf
        if strcmp(f(indx), 'Jacobians')
            jacobians = outputs.(f{indx});
            continue
        end
        params(indx) = outputs.(f{indx});
    end

    % ...
    if ~exist('resnorm')
        resnorm = 0;
    end

    if ~exist('jacobians')
        jacobians = 0;
    end


end