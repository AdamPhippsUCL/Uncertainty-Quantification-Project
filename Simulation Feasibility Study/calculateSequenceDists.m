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