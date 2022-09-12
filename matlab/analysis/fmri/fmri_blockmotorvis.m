function tscore_act = fmri_blockmotorvis(im,varargin)

    % Define default arguments
    defaults = struct(...
        'tstart',       0, ... % Start time (s) of regression
        'toff',         30, ... % Off time (s) for block stimulation
        'ton',          30, ... % On time (s) for block stimulation
        'tr',           [], ... % TR (s) (if im is not read from file)
        'fov',          [], ... % FOV (cm) (if im is not read from file)
        'show',         0, ... % Option to show regressors
        'scaleoutput',  1 ... % Option to scale output to full dynamic range
        );

    % Start timer
    t = tic;
    
    % Parse through variable inputs using matlab's built-in input parser
    args = vararginparser(defaults,varargin{:});
    
    % Define default for im
    if nargin < 1 || isempty(im)
        im = 'sub';
    elseif iscomplex(im)
        warning('Complex images are not supported, using absolute value');
        im = abs(im);
    end
    
    % If im is a nii file name, read in from file
    if ischar(im)
        [im,h] = readnii(im);
        args.tr = h.pixdim(5) * 1e-6;
        args.fov = h.dim(2:4) .* h.pixdim(2:4);
    elseif (isempty(args.tr) || isempty(args.fov))
        error('Must specify tr and fov if image is not read from file');
    end
    
    % Get dimensions and convert
    nframes = size(im,4);
    
    % Construct hi-res timing
    tend = args.tr * (nframes - 1);
    t_hires = 0:1e-3:tend;
    t_lowres = 0:args.tr:tend;
    
    % Determine stimulation signal
    stim = 1 * (t_hires >= args.tstart) .* ...
        (mod(t_hires - args.tstart, args.ton + args.toff) > args.toff);
    
    % Get hemodynamic response function using spm
    [~, p] = spm_hrf(args.tr * 1e-3);
    p(5) = 100;
    [hrf, ~] = spm_hrf(args.tr * 1e-3, p);
    
    % Create regressors
    x_flow = ones(nframes,1); % Baseline flow
    x_act_hires = conv(hrf,stim); % Activation (stim convolved with hrf)
    x_act_hires = x_act_hires(1:length(t_hires));
    x_act = interp1(t_hires, x_act_hires, t_lowres)';
    
    % Create design matrix and known vector for lsq regression 
    A = [x_flow, x_act];
    y = reshape(permute(im, [4,1:3]), h.dim(5), prod(h.dim(2:4)));
    
    % Estimate beta for lsq problem y = Ab
    beta = pinv(A) * y;
    
    % Estimate t-scores and z-scores for beta maps
    tscore = zeros(size(beta));
    for c = 1:2
        tscore(c,:) = beta(c,:) ./ std(beta(c,:),[],2);
    end
    
    % Reshape beta and t-score maps
    beta_flow = reshape(beta(1,:), h.dim(2:4));
    beta_act = reshape(beta(2,:), h.dim(2:4));
    tscore_flow = reshape(tscore(1,:),h.dim(2:4));
    tscore_act = reshape(tscore(2,:),h.dim(2:4));
    
    % Save data to files
    if nargout < 1
        % Save beta, tscore, and zscore images to file
        writenii('./beta_flow.nii',beta_flow,args.fov,args.tr,args.scaleoutput);
        fprintf('\nBeta map for baseline flow saved to beta_flow.nii');
        writenii('./beta_act.nii',beta_act,args.fov,args.tr,args.scaleoutput);
        fprintf('\nBeta map for activation saved to beta_act.nii');
        writenii('./tscore_flow.nii',tscore_flow,args.fov,args.tr,args.scaleoutput);
        fprintf('\nT-score map for baseline flow saved to tscore_flow.nii');
        writenii('./tscore_act.nii',tscore_act,args.fov,args.tr,args.scaleoutput);
        fprintf('\nT-score map for activation saved to tscore_act.nii');
        
        % Clear tscore_act so it won't be returned as ans
        clear tscore_act;
    else
        fprintf('\nImages will not be saved to file since sub image is returned');
    end
    
    % Show regressors
    if args.show
        cfigopen('fmri block motor vis analysis');
        plot(t_hires,x_act_hires,'-b'), hold on
        scatter(t_lowres,x_act,'ob'),
        plot(t_hires,ones(length(t_hires),1),'-g'),
        scatter(t_lowres,ones(length(t_lowres),1),'og'); hold off
        xlabel('time (s)');
        yticks([]);
        title('Regressors');
    end
    
    % Save and print elapsed time
    t = toc(t);
    fprintf('\nfmri block motor vis analysis completed. Elapsed time: %.2fs\n',t);
    
end

