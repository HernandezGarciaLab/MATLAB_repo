function im = recon_vsasl3dflex(varargin)
% function im = recon_vsasl3dflex(varargin)
%
% Part of umasl project by Luis Hernandez-Garcia and David Frey
% @ University of Michigan 2022
%
% Description: Function to reconstruct images from vsasl3dflex ASL
%   sequence using 3D NUFFT with Pipe & Menon Density compensation
%
%
% Notes:
%   - if output is returned, nii files will not be saved to conserve space
%       and prevent overwriting, and vice verse for when output is not
%       returned (see 'im' under 'Function output')
%   - default values in help message may not be up to date - check defaults
%       structure under the function header
%
% Path dependencies:
%   - matlab default path
%       - can be restored by typing 'restoredefaultpath'
%   - umasl
%       - github: fmrifrey/umasl
%       - umasl/matlab/ and subdirectories must be in current path
%   - mirt (matlab version)
%       - github: JeffFessler/mirt
%       - mirt setup must have successfully ran
%
% Variable input arguments (type 'help varargin' for usage info):
%   - 'raw':
%       - raw data structure
%       - complex double/float array of dimensions (nframes x ndat x
%           nleaves x nslices x ncoils)
%       - if either raw or info is left empty, function will read data from
%           the Pfile instead (see 'pfile' or type 'help readpfile' for
%           more information)
%       - default is empty (reads raw data from file)
%   - 'info'
%       - pfile information structure
%       - structure containing fields: ndat, nleaves, nframes, nslices,
%           ncoils, tr, te, dim, fov, slthick
%       - type 'help readpfile' for more information on info structure
%           format and field data types
%       - if either raw or info is left empty, function will read data from
%           the Pfile instead (see 'pfile' or type 'help readpfile' for
%           more information)
%       - default is empty (reads info from file)
%   - 'pfile'
%       - search string for pfile
%       - string containing search path for desired pfile
%       - will only be used if either 'raw' or 'info' is left blank
%       - type 'help readpfile' for more information
%       - default is 'P*.7' (uses first pfile from directory to read in
%           raw data and info)
%   - 'tol'
%       - kspace distance roundoff tolerance
%       - double/float
%       - corrects issues in generating ktraj_cart.txt float values
%       - default is 1e-5
%   - 'itrmax'
%       - maximum number of iterations for iterative operations
%       - integer describing number of iterations
%       - used in iterative psf recon to optimize density compensation
%           (Pipe & Menon method)
%       - default is 15 (optimal value described by Pipe & Menon)
%   - 'frames'
%       - frames in timeseries to recon
%       - integer array containing specific frames to recon
%       - if 'all' is passed, all frames will be reconned
%       - default is 'all'
%   - 'ndel'
%       - gradient sample delay compensation
%       - integer describing number of samples to shift signal by in each
%           echo
%       - default is 0
%   - 'nramp'
%       - number of ramp points in spiral to delete from data
%       - integer describing number of points in ramp
%       - if 'auto' is passed, ramp points will be determined based on
%           trajectory envelope
%       - make sure nramp > ndel to clip out misplaced data
%       - default is 'auto'
%   - 'pdorder'
%       - order of least squares polynomial fit for phase drift
%           compensation
%       - integer describing highest order of polynomial for lsq fit
%       - if -1 is passed, no phase detrending will be done
%       - default is -1
%   - 'nworkers'
%       - number of workers to use in parallel pool
%       - integer describing number of workers
%       - default is output of feature('numcores') (number of available
%           cores)
%   - 'scaleoutput'
%       - option to scale nii files to full dynamic range
%       - boolean integer (0 or 1) to use or not
%       - type 'help writenii' for more information
%       - default is 1
%
% Function output:
%   - im:
%       - output timeseries image
%       - complex array of image dimension
%       - if im is not returned, timeseries image and psf will be saved to
%           nii files
%       - if im is returned, timeseries image and psf will not be saved to
%           nii files
%

    % Define default arguments
    defaults = struct(...
        'raw',          [], ... % Raw data
        'info',         [], ... % Info structure
        'pfile',        'P*.7', ... % Search string for Pfile
        'tol',          1e-5, ... % Kspace distance tolerance
        'itrmax',       15, ... % Max number of iterations for IR
        'frames',       'all', ... % Frames to recon
        'ndel',         0, ... % Gradient sample delay
        'nramp',        [], ... % Number of ramp points in spiral traj
        'pdorder',      -1, ... % Order of phase detrending poly fit
        'nworkers',     feature('numcores'), ... % Number of workers to use in parpool
        'scaleoutput',  1 ... % Option to scale output to full dynamic range
        );

    % Start timer
    t = tic;
    
%% Set up recon
    % Parse through variable inputs using matlab's built-in input parser
    args = vararginparser(defaults,varargin{:});
    
    % Get raw data and info
    if (isempty(args.raw) || isempty(args.info)) % If reading from Pfile
        [raw,info] = readpfile(args.pfile);
    else % If raw/scaninfo is user specified
        raw = args.raw;
        info = args.info;
    end
    
    % Determine frames to recon 
    if strcmpi(args.frames,'all')
        args.frames = 1:info.nframes;
    elseif max(args.frames) > info.nframes
        fprintf('\nWarning: Frames array exceeds total number of frames');
        fprintf('\n\t--> reconning all frames');
        args.frames = 1:info.nframes;
    end
    
    % Get un-transformed kspace trajectory and view transformations
    ks_0 = load('ktraj_cart.txt');
    kviews = load('kviews.txt');
    
    % Start parpool
    if args.nworkers>0 && isempty(gcp('nocreate'))
        fprintf('\n');
        parpool(args.nworkers);
    end
    
%% Process Trajectory
    % Determine nramp and navpoints from 1st platter envelope
    navpts = find(vecnorm(ks_0,2,2) < args.tol);
    navpts(navpts > info.ndat*3/4) = []; % reject navpts in first 1/4
    navpts(navpts < info.ndat*1/4) = []; % reject navpts in last 1/4
    if isempty(args.nramp)
        [~,args.nramp] = max(vecnorm(ks_0(1:round(info.ndat/2),:),2,2));
        args.nramp = args.nramp + 2; % add a small sample buffer
    end
    
    % Transform trajectory to entire trajectory
    ks = zeros(info.ndat, 3, info.nleaves, info.nslices);
    for leafn = 1:info.nleaves
        for slicen = 1:info.nslices
            % Determine indexed transformation matrix
            viewn = (leafn - 1) * info.nslices + slicen;
            T = reshape(kviews(viewn, end-8:end), 3, 3)';
            
            % Apply transformation to indexed kspace view
            ks(:,:,leafn,slicen) = (T*ks_0')';
        end
    end
    
%% Apply corrections/filters
    % Perform phase detrending
    if args.pdorder > -1
        raw = phasedetrend(raw,navpts,args.pdorder);
    end
    
    % Correct for gradient sample delay
    if abs(args.ndel) > 0
        raw = circshift(raw,args.ndel,2);
    end
    
    % Remove ramp points
    fprintf('\nRemoving %d ramp points from data...', args.nramp);
    ks([1:args.nramp info.ndat-args.nramp:info.ndat],:,:,:) = [];
    raw(:,[1:args.nramp info.ndat-args.nramp:info.ndat],:,:,:) = [];

%% Set up NUFFT
    % Extract dim and fov from info so info isn't broadcasted to parfor
    dim = info.dim;
    fov = info.fov;

    % Reshape and scale k from -pi to pi
    omega = 2*pi * fov / dim * ...
        [reshape(ks(:,1,:,:),[],1), ...
        reshape(ks(:,2,:,:),[],1), ...
        reshape(ks(:,3,:,:),[],1)];
    
    % Create NUFFT object
    nufft_args = {dim*ones(1,3),...
        6*ones(1,3),...
        2*dim*ones(1,3),...
        dim/2*ones(1,3),...
        'table',...
        2^10,...
        'minmax:kb'};
    G = Gnufft(true(dim*ones(1,3)), [{omega}, nufft_args(:)']);
    
    % Create density compensation using pipe algorithm
    dcf = pipedcf(G,args.itrmax);
    
%% Reconstruct image using NUFFT
    % Reconstruct point spread function
    psf = reshape(G' * ones(numel(ks(:,1,:,:)),1), dim*ones(1,3));
    
    % Initialize output array and progress string
    im = zeros([dim*ones(1,3),info.ncoils,length(args.frames)]);
    fprintf('\nReconning image... ');
    msg_fprog = '';
    
    % Loop through frames
    for framen = args.frames
        
        % Print progress
        if ~isempty(msg_fprog)
            fprintf(repmat('\b',1,length(msg_fprog)));
        end
        msg_fprog = sprintf('(frame %d/%d)',framen,info.nframes);
        fprintf(msg_fprog);
        
        % Loop through coils
        parfor (coiln = 1:info.ncoils, args.nworkers)
            % Recon using adjoint NUFFT operation
            data = reshape(raw(framen,:,:,:,coiln),[],1);
            im(:,:,:,coiln,framen) = reshape(G' * (dcf.*data), dim*ones(1,3));
        end
        
    end
    
    % Combine coils using RMS
    im = sqrt( squeeze( mean( im.^2, 4) ) );
    
    % Save results that aren't returned to nifti:
    if nargout < 1
        % Save timeseries
        writenii('./timeseries.nii', abs(im), ...
            fov*ones(1,3), info.tr, args.scaleoutput);
        fprintf('\nTimeseries saved to timeseries.nii');
        
        % Save point spread function
        writenii('./psf.nii', abs(psf), ...
            fov*ones(1,3), info.tr, args.scaleoutput);
        fprintf('\nPoint spread function saved to psf.nii');
        
        % Clear im so it won't be returned
        clear im;
    else
        fprintf('\nImages will not be saved to file since timeseries is returned');
    end
    
    % Save and print elapsed time
    t = toc(t);
    fprintf('\nRecon complete. Total elapsed time: %.2fs\n',t);
    
end

%% phasedetrend function definition
function raw_corr = phasedetrend(raw,navpts,pdorder)

    % Get dimensions
    nframes = size(raw,1);
    ndat = size(raw,2);
    nleaves = size(raw,3);
    nslices = size(raw,4);
    ncoils = size(raw,5);
    
    % Define design matrix for lsq poly fit
    A = (navpts(:) - round(ndat/2)).^(pdorder:-1:0);
    
    % Intitialize fits array
    fits = zeros(nframes,ndat,nleaves,nslices,ncoils);
    
    % Loop through all echos
    for framen = 1:nframes
        for leafn = 1:nleaves
            for slicen = 1:nslices
                for coiln = 1:ncoils
                    % Fit polynomial to center phase of indexed echo
                    echo = raw(framen,:,leafn,slicen,coiln);
                    y = unwrap(angle(echo(navpts)));
                    betas = pinv(A)*y(:);
                    fits(framen,:,leafn,slicen,coiln) = ...
                        (1:ndat)'.^(pdorder:-1:0) * betas;
                end
            end
        end
    end
    
    % Correct echo by subtracting out fits from phase
    raw_corr = raw.*exp(-1i*fits);

end

%% pipedcf function definition
function Wi = pipedcf(G,itrmax)

    % Initialize weights to 1 (psf)
    Wi = ones(size(G,1),1);
    fprintf('\nCreating density compensation... ');
    msg_iprog = '';
    
    % Loop through iterations
    for itr = 1:itrmax
        
        % Print progress
        if ~isempty(msg_iprog)
            fprintf(repmat('\b',1,length(msg_iprog)));
        end
        msg_iprog = sprintf('(itr %d/%d)',itr,itrmax);
        fprintf(msg_iprog);
        
        % Pipe algorithm: W_{i+1} = W_{i} / (G * (G' * W_{i}))
        d = real( G.arg.st.interp_table(G.arg.st, ...
            G.arg.st.interp_table_adj(G.arg.st, Wi) ) );
        Wi = Wi ./ d;
        
    end
    
end