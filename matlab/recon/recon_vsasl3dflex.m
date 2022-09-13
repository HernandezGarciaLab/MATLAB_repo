function im = recon_vsasl3dflex(varargin)
    
    % Define default arguments
    defaults = struct(...
        'raw',          [], ... % Raw data
        'info',         [], ... % Info structure
        'pfile',        'P*.7', ... % Search string for Pfile
        'tol',          1e-5, ... % Kspace distance tolerance
        'itrmax',       15, ... % Max number of iterations for IR
        'frames',       'all', ... % Frames to recon
        'ndel',         0, ... % Gradient sample delay
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
    % Determine trajectory type from kviews
    isSOS = (all(kviews(:,4)==0) && all(kviews(:,5)==0));
    
    % Determine nramp and navpoints from 1st platter envelope
    navpts = find(vecnorm(ks_0,2,2) < args.tol);
    navpts(navpts > info.ndat*3/4) = []; % reject navpts in first 1/4
    navpts(navpts < info.ndat*1/4) = []; % reject navpts in last 1/4
    [~,nramp] = max(vecnorm(ks_0(1:round(info.ndat/2),:),2,2));
    nramp = nramp + 2; % add a small sample buffer
    
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
    fprintf('\nRemoving %d ramp points from data...', nramp);
    ks([1:nramp info.ndat-nramp:info.ndat],:,:,:) = [];
    raw(:,[1:nramp info.ndat-nramp:info.ndat],:,:,:) = [];

%% Set up NUFFT
    % Determine FOV/dim/Nneighbors based on trajectory type
    fov = info.xyfov*ones(1,3);
    dim = info.xydim*ones(1,3);
    nneighbors = 6*ones(1,3);
    if isSOS
        % For SOS, fix the z fov, dim, and number of neighbors
        fov(3) = info.slthick*info.nslices;
        dim(3) = info.nleaves*info.nslices;
        nneighbors(3) = 2;
    end
    
    % Reshape and scale k from -pi to pi
    omega = 2*pi * fov./dim .* ...
        [reshape(ks(:,1,:,:),[],1), ...
        reshape(ks(:,2,:,:),[],1), ...
        reshape(ks(:,3,:,:),[],1)];
    
    % Create NUFFT object
    nufft_args = {dim,...
        nneighbors,...
        2*dim,...
        dim/2,...
        'table',...
        2^10,...
        'minmax:kb'};
    G = Gnufft(true(dim), [{omega}, nufft_args(:)']);
    
    % Create density compensation using pipe algorithm
    dcf = pipedcf(G,args.itrmax);
    
%% Reconstruct image using NUFFT
    % Reconstruct point spread function
    psf = reshape(G' * ones(numel(ks(:,1,:,:)),1), dim);
    
    % Initialize output array and progress string
    im = zeros([dim,info.ncoils,length(args.frames)]);
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
            im(:,:,:,coiln,framen) = reshape(G' * (dcf.*data), dim);
        end
        
    end
    
    % Combine coils using RMS
    im = sqrt( squeeze( mean( im.^2, 4) ) );
    
    % Save results that aren't returned to nifti:
    if nargout < 1
        % Save timeseries
        writenii('./timeseries.nii', abs(im), ...
            fov, info.tr, args.scaleoutput);
        fprintf('\nTimeseries saved to timeseries.nii');
        
        % Save point spread function
        writenii('./psf.nii', abs(psf), ...
            fov, info.tr, args.scaleoutput);
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