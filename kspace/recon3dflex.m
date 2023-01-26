function im = recon3dflex(varargin)
% function im = recon3dflex(varargin)
%
% Part of umasl project by Luis Hernandez-Garcia and David Frey
% @ University of Michigan 2022
%
% Description: Function to reconstruct images from *3dflex ASL
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
%   - 'niter'
%       - maximum number of iterations for model-based recon
%       - integer describing number of iterations
%       - if 0 is passed, conjugate phase recon will be used
%       - default is 0
%   - 'frames'
%       - frames in timeseries to recon
%       - integer array containing specific frames to recon
%       - if 'all' is passed, all frames will be reconned
%       - default is 'all'
%   - 'clipechoes'
%       - number of echoes (slices) to remove from echo train
%       - 1 x 2 integer array describing number of echoes to clip from
%           start and beginning
%       - default is [0 0]
%   - 'resfactor'
%       - resolution (dimension) upsampling factor
%       - float/double describing factor
%       - passing a value > 1 will result in higher output image dimension
%       - default is 1
%   - 'zoomfactor'
%       - field of view zoom factor
%       - float/double describing factor
%       - passing a value > 1 will result in a smaller output image fov
%       - default is 1
%   - 'smap'
%       - coil sensitivity map for multi-coil datasets
%       - complex double/float array of dim x ncoils representing
%           sensitivity for each coil, or 'estimate'
%       - if left empty, recon will write 'coils_*.nii' or return it as
%           output so user can make a sense map
%       - default is empty
%   - 'isovox'
%       - option to use isotropic voxel sizes
%       - boolean integer (0 or 1) describing whether or not voxels are
%           isotropic
%       - some old data has a different z fov/dim, which would require
%           isovox = 0 to properly reconstruct, ideal for vsasl3dflex data
%           acquired before august 2022
%       - default is 1
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
%   - 't2'
%       - T2 value for T2 decay compensation
%       - float/double describing T2 in ms
%       - if 'fit' is passed, recon will try to estimate an average t2
%       - if left empty, no t2 will not be considered in signal model
%       - default is empty
%   - 'scaleoutput'
%       - option to scale nii files to full dynamic range
%       - boolean integer (0 or 1) to use or not
%       - type 'help writenii' for more information
%       - default is 1
%   - 'outputtag'
%       - tag to append on output file names
%       - string containing tag
%       - default is empty (no tag)
%
% Function output:
%   - im:
%       - output timeseries image, or coil-wise image if smap is not passed
%       - complex array of image dimension
%       - if im is not returned, all images will be saved to nii files
%       - if im is returned, no images will be saved to nii files
%

    % Define default arguments
    defaults = struct(...
        'raw',          [], ... % Raw data
        'info',         [], ... % Info structure
        'pfile',        'P*.7', ... % Search string for Pfile
        'tol',          1e-5, ... % Kspace distance tolerance
        'niter',        0, ... % Max number of iterations for IR
        'frames',       'all', ... % Frames to recon
        'clipechoes',   [0 0], ... % Number of echoes to remove
        'resfactor',    1, ... % Resolution upsampling factor
        'zoomfactor',   1, ... % FOV zoom factor
        'smap',         [], ... % Sensitivity map for coil combination
        'isovox',       1, ... % flag for isotropic voxels between x&y / z
        'ndel',         0, ... % Gradient sample delay
        'nramp',        [], ... % Number of ramp points in spiral traj
        'pdorder',      -1, ... % Order of phase detrending poly fit
        't2',           [], ... % Option to perform t2 compensation
        'scaleoutput',  1, ... % Option to scale output to full dynamic range
        'outputtag',    [] ... % Output filename tag
        );

    % Start timer
    t = tic;
    
%% Set up recon
    % Parse through variable inputs using matlab's built-in input parser
    args = vararginparser(defaults,varargin{:});
    
    % Get raw data and info
    if (isempty(args.raw) || isempty(args.info)) % If reading from Pfile
        [raw,info] = formatpfile(args.pfile);
    else % If raw/scaninfo is user specified
        raw = args.raw;
        info = args.info;
    end
    
    % Append output tag
    if ~isempty(args.outputtag)
        args.outputtag = ['_',args.outputtag];
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
    ks_0 = load([info.wd '/ktraj_cart.txt']);
    kviews = load([info.wd '/kviews.txt']);
    
    % Determine trajectory type based on kz encoding fractions
    isSOS = any(kviews(:,3) < 1);
    
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
    
    % Clip echoes from end of echo train if specified
    raw = raw(:,:,:,1+args.clipechoes(1):info.nslices-args.clipechoes(2),:);
    ks = ks(:,:,:,1+args.clipechoes(1):info.nslices-args.clipechoes(2),:);
    info.nslices = info.nslices - sum(args.clipechoes);
    
    % Make timing array
    ti = info.dt*(0:info.ndat-1)' + round((info.te - info.dt*info.ndat)/2);
    ti = repmat(ti,1,info.nleaves,info.nslices);
    ti = ti + info.dt*info.ndat*permute(0:info.nslices-1,[1 3 2]);
        
%% Apply corrections/filters
    % Perform phase detrending
    if args.pdorder > -1
        raw = phasedetrend(raw,navpts,args.pdorder);
    end
    
    if strcmpi(args.t2, 'fit')
        % Make t2map in kspace
        args.t2 = fitt2(raw,navpts,info.te,info.dt);
        args.t2 = args.t2 * info.dt * 1e-3;
        fprintf('\nFitted T2 value: %.1fms', args.t2);
    end
    
    % Correct for gradient sample delay
    if abs(args.ndel) > 0
        raw = circshift(raw,args.ndel,2);    
    end
    
    % Remove ramp points
    fprintf('\nRemoving %d ramp points from data...', args.nramp);
    ks([1:args.nramp info.ndat-args.nramp:info.ndat],:,:,:) = [];
    raw(:,[1:args.nramp info.ndat-args.nramp:info.ndat],:,:,:) = [];
    ti([1:args.nramp info.ndat-args.nramp:info.ndat],:,:) = [];
    info.ndat = info.ndat - 2*args.nramp;
    
%% Set up NUFFT
    % Vectorize fov/dim and apply interpolation factors
    fov = info.fov*ones(1,3) / args.zoomfactor;
    dim = round(info.dim*ones(1,3) * args.resfactor);
    if isSOS && ~args.isovox
        fov(3) = info.nslices*info.slthick / args.zoomfactor;
        dim(3) = info.nslices * args.resfactor;
    end
    
    % Create Gmri object
    kspace = [reshape(ks(:,1,:,:),[],1), ...
        reshape(ks(:,2,:,:),[],1), ...
        reshape(ks(:,3,:,:),[],1)];
    nufft_args = {dim,...
        6*ones(1,3),...
        2*dim,...
        dim/2,...
        'table',...
        2^10,...
        'minmax:kb'};
    Gm = Gmri(kspace, true(dim), ...
        'fov', fov, 'basis', {'rect'}, 'nufft', nufft_args(:)');
    
    % Create density compensation using pipe algorithm
    fprintf('\nCreating density compensation...');
    dcf = pipedcf(Gm.Gnufft,15);
    W = Gdiag(dcf(:)./Gm.arg.basis.transform);
    
    % Build t2 relaxation map into Gmri object
    if ~isempty(args.t2)
        fprintf('\n');
        R2 = 1/(args.t2*1e3); % convert to 1/usec
        Gm = feval(Gm.arg.new_zmap,Gm,ti(:),R2,1);
    end
    
%% Reconstruct image using NUFFT
    % Reconstruct point spread function
    psf = Gm' * reshape(W * ones(numel(ks(:,1,:,:)),1), [], 1);
    
    % If smap is empty, recon frame 1 coil images for smap construction
    if isempty(args.smap)
        
        % Initialize image
        im = zeros([dim, info.ncoils]);
        
        % Loop through coils
        for coiln = 1:info.ncoils
            
            % Get data for first frame and current coil
            data = reshape(raw(1,:,:,:,coiln),[],1);
            
            % Recon using cp without sensitivity encoding
            im_cp = Gm' * W * data;

            % Recon using pcg (iterative) without sensitivity encoding
            if args.niter > 0
                im_pcg = ir_wls_init_scale(Gm, data, im_cp);
                im_pcg = qpwls_pcg1(im_pcg, Gm, 1, data_k, C, ...
                    'niter', args.niter);
                im(:,:,:,coiln) = reshape(im_pcg, dim);
            
            % ...or save image with CP recon
            else
                im(:,:,:,coiln) = reshape(im_cp, dim);
            end
        end
        
        % Save to file if image is not outputted
        if nargout < 1
            writenii([info.wd,'/coils_mag_',args.outputtag], abs(im), ...
                'fov', fov, 'doScl', args.scaleoutput);
            writenii([info.wd,'/coils_ang_',args.outputtag], angle(im), ...
                'fov', fov, 'doScl', args.scaleoutput);
            fprintf('\nCoil images (frame 1) saved to coil_*%s.nii',args.outputtag);
        else
            fprintf('\nImages will not be saved to file since coils image is returned');
        end
        
    % ...or recon with sensitivity encoding
    else
        
        % Correct size of W
        W = Gdiag(repmat(dcf(:),1,info.ncoils));
        
        % Incorporate sensitivity encoding into system matrix
        Ac = repmat({[]},info.ncoils,1);
        for coiln = 1:info.ncoils
            tmp = args.smap(:,:,:,coiln);
            tmp = Gdiag(tmp);
            Ac{coiln} = Gm * tmp;
        end
        A = block_fatrix(Ac, 'type', 'col');
       
        % Save point spread function to file if specified
        if nargin < 1
            writenii([info.wd,'/psf',args.outputtag], abs(psf), ...
                'fov', fov, 'tr', 1, 'doscl', args.scaleoutput);
            fprintf('\nPoint spread function saved to psf%s.nii',args.outputtag);
        end
        
        % Initialize image
        im = zeros([dim, length(args.frames)]);
        fprintf('\nReconning images...');
        
        % make quadratic regularizer?
        R = Reg1(ones(dim), 'beta', 2^-12 * numel(ks(:,1,:,:)));
        C = block_fatrix({R.C}, 'type', 'col');
        
        % Loop through frames
        for framen = args.frames
            
            % Get data for current frame
            data = reshape(raw(framen,:,:,:,:),[],info.ncoils);
                      
            % Recon using cp
            im_cp = A' * reshape(W * data, [], 1);
            
            % Recon using pcg (iterative) without sensitivity encoding
            if args.niter > 0
                im_pcg = ir_wls_init_scale(A, data(:), im_cp);
                im_pcg = qpwls_pcg1(im_pcg, A, 1, data(:), C, ...
                    'niter', args.niter);
                im(:,:,:,savef) = reshape(im_pcg, dim);
            
            % ...or save image with CP recon
            else
                im(:,:,:,savef) = reshape(im_cp, dim);
            end
            
            % Update saved frame index
            savef = savef + 1;
        end
        
        % Save results that aren't returned to nifti:
        if nargout < 1
            % Save timeseries
            writenii([info.wd,'/timeseries_mag',args.outputtag], abs(im), ...
                'fov', fov, 'tr', info.tr, 'doscl', args.scaleoutput);
            writenii([info.wd,'/timeseries_ang',args.outputtag], angle(im), ...
                'fov', fov, 'tr', info.tr, 'doscl', args.scaleoutput);
            fprintf('\nTimeseries saved to timeseries_mag%s.nii',args.outputtag);
            
            % Clear im so it won't be returned
            clear im;
            
        % or remind user that output will not be returned...
        else
            fprintf('\nImages will not be saved to file since timeseries is returned');
        end
        
    end
    
    % Save and print elapsed time
    t = toc(t);
    fprintf('\nRecon complete. Total elapsed time: %.2fs\n',t);
    
end

%% readpfile function definition
function [raw,info] = formatpfile(searchstr)

    [raw,h] = readpfile(searchstr);
    d = dir(searchstr);
    
    % Create info struct based on header info
    info = struct(...
        'ndat',     h.rdb.da_xres, ... % Number of points / echo
        'nleaves',  h.image.user0, ... % Number of interleaves
        'nframes',  h.image.user1, ... % Number of temporal frames
        'nslices',  h.image.slquant, ... % Number of slices
        'ncoils',   h.rdb.dab(2) - h.rdb.dab(1) + 1, ... % Number of coils
        'tr',       h.image.tr*1e-3, ... % TR (ms)
        'te',       h.image.te, ... % TE (usec)
        'dt',       4, ... % Sampling period (usec)
        'dim',      h.image.dim_X, ... % Image x/y dimension
        'fov',      h.image.dfov/10, ... % FOV (cm)
        'slthick',  h.image.slthick/10, ... % Slice Thickness (cm)
        'wd',       d(1).folder ... % Working directory
        );
    
    % Reshape data
    raw = reshape(raw, ...
        info.ndat,info.nframes*info.nleaves+1,info.nslices,info.ncoils);
    raw(:,1,:,:) = []; % cut out baseline
    raw = reshape(raw, ...
        info.ndat,info.nleaves,info.nframes,info.nslices,info.ncoils);
    raw = permute(raw, [3,1,2,4,5]);

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
                    y = angle(echo(navpts));
                    betas = pinv(A)*y(:);
                    fits(framen,:,leafn,slicen,coiln) = ...
                        ((1:ndat) - round(ndat/2))'.^(pdorder:-1:0) * betas;
                end
            end
        end
    end
    
    % Correct echo by subtracting out fits from phase
    raw_corr = raw.*exp(-1i*fits);

end

%% t2comp function definition
function t2 = fitt2(raw,navpts,te,dt)

    % Get dimensions
    ndat = size(raw,2);
    ndat_all = round(te/dt);
    nleaves = size(raw,3);
    nslices = size(raw,4);
    ncoils = size(raw,5);
    
    % Define design matrix for lsq poly fit
    x = round((ndat_all - ndat)/2) + (navpts + ndat_all * (0:nslices-1));
    A = x(:).^[1 0];
    
    % Initialize R2
    t2 = zeros(nleaves,ncoils);
    
    % Loop through all echos
    for leafn = 1:nleaves
        for coiln = 1:ncoils

            % Get current echo train with full echo time
            echo = padarray(raw(1,:,leafn,:,coiln), ...
                [0 floor((ndat_all - ndat)/2) 0 0 0],'pre');
            echo = padarray(echo, ...
                [0 ceil((ndat_all - ndat)/2) 0 0 0],'post');
            echo = abs(reshape(echo,ndat_all*nslices,1));

            % Fit the decay using least squares
            y = echo(x(:));
            b = pinv(A) * log(y);

            % Save R2 values
            t2(leafn,coiln) = 1/abs(b(1));

        end
    end
    
    % Get average R2 value
    t2 = mean(t2(~isinf(t2(:))));

end