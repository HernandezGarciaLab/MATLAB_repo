function [im,imc] = recon3dflex(varargin)
% function [im,imc] = recon3dflex(varargin)
%
% Part of umasl project by Luis Hernandez-Garcia and David Frey
% @ University of Michigan 2022
%
% Description: Function to reconstruct images from *3dflex ASL sequence
%   using Model-based SENSE recon with Pipe & Menon Density compensation
%
%
% Notes:
%   - this function does not write images, use writenii to write the
%       outputs to file
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
%           nechos x ntrains x ncoils)
%       - if either raw or info is left empty, function will read data from
%           the Pfile instead (see 'pfile' or type 'help readpfile' for
%           more information)
%       - default is empty (reads raw data from file)
%   - 'info'
%       - pfile information structure
%       - structure containing fields: ndat, nechos, nframes, ntrains,
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
%   - 'niter'
%       - maximum number of iterations for model-based recon
%       - integer describing number of iterations
%       - if 0 is passed, conjugate phase recon will be used
%       - if a value less than 0 is passed, NUFFT recon will be used
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
%           sensitivity for each coil, or 'espirit' to estimate
%       - default is empty
%   - 'nramp'
%       - number of ramp points in spiral to delete from data
%       - integer describing number of points in ramp
%       - default is 10
%   - 'nnav'
%       - number of navigator points in spiral to use in phase detrend
%       - integer describing number of points in navigator
%       - default is 10
%   - 'ndel'
%       - number of points in gradient/sampling delay
%       - integer describing number of points
%       - default is 0
%   - 'pdorder'
%       - maximum order of polynomial fit for phase detrending
%       - integer describing order of polynomial
%       - if pdorder < 0, no phase detrending will occur
%       - default is -1
%
% Function output:
%   - im:
%       - output timeseries image (coil combined)
%       - complex array of image dimension
%   - imc:
%       - output timeseries image of each coil
%       - complex array of image dimension
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
    'girf',         '20230206UHP3T_rbw125.girf', ... % girf file to use
    'ndel',         0, ... % Gradient sample delay
    'despike',      {[]}, ... % linked frames to despike
    'nramp',        [], ... % Number of ramp points in spiral traj
    'pdorder',      -1, ... % Order of phase detrending poly fit
    't2',           [], ... % Option to perform t2 compensation
    'scaleoutput',  1, ... % Option to scale output to full dynamic range
    'dt',           4, ... % Sampling period (usec)
    'outputtag',    [] ... % Output filename tag
    );

% Start timer
t = tic;

%% Set up data, arguments and info

% Parse through variable inputs using matlab's built-in input parser
args = vararginparser(defaults,varargin{:});

%% Read in the data
%{
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
%}

%parpool(feature('numcores'),'IdleTimeout',90);

%% Load data (new)
% Load in raw data
[raw,phdr] = readpfile(args.pfile);

% data is ordered differently in asl3dflex and vsasl3dflex
%if strcmpi(phdr.image.psdname,'asl3dflex') % new sequence
if 1
    ndat = phdr.rdb.frame_size;
    nechoes = phdr.rdb.user2;
    ncoils = phdr.rdb.dab(2) - phdr.rdb.dab(1) + 1;
    ntrains = phdr.rdb.user1;
    nframes = phdr.rdb.user0;
    tr = phdr.image.tr*1e-3;
    dim = phdr.image.dim_X;
    fov = phdr.image.dfov/10;
else % old sequence
    ndat = phdr.rdb.frame_size;
    nechoes = phdr.rdb.nslices;
    ncoils = phdr.rdb.dab(2) - phdr.rdb.dab(1) + 1;
    ntrains = phdr.image.user0;
    nframes = phdr.rdb.nframes / ntrains;
    tr = phdr.image.tr*1e-3;
    dim = phdr.image.dim_X;
    fov = phdr.image.dfov/10;
end

%{
info.nleaves --> ntrains
info.nslices --> nechoes
%}

% reshape: ndat x ntrains*nframes x nechoes x 1 x ncoils
%           --> ndat x ntrains x nframes x nechoes x ncoils
raw = reshape(raw,ndat,ntrains,nframes,nechoes,ncoils);
% permute: ndat x ntrains x nframes x nechoes x ncoils
%           --> nframes x ndat x nechoes x ntrains x ncoils
raw = permute(raw,[3,1,4,2,5]);
dt = args.dt;

if strcmp(args.frames,'all')
    args.frames = 1:nframes;
else
    nframes = length(args.frames);
end

workDir = pwd;

%% Process Trajectory (new)
% Load in kspace trajectory & view transformation matrices
ktraj = dir('ktraj*.txt');
ktraj = load(ktraj(1).name);
kviews = dir('kviews*.txt');
kviews = load(kviews(1).name);

% Reshape transformation matrices as an array of 3x3 matrices
T = permute(reshape(kviews(:,end-8:end)',3,3,[]),[2,1,3]);

% Allocate space for entire trajectory
ks = zeros(ndat,3,nechoes,ntrains);

% Transform each view
for trainn = 1:ntrains
    for echon = 1:nechoes
        % Index the transformation matrix for current view
        mtxi = (trainn-1)*nechoes + echon;

        % Transform the trajectory
        ks(:,:,echon,trainn) = ktraj*T(:,:,mtxi)';
    end
end

% (old) Determine nramp and navpoints from 1st platter envelope
navpts = find(vecnorm(ktraj,2,2) < args.tol);
navpts(navpts > ndat*3/4) = []; % reject navpts in first 1/4
navpts(navpts < ndat*1/4) = []; % reject navpts in last 1/4

if isempty(args.nramp)
    args.nramp = 10;
    %[~,args.nramp] = max(vecnorm(ktraj(1:round(ndat/2),:),2,2));
    %args.nramp = round(0.1*args.nramp) ; % reject only up to 10% of the max
end

%%
%% Process Trajectory  (old section)
%{
    % Get un-transformed kspace trajectory and view transformations
    ks_0 = load([info.wd '/ktraj_cart.txt']);
    kviews = load([info.wd '/kviews.txt']);
    
    % Determine trajectory type based on kz encoding fractions
    isSOS = any(kviews(:,3) < 1);


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
%}

% Correct for gradient distortion if GIRF is passed
if ~isempty(args.girf)
    fprintf('\nCorrecting gradient distortion based on GIRF file: %s...', ...
        args.girf);
    girf = readmatrix(args.girf,'FileType','text');
    grad = padarray(diff(ks,1),[1,0,0,0],0,'pre');
    grad_corr = zeros(size(grad));
    for axis = 1:3
        tmp = convn(grad(:,axis,:,:),girf(:,axis),'full');
        grad_corr(:,axis,:,:) = tmp(1:ndat,1,:,:);
    end
    ks = cumsum(grad_corr,1);
    fprintf(' Done.')
end

% Clip echoes from end of echo train if specified
raw = raw(:, :, 1+args.clipechoes(1):nechoes-args.clipechoes(2), :, :);
ks = ks(:, :, 1+args.clipechoes(1):nechoes-args.clipechoes(2), :);
nechoes = nechoes - sum(args.clipechoes);

% Make timing array
ti = dt*(0:ndat-1)' + round((phdr.rdb.te - dt*ndat)/2);
ti = repmat(ti,1,ntrains,nechoes);
ti = ti + dt*ndat*permute(0:nechoes-1,[1 3 2]);

%% Apply corrections/filters
% Despike kspace
if ~isempty(args.despike)
    %{
	ttmp = tic;
        fprintf('\nDespiking kspace samples along timecourse...')
        [raw,nspikes] = despike1d(raw,1,'linked',args.despike);
        fprintf(' Done. %d/%d (%.2f%%) of data replaced. %.2fs elapsed.', ...
            nspikes, 2*numel(raw), 100*nspikes/(2*numel(raw)), toc(ttmp));
    %}
    fprintf('\nDespiking kspace...')
    sr = size(raw);


    Nspikes=0;
    Ndata = prod(size(raw));
    SDlevel = 2;  % Threshold for despiking : N of standard deviations

    for n=1:sr(3) % shots
        for m=1:sr(4) % N echoes per train (slices)
            for p=1:sr(5) % coils
                % odd time frames
                tmp = squeeze(raw(args.despike{1},:,n,m,p));
                [raw(args.despike{1},:,n,m,p) spikes] = ...
                    smoothOutliers(tmp,SDlevel);
                Nspikes = Nspikes + spikes;
                % even time frames
                tmp = squeeze(raw(args.despike{2},:,n,m,p));
                [raw(args.despike{2},:,n,m,p) spikes] = ...
                    smoothOutliers(tmp, SDlevel);
                Nspikes = Nspikes + spikes;

            end
        end
    end
    fprintf('\rFound %d spikes in %d data.  (%f percent)', ...
        Nspikes, Ndata, Nspikes/Ndata*100);
    fprintf('\n...Despiking Done')
end

% Perform phase detrending
if args.pdorder > -1
    fprintf('\nPerforming %d th order detrending in each echo...', args.pdorder);
    raw = phasedetrend(raw,navpts,args.pdorder);
end

if strcmpi(args.t2, 'fit')
    % Make t2map in kspace
    fprintf('\nFitting average T2 to raw data... ')
    args.t2 = fitt2(raw,navpts,phdr.rdb.te,dt);
    args.t2 = args.t2 * dt * 1e-3;
    fprintf(' Done. Fitted T2 value: %.1fms', args.t2);
end

% Correct for gradient sample delay
if abs(args.ndel) > 0
    fprintf('\nCorrecting for delay : %d samples', args.ndel)
    raw = circshift(raw,args.ndel,2);
end

% Remove ramp points : Not as many as before!!
fprintf('\nRemoving %d ramp points from data...\n', args.nramp);
ks([1:args.nramp ndat-args.nramp:ndat],:,:,:) = [];
raw(:,[1:args.nramp ndat-args.nramp:ndat],:,:,:) = [];
ti([1:args.nramp ndat-args.nramp:ndat],:,:) = [];
ndat = ndat - 2*args.nramp;

%% Set up reconstruction
% Vectorize fov/dim and apply interpolation factors
fov = fov*ones(1,3) / args.zoomfactor;
dim = round(dim*ones(1,3) * args.resfactor);
%if isSOS && ~args.isovox
%    fov(3) = nechoes*info.slthick / args.zoomfactor;
%    dim(3) = nechoes * args.resfactor;
%end

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
niter = args.niter; % extract to avoid broadcasting args into parfor

% Create density compensation using pipe algorithm
fprintf('\nCalculating density compensation (Pipe-Menon method, 15 itr)...');
dcf = pipedcf(Gm.Gnufft,15);
W = Gdiag(dcf(:)./Gm.arg.basis.transform);
fprintf(' Done.')

% Reconstruct point spread function
fprintf('\nCalculating PSF ...');
psf = Gm' * reshape(W * ones(numel(ks(:,1,:,:)),1), [], 1);
fprintf(' Done.');

% Build t2 relaxation map into Gmri object
if ~isempty(args.t2) && args.niter >= 0
    fprintf('\nCreating R2 Gmri zmap for reconstruction model...');
    R2 = 1/(args.t2*1e3); % convert to 1/usec
    Gm = feval(Gm.arg.new_zmap,Gm,ti(:),R2*ones(dim),1);
    fprintf(' Done.');
end

% Correct smap for single-coil data
if isempty(args.smap) && ncoils == 1
    args.smap = ones(dim);
end

% Make quadratic regularizer for PCG recon
R = Reg1(ones(dim), 'beta', 2^-12 * numel(ks(:,1,:,:)));
C = block_fatrix({R.C}, 'type', 'col');

% Use only nufft for fourier encoding if niter < 0
if args.niter < 0
    fprintf('\nSetting up NUFFT reconstruction (since niter < 0) for')
    Gm = Gm.Gnufft;
elseif args.niter == 0
    fprintf('\nSetting up CP reconstruction model for')
else
    fprintf('\nSetting up %d-iteration PCG reconstruction model for',args.niter)
end

if ~isempty(args.smap) || ncoils == 1
    fprintf(' image timeseries (with sensitivity encoding)...');

    % Correct size of W (weights for density compensation)
    W = Gdiag(repmat(dcf(:),1,ncoils));

    % Incorporate sensitivity encoding into system matrix
    Ac = repmat({[]},ncoils,1);
    for coiln = 1:ncoils
        tmp = args.smap(:,:,:,coiln);
        tmp = Gdiag(tmp);
        Ac{coiln} = Gm * tmp;
    end
    A = block_fatrix(Ac, 'type', 'col');

    % Reshape data for recon (ndata x ncoils x nframes)
    tmp = permute(raw(:,:,:,:),[2:5,1]);
    data_all = reshape(tmp,[],ncoils,nframes);
    clear tmp;
        

    % Initialize output image (dim x nframes)
    im = zeros([dim,nframes]);

else
    fprintf(' coil-wise images for frame 1 (no sensitivity encoding)...');

    % System matrix is only fourier encoding
    A = Gm;

    % Reshape data for recon (ndata x ncoils x nframes)
    data_all = reshape(raw(1,:,:,:,:),[],1,ncoils);

    % Initialize output image (dim x ncoils)
    im = zeros([dim,ncoils]);


end
fprintf(' Done.')

%% Reconstruct images (OLD CODE)
%{
    fprintf('\nReconstructing...');

    % Loop through image series (frames or coils index)
    fprintf('\nExecuting the Reconstruction. %d Iteration(s)...', niter);
    for n = 1:size(im,4)
        
        % Get indexed data
        data = data_all(:,:,n);
        
        % Recon preconditioner using conjugate-phase
        im_cp = A' * reshape(W * data, [], 1);
        im_cp = embed(im_cp,true(dim*ones(1,3)));
        im_cp = ir_wls_init_scale(A, data(:), im_cp);
        
        % Recon using preconditioned conjugate gradient (iterative)
        if niter > 0
            im_pcg = qpwls_pcg1(im_cp(true(dim*ones(1,3))), A, 1, data(:), C, ...
                'niter', niter);
            im(:,:,:,framen) = embed(im_pcg,true(dim*ones(1,3)));
            
        % ...or save image with CP recon
        else
            im(:,:,:,n) = reshape(im_cp, dim);
        end
        
    end
    fprintf(' Done.');
%}


%% Reconstruct images (NEW CODE)

% Initialize output image (dim x nframes)
im = zeros([dim,nframes]);
if ncoils==1
    args.smap=[];
end

if isempty(args.smap) % Coil-wise recon with RMS coil combo

    for framen = 1:length(args.frames)

        if ncoils==1  % single coil- keep the phase info
            data = reshape(raw(args.frames(framen),:,:,:,1),[],1);
            im(:,:,:,framen) = reshape(Gm' * (W*data(:)),dim);

        else  % multi-coil commbine with RMS
            data = zeros(dim);

            parfor coiln = 1:ncoils
                % Recon data
                data = reshape(raw(args.frames(framen),:,:,:,coiln),[],1);
                imc(:,:,:,coiln,framen) = reshape(Gm' * (W*data(:)),dim);

                % rescale to conserve energy
                tmp = raw(args.frames(framen), :,:,:, coiln);
                nrg_in = norm(tmp(:));
                tmp = imc(:,:,:,coiln,framen);
                nrg_out = norm(tmp(:));
                SCL = nrg_in / nrg_out;

                imc(:,:,:,coiln, framen) = imc(:,:,:,coiln,framen) * SCL;

            end
            im(:,:,:,framen) = sqrt(mean(imc(:,:,:,:,framen).^2,4));

        end

    end

else % CP-SENSE recon

    % Incorporate sensitivity encoding into system matrix
    Ac = repmat({[]},ncoils,1);
    for coiln = 1:ncoils
        tmp = args.smap(:,:,:,coiln);
        tmp = Gdiag(tmp(true(dim.*ones(1,3))),'mask',true(dim.*ones(1,3)));
        Ac{coiln} = Gm * tmp;
    end
    A = block_fatrix(Ac, 'type', 'col');

    % Reshape density weighting matrix
    W = Gdiag(repmat(dcf(:)./Gm.arg.basis.transform,1,ncoils));

    % Recon data
    parfor framen = args.frames
        data = reshape(raw(args.frames(framen),:,:,:,:),[],1);
        im(:,:,:,framen) = embed(A' * reshape(W * data, [], 1), ...
            true(dim.*ones(1,3)));
    end


end



%% Save and return output data

if nargout < 1
    % Save timeseries

    if (ncoils>1) && isempty(args.smap)
        % Save coil images
        writenii([workDir,'/coils_mag',args.outputtag], abs(imc), ...
            'fov', fov, 'doScl', args.scaleoutput);
        writenii([workDir,'/coils_ang',args.outputtag], angle(imc), ...
            'fov', fov, 'doScl', args.scaleoutput);
        fprintf('\nCoil images (frame 1) saved to coil_*%s.nii',args.outputtag);
    end

    writenii([workDir,'/timeseries_mag',args.outputtag], abs(im), ...
        'fov', fov, 'tr', tr*ntrains, 'doscl', args.scaleoutput);
    writenii([workDir,'/timeseries_ang',args.outputtag], angle(im), ...
        'fov', fov, 'tr', tr*ntrains, 'doscl', args.scaleoutput);
    fprintf('\nTimeseries saved to timeseries_mag%s.nii',args.outputtag);

    % Save point spread function
    writenii([workDir,'/psf',args.outputtag], abs(psf), ...
        'fov', fov, 'tr', 1, 'doscl', args.scaleoutput);
    fprintf('\nPoint spread function saved to psf%s.nii',args.outputtag);

    % Destroy im variable so it will not be returned
    clear im
else
    % Notify user that images will not be saved
    fprintf('\nImages will not be saved to file since output is returned');
end

% Save and print elapsed time
t = toc(t);
fprintf('\nRecon complete. Total elapsed time: %.2fs\n',t);

end

%% readpfile function definition (old - superceded)
%{
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
%}

%% phasedetrend function definition
function raw_corr = phasedetrend(raw,navpts,pdorder)
% In the ideal case the phase should be flat and constant from
% echo to echo in the navigator.
% the idea is to fit the signal in the navigator data segment to a
% polynomial trend. We assume that this trend is present in the whole echo and
% that it is undesirable, so we subtract this fitted polynomial from the
% echo

% Get dimensions
nframes = size(raw,1);
ndat = size(raw,2);
nechos = size(raw,3);
ntrains = size(raw,4);
ncoils = size(raw,5);

% The linear model for the noise is asum of polynomials:
% y = x0 + x + x^2 + ....x^pdorder
A = (navpts(:) - round(ndat/2)).^(pdorder:-1:0);

% Intitialize fits array
fits = zeros(nframes,ndat,nechos,ntrains,ncoils);

% Loop through all echos
for framen = 1:nframes
    for echon = 1:nechos
        for trainn = 1:ntrains
            for coiln = 1:ncoils
                % Fit polynomial to phase of indexed echo using
                % simple linear least squares.  The fit is only in the
                % center region : the navigator points.
                echo = raw(framen,:,echon,trainn,coiln);
                y = angle(echo(navpts));
                betas = pinv(A)*y(:);

                % now remove the fitted polynomial from the WHOLE echo,
                % noting that the center of the echo is the same data
                % as the center of the navigators
                fits(framen,:,echon,trainn,coiln) = ...
                    ((1:ndat) - round(ndat/2))'.^(pdorder:-1:0) * betas;
            end
        end
    end
end

% Correct echo by subtracting out fits from phase
raw_corr = raw.*exp(-1i*fits);

end
%% smoothOutliers  function definition
function [out Nspikes]= smoothOutliers(raw, threshold)
% function [out Nspikes_found] = smoothOutliers(raw, threshold);
%
% Given a 2D data set, smooth out the outliers along the row direction
% This means replacing outliers with the average of their neighbors
% along the row dimension
%
% outliers are defined as being threshold * stddev(row) away from the mean

% detrend the data, but preserve its mean value
mu = mean(raw, 1);
raw = detrend(raw,2);
raw = raw + mu; % preserve the mean value

% allocate space for output:
out = raw;
mu = mean(raw(2:end-1, :), 1);
sigma = std(raw(2:end-1, :), 1);
Ndata = prod(size(raw));
Nspikes = 0;
% go through all the columns
for col=1:size(out,2)
    % identify spikes at each k-space location:
    % first we detrend the data, then we identify the spikes
    % datacolumn = detrend(raw(:,col));
    datacolumn = (raw(:,col));
    outlier_inds = find( abs(datacolumn-mu(col)) > sigma(col)*threshold );
    Nspikes = Nspikes + length(outlier_inds);

    % now replace the spikes by the mean of their neighbors
    for r=1:length(outlier_inds)
        % Now replace spikes by a mean of its neighbors (in time) -
        % (linear interpolation)
        % make sure we're not overwriting the field map (the first row)!
        if (outlier_inds(r)>1 &&   outlier_inds(r) < length(datacolumn) -1)

            row = outlier_inds(r);
            out(row, col) = 0.5*(raw(row-1, col) + raw(row+1, col)) ;

        end
    end
end

% fprintf('\rFound %d spikes in %d data.  (%f percent)', ...
%     Nspikes, Ndata, Nspikes/Ndata*100);
end


%% t2comp function definition
function t2 = fitt2(raw,navpts,te,dt)

% Get dimensions
ndat = size(raw,2);
ndat_all = round(te/dt);
nechos = size(raw,3);
ntrains = size(raw,4);
ncoils = size(raw,5);

% Define design matrix for lsq poly fit
x = round((ndat_all - ndat)/2) + (navpts + ndat_all * (0:ntrains-1));
A = x(:).^[1 0];

% Initialize R2
t2 = zeros(nechos,ncoils);

% Loop through all echos
for echon = 1:nechos
    for coiln = 1:ncoils

        % Get current echo train with full echo time
        echo = padarray(raw(1,:,echon,:,coiln), ...
            [0 floor((ndat_all - ndat)/2) 0 0 0],'pre');
        echo = padarray(echo, ...
            [0 ceil((ndat_all - ndat)/2) 0 0 0],'post');
        echo = abs(reshape(echo,ndat_all*ntrains,1));

        % Fit the decay using least squares
        y = echo(x(:));
        b = pinv(A) * log(y);

        % Save R2 values
        t2(echon,coiln) = 1/abs(b(1));

    end
end

% Get average R2 value
t2 = mean(t2(~isinf(t2(:))));

end
