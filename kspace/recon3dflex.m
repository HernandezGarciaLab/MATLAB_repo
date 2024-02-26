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
%   - 'kdespike'
%       - frame groups to despike
%       - cell array containing vectors of frame indicies to perform
%           despiking on
%       - for example, {1:2:100, 2:2:100} would despike even and odd frames
%           seperately
%       - default is an empty cell (no despiking)
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
%   - 'ccfac'
%       - coil compression factor
%       - float from 0 to 1 describing factor of reduction in # of coils
%       - default is 1
%
% Function output:
%   - im:
%       - output timeseries image (coil combined)
%       - complex array of image dimension
%   - imc:
%       - output timeseries image of each coil
%       - complex array of image dimension
%

% Start parallel pool
if isempty(gcp('nocreate'))
    parpool(feature('numcores'),'IdleTimeout',90);
end

% Define default arguments
defaults = struct( ...
    'pfile',        './P*.7', ... % Search string for Pfile
    'smap',         [], ... % Sensitivity maps
    'niter',        0, ... % Number of iterations for pcg recon
    'frames',       'all', ... % Frames to recon
    'ndel',         0, ... % Number of points in gradient/sampling delay
    'nramp',        10, ... % Number of ramp points in trajectory to remove
    'nnav',         100, ... % Number of navigator points to use in PD
    'clipechoes',   [0,0], ... % Number of echoes to remove from each end
    'pdorder',      -1, ... % Order of polynomial to fit for phase detrend
    'fdespike',     {[]}, ... % linked frames to despike
    'zoomfac',      1, ... % FOV downscaling factor
    'resfac',       1, ... % Resolution upscaling factor
    'ccfac',        1 ... % Coil compression factor
    );


% Parse through variable inputs using matlab's built-in input parser
args = vararginparser(defaults,varargin{:});

% Get target pfile
pfile = dir(args.pfile);
if isempty(pfile)
    error('no pfiles found with search string: %s', args.pfile);
end
pdir = pfile(1).folder;
pfile = pfile(1).name;

% Load in kspace trajectory & view transformation matrices
ktraj = dir([pdir,'/ktraj*.txt']);
ktraj = load([ktraj(1).folder,'/',ktraj(1).name]);
kviews = dir([pdir,'/kviews*.txt']);
kviews = load([kviews(1).folder,'/',kviews(1).name]);
fprintf('successfully loaded trajectory files\n');

% Reshape transformation matrices as an array of 3x3 matrices
T = permute(reshape(kviews(:,end-8:end)',3,3,[]),[2,1,3]);

% Load in raw data
[raw,phdr] = readpfile([pdir,'/',pfile]);
fprintf('successfully loaded pfile data\n');
ndat = phdr.rdb.frame_size;
nechoes = phdr.rdb.user2;
ncoils = phdr.rdb.dab(2) - phdr.rdb.dab(1) + 1;
ntrains = phdr.rdb.user1;
nframes = phdr.rdb.user0;
dim = phdr.image.dim_X*args.resfac;
fov = phdr.image.dfov/10/args.zoomfac;
tr = phdr.image.tr*1e-3;

% reshape: ndat x ntrains*nframes x nechoes x 1 x ncoils
%           --> ndat x ntrains x nframes x nechoes x ncoils
raw = reshape(raw,ndat,ntrains,[],nechoes,ncoils);
% permute: ndat x ntrains x nframes x nechoes x ncoils
%           --> nframes x ndat x nechoes x ntrains x ncoils
raw = permute(raw,[3,1,4,2,5]);

% Compress coils using mirt's SVD coil compression
if args.ccfac < 1
    nvcoils = round(args.ccfac*ncoils);
    fprintf('compressing %d coils --> %d virtual coils...\n', ncoils, nvcoils);
    raw = ir_mri_coil_compress(raw, 'ncoil', nvcoils);
    ncoils = nvcoils;
end

% Allocate space for entire trajectory
ktraj_all = zeros(ndat,3,nechoes,ntrains);

% Transform each view
fprintf('transforming trajectory...\n');
for trainn = 1:ntrains
    for echon = 1:nechoes
        % Index the transformation matrix for current view
        mtxi = (trainn-1)*nechoes + echon;
        
        % Transform the trajectory
        ktraj_all(:,:,echon,trainn) = ktraj*T(:,:,mtxi)';
    end
end

% Apply gradient/sample delay
if abs(args.ndel) > 0
    fprintf('applying %d sample gradient delay...\n', args.ndel);
    raw = circshift(raw,[0,args.ndel,0,0,0]);
end

% Despike the data in kspace
if ~isempty(args.fdespike)
    fprintf('despiking data...\n');
    raw = kdespike(raw, args.fdespike);
end

% Remove unwanted frames
if strcmpi(args.frames,'all')
    args.frames = 1:nframes;
end
raw = raw(args.frames,:,:,:,:);
nframes = length(args.frames);

% Remove ramp points
if args.nramp > 0
    fprintf('removing %d ramp points...\n', args.nramp);
    raw = raw(:,args.nramp+1:end-args.nramp,:,:,:);
    ktraj_all = ktraj_all(args.nramp+1:end-args.nramp,:,:,:);
end

% Remove unwanted echoes
if any(args.clipechoes > 0)
    fprintf('removing %d echoes from beginning of train\n', args.clipechoes(1));
    fprintf('and %d echoes from end of train...\n', args.clipechoes(2));
    raw = raw(:,:,args.clipechoes(1)+1:end-args.clipechoes(2),:,:);
    ktraj_all = ktraj_all(:,:,args.clipechoes(1)+1:end-args.clipechoes(2),:);
end

% Perform phase detrending
if args.pdorder >= 0
    fprintf('performing %s order phase detrending...\n', iptnum2ordinal(args.pdorder));
    raw = phasedetrend(raw,args.nnav,args.pdorder);
end

% Format kspace into column vectors
kspace = [reshape(ktraj_all(:,1,:,:),[],1), ...
    reshape(ktraj_all(:,2,:,:),[],1), ...
    reshape(ktraj_all(:,3,:,:),[],1)];

% Create NUFFT object
fprintf('creating Gmri object...\n');
nufft_args = {dim*ones(1,3),...
    6*ones(1,3),...
    2*dim*ones(1,3),...
    dim*ones(1,3)/2,...
    'table',...
    2^10,...
    'minmax:kb'};
Gm = Gmri(kspace, true(dim*ones(1,3)), ...
    'fov', fov, 'basis', {'rect'}, 'nufft', nufft_args(:)');

% Calculate density weighting matrix
fprintf('calculating density weighting matrix...\n');
dcf = pipedcf(Gm.Gnufft);
W = Gdiag(dcf(:)./Gm.arg.basis.transform);

% Initialize image matrices
im = zeros([dim*ones(1,3),nframes]);
imc = zeros([dim*ones(1,3),ncoils,nframes]);

% Set smap to uniform if single coil
if ncoils == 1
    args.smap = ones(dim*ones(1,3));
end

% Coil-wise recon with RMS coil combo
if isempty(args.smap) || nargout == 2
    fprintf('performing coil-wise nufft recon...\n');
    
    for framen = 1:nframes
        parfor coiln = 1:ncoils
            % Recon data
            data = reshape(raw(framen,:,:,:,coiln),[],1);
            imc(:,:,:,coiln,framen) = reshape(Gm' * (W*data(:)),dim,dim,dim);
        end
        
        if ncoils == 1 % Save single coil data
            im(:,:,:,framen) = imc(:,:,:,1,framen);
        else % RMS coil combo
            im(:,:,:,framen) = sqrt(mean(imc(:,:,:,:,framen).^2,4));
        end
    end
end

% PCG/CP-SENSE recon
if ~isempty(args.smap)
    fprintf('performing PCG/CP-SENSE recon (%d iterations)...\n', args.niter);
    
    % Make regularizer
    R = Reg1(ones(dim*ones(1,3)), 'beta', 2^-12 * numel(raw(1,:,:,:,:))/3, ...
        'mask', true(dim*ones(1,3)));
    C = R.C;
    
    % Incorporate sensitivity encoding into system matrix
    Ac = repmat({[]},ncoils,1);
    for coiln = 1:ncoils
        tmp = args.smap(:,:,:,coiln);
        tmp = Gdiag(tmp(true(dim*ones(1,3))),'mask',true(dim*ones(1,3)));
        Ac{coiln} = Gm * tmp;
    end
    A = block_fatrix(Ac, 'type', 'col');
    
    % Reshape density weighting matrix
    W = Gdiag(repmat(dcf(:)./Gm.arg.basis.transform,1,ncoils));
    
    % Extract variables to avoid broadcasting args to parpool
    niter = args.niter;
    
    % Recon data
    parfor framen = 1:nframes
        data = reshape(raw(framen,:,:,:,:),[],1);
        
        % Recon preconditioner using conjugate-phase
        im_cp = A' * reshape(W * data, [], 1);
        im_cp = embed(im_cp,true(dim*ones(1,3)));
        im_cp = ir_wls_init_scale(A, data(:), im_cp);
        
        % Recon using preconditioned conjugate gradient (iterative)
        if niter > 0
            im_pcg = qpwls_pcg1(im_cp(true(dim*ones(1,3))), A, 1, data(:), C, ...
                'niter', niter);
            im(:,:,:,framen) = embed(im_pcg,true(dim*ones(1,3)));
            
        else % ...or save image with CP recon
            im(:,:,:,framen) = im_cp;
        end
        
    end
    
end

% Save and return output data
if nargout < 1
    writenii([pdir,'/im_mag'], abs(im), ...
        'fov', fov*ones(1,3), 'tr', tr);
    writenii([pdir,'/im_ang'], angle(im), ...
        'fov', fov*ones(1,3), 'tr', tr);
    fprintf('timeseries saved to file: im*.nii\n');
    clear im imc
else
    fprintf('timeseries not saved to file since data was returned\n'); 
end

end


function raw_corr = kdespike(raw,flinked)

    raw_corr = raw;

    for i = 1:length(flinked) % Loop through linked frame groups
        
        f = flinked{i}; % Get linked frame group
        nspikes = 0;
        
        for j = 1:size(raw,3)*size(raw,4)*size(raw,5) % Loop through leaves, echoes, coils
            
            % Get mean and standard deviation for determining outliers
            mu = mean(abs(raw(f,:,j)),1);
            sig = std(abs(raw(f,:,j)),[],1);
            
            for k = 1:size(raw,2) % Loop through kspace locs
                % Determine indicies of outlier frames
                fo = f(abs(raw(f,k,j)) > mu(k) + 2*sig(k));
                nspikes = nspikes + length(fo);
                
                for ifo = 2:length(fo)-1 % Loop through outlier frames
                    % Replace outlier data with mean of neighboring frames
                    raw_corr(fo(ifo),k,j) = mean(raw(fo(ifo+[-1,1]),k,j),1);
                end
            end
            
        end
        
        fprintf('despiked %d/%d data points from frame group %d\n', ...
            nspikes, numel(raw(f,:)), i);
        
    end

end

function raw_corr = phasedetrend(raw,nnav,pdorder)

% Get dimensions
nframes = size(raw,1);
ndat = size(raw,2);
nechos = size(raw,3);
ntrains = size(raw,4);
ncoils = size(raw,5);

% Determine navigator points to use in lsq fit
navpts = round((ndat - nnav)/2):round((ndat + nnav)/2);

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
                
                % Calculate least squares fit of navigator phase
                echo = raw(framen,:,echon,trainn,coiln);
                y = angle(echo(navpts));
                betas = pinv(A)*y(:);
                
                % Store fit for current index
                fits(framen,:,echon,trainn,coiln) = ...
                    ((1:ndat) - round(ndat/2))'.^(pdorder:-1:0) * betas;
            end
        end
    end
end

% Correct echo by subtracting out fits from phase
raw_corr = raw.*exp(-1i*fits);
    
end