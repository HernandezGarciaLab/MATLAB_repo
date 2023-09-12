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
    'zoomfac',      1, ... % FOV downscaling factor
    'resfac',       1 ... % Resolution upscaling factor
    );

% Parse through variable inputs using matlab's built-in input parser
args = vararginparser(defaults, varargin{:});

% Get target pfile
pfile = dir(args.pfile);
if isempty(pfile)
    error('no pfiles found with search string: %s', args.pfile);
end
pdir = pfile(1).folder;
pfile = pfile(1).name;

% Load in kspace trajectory & view transformation matrices
ktraj = dir([pdir,'/ktraj*.txt']);
ktraj = load(ktraj(1).name);
kviews = dir([pdir,'/kviews*.txt']);
kviews = load(kviews(1).name);

% Reshape transformation matrices as an array of 3x3 matrices
T = permute(reshape(kviews(:,end-8:end)',3,3,[]),[2,1,3]);

% Load in raw data
[raw,phdr] = readpfile(pfile);
ndat = phdr.rdb.frame_size;
nechoes = phdr.rdb.user2;
ncoils = phdr.rdb.dab(2) - phdr.rdb.dab(1) + 1;
ntrains = phdr.rdb.user1;
nframes = phdr.rdb.user0;
dim = phdr.image.dim_X*args.resfac;
fov = phdr.image.dfov/10/args.zoomfac;

% reshape: ndat x ntrains*nframes x nechoes x 1 x ncoils
%           --> ndat x ntrains x nframes x nechoes x ncoils
raw = reshape(raw,ndat,ntrains,[],nechoes,ncoils);
% permute: ndat x ntrains x nframes x nechoes x ncoils
%           --> nframes x ndat x nechoes x ntrains x ncoils
raw = permute(raw,[3,1,4,2,5]);

% Allocate space for entire trajectory
ktraj_all = zeros(ndat,3,nechoes,ntrains);

% Transform each view
for trainn = 1:ntrains
    for echon = 1:nechoes
        % Index the transformation matrix for current view
        mtxi = (trainn-1)*nechoes + echon;
        
        % Transform the trajectory
        ktraj_all(:,:,echon,trainn) = ktraj*T(:,:,mtxi)';
    end
end

% Apply gradient/sample delay
raw = circshift(raw,[0,args.ndel,0,0,0]);

% Remove unwanted frames
if strcmpi(args.frames,'all')
    args.frames = 1:nframes;
end
raw = raw(args.frames,:,:,:,:);
nframes = length(args.frames);

% Remove ramp points
raw = raw(:,args.nramp+1:end-args.nramp,:,:,:);
ktraj_all = ktraj_all(args.nramp+1:end-args.nramp,:,:,:);

% Remove unwanted echoes
raw = raw(:,:,args.clipechoes(1)+1:end-args.clipechoes(2),:,:);
ktraj_all = ktraj_all(:,:,args.clipechoes(1)+1:end-args.clipechoes(2),:);

% Perform phase detrending
if args.pdorder >= 0
    raw = phasedetrend(raw,args.nnav,args.pdorder);
end

% Format kspace into column vectors
kspace = [reshape(ktraj_all(:,1,:,:),[],1), ...
    reshape(ktraj_all(:,2,:,:),[],1), ...
    reshape(ktraj_all(:,3,:,:),[],1)];

% Create NUFFT object
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
dcf = pipedcf(Gm.Gnufft);
W = Gdiag(dcf(:)./Gm.arg.basis.transform);

% Initialize image matrices
im = zeros([dim*ones(1,3),nframes]);
imc = zeros([dim*ones(1,3),ncoils,nframes]);

% Set smap to uniform if single coil
if ncoils == 1
    smap = ones(dim*ones(1,3));
end

if isempty(args.smap) % Coil-wise recon with RMS coil combo
    
    for framen = 1:nframes
        parfor coiln = 1:ncoils
            % Recon data
            data = reshape(raw(framen,:,:,:,coiln),[],1);
            imc(:,:,:,coiln,framen) = reshape(Gm' * (W*data(:)),dim,dim,dim);
        end
        im(:,:,:,framen) = sqrt(mean(imc(:,:,:,:,framen).^2,4));
    end
    
else % PCG/CP-SENSE recon
    
    % Make regularizer
    R = Reg1(ones(dim*ones(1,3)), 'beta', 2^-12 * numel(raw(1,:,:,:,:))/3, ...
        'mask', true(dim*ones(1,3)));
    C = R.C;
    
    % Incorporate sensitivity encoding into system matrix
    Ac = repmat({[]},ncoils,1);
    for coiln = 1:ncoils
        tmp = smap(:,:,:,coiln);
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
            
            % ...or save image with CP recon
        else
            im(:,:,:,framen) = im_cp;
        end
        
    end
    
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