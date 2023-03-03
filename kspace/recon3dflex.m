function im_all = recon3dflex(varargin)
% function im_all = recon3dflex(varargin)
%
% Part of umasl project by Luis Hernandez-Garcia and David Frey
% @ University of Michigan 2022
%
% Description: Function to reconstruct images from *3dflex ASL sequence
%   using Model-based SENSE recon with Pipe & Menon Density compensation
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
%   - 'pfile'
%       - search string for pfile
%       - string containing search path for desired pfile
%       - will only be used if either 'raw' or 'info' is left blank
%       - type 'help readpfile' for more information
%       - default is 'P*.7' (uses first pfile from directory to read in
%           raw data and info)
%   - 'niter_pcg'
%       - maximum number of iterations for model-based recon
%       - integer describing number of iterations
%       - if 0 is passed, conjugate phase recon will be used
%       - if a value less than 0 is passed, NUFFT recon will be used
%       - default is 0
%   - 'niter_dcf'
%       - maximum number of iterations for dcf calculation
%       - integer describing number of iterations
%       - default is 15
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
%       - if left empty, recon will write 'coils_*.nii' or return it as
%           output so user can make a sense map
%       - default is empty
%   - 'girf'
%       - name of girf file for gradient distortion correction
%       - string describing name of text-based file (can be any extension)
%       - if left empty, no gradient distortion correction will be
%           performed
%       - default is empty
%   - 'despike'
%       - linked frames for despiker in kspace along frames dimension
%       - cell array containing linked frame arrays to depsike
%       - see 'linked' argument of despike1d() for more info
%       - default is empty
%   - 'nramp'
%       - number of ramp points in spiral to delete from data
%       - integer describing number of points in ramp
%       - if 'auto' is passed, ramp points will be determined based on
%           trajectory envelope
%       - default is 'auto'
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
%   - im_all:
%       - output timeseries image, or coil-wise image if smap is not passed
%       - complex array of image dimension
%       - if im is not returned, all images will be saved to nii files
%       - if im is returned, no images will be saved to nii files
%

    % Start parallel pool
    if isempty(gcp('nocreate'))
        parpool(feature('numcores'),'IdleTimeout',90);
    end
    
    % Define default arguments
    defaults = struct(...
        'pfile',        'P*.7', ... % Search string for Pfile
        'niter_pcg',    0, ... % Max number of iterations for IR
        'niter_dcf',    15, ... % Max number of itr for dcf
        'frames',       'all', ... % Frames to recon
        'clipechoes',   [0 0], ... % Number of echoes to remove
        'resfactor',    1, ... % Resolution upsampling factor
        'zoomfactor',   1, ... % FOV zoom factor
        'smap',         [], ... % Sensitivity map for coil combination
        'girf',         '20230206UHP3T_rbw125.girf', ... % girf file to use
        'despike',      {[]}, ... % linked frames to despike
        'nramp',        50, ... % Number of ramp points in spiral traj
        'scaleoutput',  1, ... % Option to scale output to full dynamic range
        'outputtag',    [] ... % Output filename tag
        );
    
    % Start timer
    t = tic;
    
%% Set up data, arguments and info
    % Parse through variable inputs using matlab's built-in input parser
    args = vararginparser(defaults,varargin{:});
    
    % Determine pfile
    d = dir(args.pfile);
    wd = d(1).folder;
    args.pfile = [wd,'/',d(1).name];
    
    % Read in raw data and header from pfile
    fprintf('\nReading raw data from: %s...',args.pfile);
    [raw,phdr] = readpfile(args.pfile);
    fprintf(' Done.');
    
    % Get info from the header:
    ndat = phdr.rdb.frame_size;
    nslices = phdr.rdb.nslices;
    ncoils = phdr.rdb.dab(2) - phdr.rdb.dab(1) + 1;
    nleaves = phdr.image.user0;
    nframes = phdr.rdb.nframes / nleaves;
    tr = phdr.image.tr*1e-3;
    dim = phdr.image.dim_X;
    fov = phdr.image.dfov/10;
    
    % Reshape and permute raw
    raw = reshape(raw,ndat,nleaves,nframes,nslices,ncoils);
    raw = permute(raw,[3,1,2,4,5]);
    
    % Append output tag
    if ~isempty(args.outputtag)
        args.outputtag = ['_',args.outputtag];
    end
    
    % Determine frames to recon 
    if strcmpi(args.frames,'all')
        args.frames = 1:nframes;
    elseif max(args.frames) > nframes
        warning('frames array exceeds total number of frames');
        args.frames = 1:nframes;
    end
    
%% Process Trajectory
    % Get un-transformed kspace trajectory and view transformations
    ks_0 = load([wd '/ktraj_cart.txt']);
    kviews = load([wd '/kviews.txt']);
    
    % Transform trajectory to entire trajectory
    ks = zeros(ndat, 3, nleaves, nslices);
    for leafn = 1:nleaves
        for slicen = 1:nslices
            % Determine indexed transformation matrix
            viewn = (leafn - 1) * nslices + slicen;
            T = reshape(kviews(viewn, end-8:end), 3, 3)';
            
            % Apply transformation to indexed kspace view
            ks(:,:,leafn,slicen) = (T*ks_0')';
        end
    end

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
    raw = raw(:,:,:,1+args.clipechoes(1):nslices-args.clipechoes(2),:);
    ks = ks(:,:,:,1+args.clipechoes(1):nslices-args.clipechoes(2),:);
    nslices = nslices - sum(args.clipechoes);
        
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
        SDlevel = 2;  % Threshold for despoking : N of standard deviations

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
    
    % Remove ramp points
    fprintf('\nRemoving %d ramp points from data...', args.nramp);
    ks([1:args.nramp ndat-args.nramp:ndat],:,:,:) = [];
    raw(:,[1:args.nramp ndat-args.nramp:ndat],:,:,:) = [];
    ndat = ndat - 2*args.nramp;
    fprintf(' Done.');
    
%% Set up reconstruction
    % Vectorize fov/dim and apply interpolation factors
    fov = fov*ones(1,3) / args.zoomfactor;
    dim = round(dim*ones(1,3) * args.resfactor);
    immask = true(dim);
    
    % Create Gmri object
    kspace = [reshape(ks(:,1,:,:),[],1), ...
        reshape(ks(:,2,:,:),[],1), ...
        reshape(ks(:,3,:,:),[],1)];
    kmask = abs(2*pi*vecnorm(kspace.*fov./dim,2,2)) <= pi;
    nufft_args = {dim,...
        6*ones(1,3),...
        2*dim,...
        dim/2,...
        'table',...
        2^10,...
        'minmax:kb'};
    Gm = Gmri(kspace(kmask,:), immask, ...
        'fov', fov, 'basis', {'rect'}, 'nufft', nufft_args(:)');
    niter_pcg = args.niter_pcg; % extract to avoid broadcasting args into parfor
    
    % Create density compensation using Pipe-Menon algorithm
    fprintf('\nCalculating density compensation (Pipe-Menon method, %d itr)...', ...
        args.niter_dcf);
    dcf = pipedcf(Gm.Gnufft,args.niter_dcf);
    W = Gdiag(dcf(:)./Gm.arg.basis.transform, 'mask', kmask);
    fprintf(' Done.')
    
    % Reconstruct point spread function
    fprintf('\nCalculating PSF ...');
    psf = Gm' * (W * ones(sum(kmask),1));
    psf = embed(psf,immask);
    fprintf(' Done.');
    
    % Correct smap for single-coil data
    if isempty(args.smap) && ncoils == 1
        args.smap = ones(dim);
    end
    
    % Make quadratic regularizer for PCG recon
    R = Reg1(ones(dim), 'beta', 2^-12 * sum(kmask)/3, 'mask', immask);
    C = R.C;

    % Use only nufft for fourier encoding if niter < 0
    if args.niter_pcg < 0
        fprintf('\nSetting up NUFFT reconstruction (since niter < 0) for')
        Gm = Gm.Gnufft;
    elseif args.niter_pcg == 0
        fprintf('\nSetting up CP reconstruction model for')
    else
        fprintf('\nSetting up %d-iteration PCG reconstruction model for',args.niter_pcg)
    end
    
    if ~isempty(args.smap) || ncoils == 1
        
        fprintf(' image timeseries (with sensitivity encoding)...');
        
        if isempty(args.smap)
            % If single coil, sensitivity is assumed to be uniform
            args.smap = ones(dim);
            
        elseif ischar(args.smap) && strcmpi(args.smap,'espirit')
            % Estimate sense map with espirit
            fprintf('Estimating sensitivity map using espirit (BART)... ');
            imc = recon3dflex(varargin{:}, ...
                'resfactor', args.resfactor*0.5,'smap',[]);
            args.smap = bart('ecalib -b0 -m1', fftc(imc,1:3));
            args.smap = regrid(args.smap, ... % no zoom factor
                'resfactor', dim./size(args.smap(:,:,:,1)));
            
            if nargout < 1
                % Save sense map to nii file
                writenii([wd,'/smap_mag',args.outputtag], ...
                    abs(args.smap), 'fov', fov, 'tr', 1, 'doscl', 1);
                writenii([wd,'/smap_ang',args.outputtag], ...
                    angle(args.smap), 'fov', fov, 'tr', 1, 'doscl', ...
                    args.scaleoutput);
                fprintf('SENSE map saved to smap_*%s.nii',args.outputtag);
            end
        else    
           % Correct size of smap and apply zoom factor if needed
            args.smap = regrid(args.smap, ...
                'zoomfactor',[args.zoomfactor*ones(1,3),1], ...
                'resfactor', dim./size(args.smap(:,:,:,1))); 
        end
        
        % Correct size of W
        W = Gdiag(repmat(dcf(:)./Gm.arg.basis.transform,1,ncoils), 'mask', kmask);
        
        % Incorporate sensitivity encoding into system matrix
        Ac = repmat({[]},ncoils,1);
        for coiln = 1:ncoils
            tmp = args.smap(:,:,:,coiln);
            tmp = Gdiag(tmp(immask),'mask',immask);
            Ac{coiln} = Gm * tmp;
        end
        A = block_fatrix(Ac, 'type', 'col');
        
        % Reshape data for recon (ndata x ncoils x nframes)
        data_all = reshape(permute(raw(args.frames,:,:,:,:),[2:5,1]), ...
            [],ncoils,length(args.frames));
        
        % Initialize output image (dim x nframes)
        im_all = zeros([dim,length(args.frames)]);
        
    else
        fprintf(' coil-wise images for frame 1 (no sensitivity encoding)...');

        % System matrix is only fourier encoding
        A = Gm;
        
        % Reshape data for recon (ndata x 1 x ncoils)
        data_all = reshape(raw(1,:,:,:,:),[],1,ncoils);
        
        % Initialize output image (dim x ncoils)
        im_all = zeros([dim,ncoils]);
        
    end
    data_all = data_all(kmask,:,:);
    
    fprintf(' Done.')
    
%% Reconstruct images
    fprintf('\nReconstructing...');

    % Loop through image series (frames or coils index)
    parfor n = 1:size(im_all,4)
        
        % Get indexed data
        data_n = data_all(:,:,n);
        
        % Recon preconditioner using conjugate-phase
        im_cp = A' * reshape(W * data_n, [], 1);
        im_cp = embed(im_cp,immask);
        im_cp = ir_wls_init_scale(A, data_n(:), im_cp);
        
        % Recon using preconditioned conjugate gradient (iterative)
        if niter_pcg > 0
            im_pcg = qpwls_pcg1(im_cp(immask), A, 1, data_n(:), C, ...
                'niter', niter_pcg);
            im_all(:,:,:,n) = embed(im_pcg,immask);
            
        % ...or save image with CP recon
        else
            im_all(:,:,:,n) = im_cp;
        end
        
    end
    fprintf(' Done.');
    
%% Save and return output data
    if nargout < 1
        if ~isempty(args.smap)
            % Save timeseries
            writenii([wd,'/timeseries_mag',args.outputtag], abs(im_all), ...
                'fov', fov, 'tr', tr*nleaves, 'doscl', args.scaleoutput);
            writenii([wd,'/timeseries_ang',args.outputtag], angle(im_all), ...
                'fov', fov, 'tr', tr*nleaves, 'doscl', args.scaleoutput);
            fprintf('\nTimeseries saved to timeseries_*%s.nii',args.outputtag);
        else
            % Save coil images
            writenii([wd,'/coils_mag',args.outputtag], abs(im_all), ...
                'fov', fov, 'doScl', args.scaleoutput);
            writenii([wd,'/coils_ang',args.outputtag], angle(im_all), ...
                'fov', fov, 'doScl', args.scaleoutput);
            fprintf('\nCoil images (frame 1) saved to coil_*%s.nii',args.outputtag);
        end
        
        % Save point spread function
        writenii([wd,'/psf',args.outputtag], abs(psf), ...
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

%%
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
