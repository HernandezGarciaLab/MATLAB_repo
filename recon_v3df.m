function [t,im,psf] = recon_v3df(varargin)
    % Define default arguments
    defaults = struct(...
        'raw',          [], ... % Raw data
        'info',         [], ... % Info structure
        'pfile',        'P*.7', ... % Search string for Pfile
        'tol',          1e-5, ... % Kspace distance tolerance
        'itrmax',       15, ... % Max number of iterations for IR
        'frames',       'all', ... % Frames to recon
        'pdorder',      -1, ... % Order of phase detrending poly fit
        'nworkers',     feature('numcores'), ... % Number of workers to use in parpool
        'scaleoutput',  1 ... % Option to scale output to full dynamic range
        );

    % Start timer
    t = tic;
    
%% Set up recon
    % Parse through variable inputs using matlab's built-in input parser
    p = inputParser;    
    parmnames = fieldnames(defaults);
    for i = 1:size(parmnames,1)
        parmname = char(parmnames{i});
        p.addParameter(parmname,defaults.(parmname),@(x)1);
    end
    p.parse(varargin{:});
    args = p.Results;
    
    % Get raw data and info
    if (isempty(args.raw) || isempty(args.info)) % If reading from Pfile
        [raw,info] = readraw_v3df(args.pfile);
        if (isnan(raw) || isnan(info)) % Exit when no Pfiles are found
            t = NaN(1); im = NaN(1); psf = NaN(1);
            return;
        end
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
    if nargout < 2
        % Save timeseries
        im2nii('./timeseries.nii', abs(im), ...
            dim, fov, info.tr, args.scaleoutput);
        fprintf('\nTimeseries saved to timeseries.nii');
    else
        fprintf('\nTimeseries will not be saved to file since it is returned');
    end
    if nargout < 3
        % Save point spread function
        im2nii('./psf.nii', abs(psf), ...
            dim, fov, info.tr, args.scaleoutput);
        fprintf('\nPSF saved to psf.nii');
    else
        fprintf('\nPSF will not be saved to file since it is returned');
    end
    
    % Save and print elapsed time
    t = toc(t);
    fprintf('\nRecon complete. Total elapsed time: %.2fs\n',t);
    
return;

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

return;

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
    
return;

%% im2nii function definition
function im2nii(niifile_name,im,dim,fov,tr,doscl)
    
    % Open nifti file for writing
    [niifile,msg_fopen] = fopen(niifile_name,'wb','ieee-le');
    if ~isempty(msg_fopen), error(msg_fopen); end

    % Determine scaling factors for saving full dynamic range
    if doscl
        DACFactor = 2^15-1;
        y = im;
        y_min = min(y,[],'all'); y_max = max(y,[],'all');
        m = 2*DACFactor/(y_max - y_min);
        x = m*y - (m*y_min + DACFactor);
        scl_inter = (m*y_min + DACFactor)/m;
        scl_slope = 1/m;
        im = x;
    else
        scl_inter = 0;
        scl_slope = 0;
    end
    
    % Define header
    h = struct (...
        ... % Relevant header info:
        'sizeof_hdr',   348, ...
        'dim',          [4, dim, size(im,4), 0, 0, 0], ...
        'pixdim',       [4, fov./dim, tr, 0, 0, 0], ...
        'datatype',     4, ...
        'bitpix',       16, ...
        'scl_inter',    scl_inter, ...
        'scl_slope',    scl_slope, ...
        ...
        ... % Less relevant header info:
        'data_type',        repmat(' ',10,1), ...
        'db_name',          repmat(' ',18,1), ...
        'extents',          0, ...
        'session_error',    1, ...
        'regular',          'r', ...
        'dim_info',         ' ', ...
        'intent_p1',        0, ...
        'intent_p2',        0, ...
        'intent_p3',        0, ...
        'intent_code',      0, ...
        'slice_start',      0, ...
        'vox_offset',       352, ...
        'slice_end',        0, ...
        'slice_code',       ' ', ...
        'xyzt_units',       ' ', ...
        'cal_max',          0, ...
        'cal_min',          0, ...
        'slice_duration',   0, ...
        'toffset',          0, ...
        'glmax',            0, ...
        'glmin',            0, ...
        'descrip',          repmat(' ',80,1), ...
        'aux_file',         repmat(' ',24,1), ...
        'qform_code',       0, ...
        'sform_code',       0, ...
        'quatern_b',        0, ...
        'quatern_c',        0, ...
        'quatern_d',        0, ...
        'qoffset_x',        0, ...
        'qoffset_y',        0, ...
        'qoffset_z',        0, ...
        'srow_x',           zeros(1,4), ...
        'srow_y',           zeros(1,4), ...
        'srow_z',           zeros(1,4), ...
        'intent_name',      repmat(' ',16,1), ...
        'magic',            repmat(' ',1,4), ...
        'originator',       zeros(1,4), ...
        'esize',            0, ...
        'ecode',            0, ...
        'edata',            ' ' ...
    );
    
    % Write header info
    fwrite(niifile, h.sizeof_hdr,       'int32');
    fwrite(niifile, h.data_type,        'ubit8');
    fwrite(niifile, h.db_name,          'ubit8');
    fwrite(niifile, h.extents,          'int32');
    fwrite(niifile, h.session_error,    'int16');
    fwrite(niifile, h.regular,          'ubit8');
    fwrite(niifile, h. dim_info,        'ubit8');
    fwrite(niifile, h.dim,              'int16');
    fwrite(niifile, h.intent_p1,        'float32');
    fwrite(niifile, h.intent_p2,        'float32');
    fwrite(niifile, h.intent_p3,        'float32');
    fwrite(niifile, h.intent_code,      'int16');
    fwrite(niifile, h.datatype,         'int16');
    fwrite(niifile, h.bitpix,           'int16');
    fwrite(niifile, h.slice_start,      'int16');
    fwrite(niifile, h.pixdim,           'float32');
    fwrite(niifile, h.vox_offset,       'float32');
    fwrite(niifile, h.scl_slope,        'float32');
    fwrite(niifile, h.scl_inter,        'float32');
    fwrite(niifile, h.slice_end,        'int16');
    fwrite(niifile, h.slice_code,       'ubit8');
    fwrite(niifile, h.xyzt_units,       'ubit8');
    fwrite(niifile, h.cal_max,          'float32');
    fwrite(niifile, h.cal_min,          'float32'); 
    fwrite(niifile, h.slice_duration,   'float32');
    fwrite(niifile, h.toffset,          'float32');
    fwrite(niifile, h.glmax,            'int32');
    fwrite(niifile, h.glmin,            'int32');
    fwrite(niifile, h.descrip,          'ubit8');
    fwrite(niifile, h.aux_file,         'ubit8');
    fwrite(niifile, h.qform_code,       'int16');
    fwrite(niifile, h.sform_code,       'int16');
    fwrite(niifile, h.quatern_b,        'float32');
    fwrite(niifile, h.quatern_c,        'float32');
    fwrite(niifile, h.quatern_d,        'float32');
    fwrite(niifile, h.qoffset_x,        'float32');
    fwrite(niifile, h.qoffset_y,        'float32');
    fwrite(niifile, h.qoffset_z,        'float32');
    fwrite(niifile, h.srow_x,           'float32');
    fwrite(niifile, h.srow_y,           'float32');
    fwrite(niifile, h.srow_z,           'float32');
    fwrite(niifile, h.intent_name,      'ubit8');
    fwrite(niifile, h.magic,            'ubit8');
    if length(h.magic)==3
        fwrite(niifile, 0,              'ubit8');
    end
    fwrite(niifile, 0.0,                'float32');
    fwrite(niifile, repmat(' ',1,13),   'ubit8');
    
    % Write data
    fseek(niifile, h.vox_offset, 'bof');
    if h.dim(5) >1
        fwrite(niifile, im(:)', 'short');
    else
        fwrite(niifile, im(:), 'short');
    end
    
    fclose(niifile);
    
return;