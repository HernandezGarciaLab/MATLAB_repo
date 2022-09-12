function [raw,info] = simsignal_vsasl3dflex(varargin)

    % Define default arguments
    defaults = struct(...
        'phantom',      [], ... % Phantom image (empty for 3D shepp-logan)
        'fov',          [], ... % FOV (cm) (if not reading from Pfile)
        'dim',          [] ... % Image dimensions (if not reading from Pfile)
        );

    % Start timer
    t = tic;
    
    % Parse through variable inputs using matlab's built-in input parser
    args = vararginparser(defaults,varargin{:});
    
    % Read in info structure from pfile
    if (isempty(args.fov) || isempty(args.dim))
        [~,info] = readpfile('./P*.7');
        info.nframes = 1;
        info.ncoils = 1;
    elseif (~isempty(args.fov) && ~isempty(args.dim))
        info = struct(...
            'ndat',     [], ... % Number of points / echo
            'nleaves',  [], ... % Number of interleaves
            'nframes',  1, ... % Number of temporal frames
            'nslices',  [], ... % Number of slices
            'ncoils',   1, ... % Number of coils
            'tr',       4.5e6, ... % TR (usec)
            'te',       [], ... % TE (usec)
            'xydim',    args.dim(1), ... % Image x/y dimension
            'xyfov',    args.fov(1), ... % FOV (cm)
            'slthick',  args.fov(3)/args.dim(3) ... % Slice Thickness (cm)
            );
    else
        error('Must either specify both fov & dim or read it in from a pfile');
    end
    
    % Get un-transformed kspace trajectory and view transformations
    ks_0 = load('ktraj_cart.txt');
    kviews = load('kviews.txt');
    
    % Determine trajectory type from kviews
    isSOS = (all(kviews(:,4)==0) && all(kviews(:,5)==0));
    
    % Update info if not reading from pfile
    if (~isempty(args.fov) && ~isempty(args.dim))
        info.ndat = size(ks_0,1);
        info.nleaves = length(unique(kviews(:,1)));
        info.nslices = length(unique(kviews(:,2)));
        info.te = info.ndat * 4;
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
    
    % Update fov & dim if reading from pfile
    if (isempty(args.fov) || isempty(args.dim))
        % Determine FOV/dim/Nneighbors based on trajectory type
        args.fov = info.xyfov*ones(1,3);
        args.dim = info.xydim*ones(1,3);
        if isSOS
            % For SOS, fix the z fov, dim, and number of neighbors
            args.fov(3) = info.slthick*info.nslices;
            args.dim(3) = info.nleaves*info.nslices;
        end
    end
    
    % Make image space grid
    [X,Y,Z] = ndgrid( ...
        linspace(-args.fov(1)/2, args.fov(1)/2, args.dim(1)), ...
        linspace(-args.fov(2)/2, args.fov(2)/2, args.dim(2)), ...
        linspace(-args.fov(3)/2, args.fov(3)/2, args.dim(3)));
    
    % Get/check phantom
    if isempty(args.phantom)
        args.phantom = phantom3d(max(args.dim));
        [X0,Y0,Z0] = ndgrid( ...
            linspace(-max(args.fov)/2, max(args.fov)/2, max(args.dim)), ...
            linspace(-max(args.fov)/2, max(args.fov)/2, max(args.dim)), ...
            linspace(-max(args.fov)/2, max(args.fov)/2, max(args.dim)));
        interp_phantom = griddedInterpolant(X0,Y0,Z0,args.phantom);
        args.phantom = interp_phantom(X,Y,Z);
    elseif ~isequal(size(args.phantom),args.dim)
        error('Phantom size must match size of image in Pfile: %dx%dx%d', ...
            args.dim(1), args.dim(2), args.dim(3));
    end
    
    % Initialize raw and progress message
    raw = zeros(1,info.ndat,info.nleaves,info.nslices,1);
    fprintf('\nSimulating signal... ');
    msg_simprog = '';
    
    % Loop through views
    for leafn = 1:info.nleaves
        for slicen = 1:info.nslices
            % Print progress
            if ~isempty(msg_simprog)
                fprintf(repmat('\b',1,length(msg_simprog)));
            end
            msg_simprog = sprintf('(leaf %d/%d, slice %d/%d)', ...
                leafn,info.nleaves,slicen,info.nslices);
            fprintf(msg_simprog);
            
            % Use signal eq. to simulate signal
            ks_view = ks(:,:,leafn,slicen);
            raw(1,:,leafn,slicen,1) = exp(-1i * 2 * pi * ks_view * [X(:),Y(:),Z(:)]') * args.phantom(:);
        end
    end
    
    % Save and print elapsed time
    t = toc(t);
    fprintf('\nSignal simulation complete. Total elapsed time: %.2fs\n',t);
    
end