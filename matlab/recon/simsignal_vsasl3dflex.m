function [raw,info] = simsignal_vsasl3dflex(varargin)
% [raw,info] = simsignal_vsasl3dflex(varargin)
%
% Part of umasl project by Luis Hernandez-Garcia and David Frey
% @ University of Michigan 2022
%
% Description: Function to simulate a raw signal based on trajectory and
%   phantom image
%
%
% Notes:
%   - default values in help message may not be up to date - check defaults
%       structure under the function header
%
% Dependencies:
%   - matlab default path
%       - can be restored by typing 'restoredefaultpath'
%   - umasl
%       - github: fmrifrey/umasl
%       - umasl/matlab/ and subdirectories must be in current path
%
% Variable input arguments (type 'help varargin' for usage info):
%   - 'phantom':
%       - [variableinput1 description]
%       - [variableinput1 notes]
%       - [variableinput1 example (if applicable)]
%       - [variableinput1 default]
%   - '[variableinput2 fieldname]':
%       - [variableinput2 description]
%       - [variableinput2 notes]
%       - [variableinput2 example (if applicable)]
%       - [variableinput2 default]
%   ...
%   - '[variableinputN fieldname]':
%       - [variableinputN description]
%       - [variableinputN notes]
%       - [variableinputN example (if applicable)]
%       - [variableinputN default]
%
% Function output:
%   - [output1]:
%       - [output1 desc]
%       - [output1 notes]
%       - [output1 example (if applicable)]
%   - [output2]:
%       - [output2 desc]
%       - [output2 notes]
%       - [output2 example (if applicable)]
%   ...
%   - [outputN]:
%       - [outputN desc]
%       - [outputN notes]
%       - [outputN example (if applicable)]
%

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
            'dim',      args.dim(1), ... % Image x/y dimension
            'fov',      args.fov(1), ... % FOV (cm)
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
        % Determine FOV/dim
        args.fov = info.fov*ones(1,3);
        args.dim = info.dim*ones(1,3);
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
        clear X0 Y0 Z0
    elseif ~isequal(size(args.phantom),args.dim)
        error('Phantom size must match size of image in Pfile: %dx%dx%d', ...
            args.dim(1), args.dim(2), args.dim(3));
    end
    
    % Initialize raw and progress message
    raw = complex(zeros(1,info.ndat,info.nleaves,info.nslices,1));
    fprintf('\nSimulating signal... ');
    msg_simprog = '';
    
    % Make r
    r = [X(:),Y(:),Z(:)];
    clear X Y Z
    
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
            raw(1,:,leafn,slicen,1) = exp(-1i * 2 * pi * ks_view * r') * args.phantom(:);
        end
    end
    
    % Save and print elapsed time
    t = toc(t);
    fprintf('\nSignal simulation complete. Total elapsed time: %.2fs\n',t);
    
end