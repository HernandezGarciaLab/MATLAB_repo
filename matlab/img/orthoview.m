function im_ortho = orthoview(im,varargin)
% function im_ortho = orthoview(im, varargin)
%
% Part of umasl project by Luis Hernandez-Garcia and David Frey
% @ University of Michigan 2022
%
% Description: Function to display 3d images in orthogonal cut view
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
% Static input arguments:
%   - im:
%       - image to display
%       - either a float/double 3D image array or name of a .nii file
%       - no default; required argument
%
% Variable input arguments (type 'help varargin' for usage info):
%   - 'fov':
%       - field of view of image
%       - 1x3 double/float array describing image fov (or at least
%           dimensional ratio, units don't matter as much here)
%       - if reading im from nii file, fov will automatically be read
%       - if empty and not reading im from file, fov will be assumed from
%           dimensions, which may be inaccurate
%       - default is empty
%   - 'frame':
%       - frame of timeseries image to display
%       - integer describing desired frame index
%       - only applicable if size of 4th dimension > 1
%       - default is 1
%   - 'offset':
%       - cut offset for each dimension
%       - 1x3 integer array describing offset from isocenter in each
%           dimension
%       - default is [0,0,0]
%   - 'logscale'
%       - option to display image in logarithmic scale
%       - boolean integer (0 or 1) describing whether or not to use
%       - default is 0
%   - 'ax'
%       - axis to draw orthoviews on
%       - if left empty, function will draw on current axis
%       - if there already exists an image on the current axis which
%           matches the size, the new image will be overlaid
%       - default is empty
%   - 'caxis':
%       - color scale axis bounds
%       - float/double 1x2 array decribing minimum and maximum intensity of
%           color scale
%       - if 'auto' is passed, caxis will use min and max values of image
%       - default is 'auto'
%   - 'colormap:
%       - color map
%       - string describing matlab color map to use
%       - default is 'gray'
%   - 'colorbar':
%       - option to include a colorbar
%       - boolean integer (0 or 1) describing whether or not to use
%       - default is 1
%
% Function output:
%   - im_ortho
%       - orthogonal concatenated image
%       - 2D array of slicewise images (size is dependent on input specs)
%       - no output will be returned unless specified (don't have to use ;
%           to avoid answer being printed to console) 
%       - if returned, ortho image will not be shown
%

    % Define default arguments
    defaults = struct(...
        'fov',      [], ...
        'frame',    1, ...
        'offset',   [0,0,0], ...
        'logscale', 0, ...
        'ax',       [], ...
        'caxis',    'auto', ...
        'colormap', 'gray', ...
        'colorbar', 1 ...
        );
    
    % Parse through variable inputs using matlab's built-in input parser
    args = vararginparser(defaults,varargin{:});
    
    % If im is a nii file name, read in from nii file
    if ischar(im)
        [im,h] = readnii(im);
        args.fov = h.dim(2:4).*h.pixdim(2:4);
    end
    
    % Warn user if image is complex
    if iscomplex(im)
        warning('Complex images are not supported, using absolute value');
        im = abs(im);
    end
    
    % Apply log scale
    if args.logscale
        im = log(im - min(im(:)) + eps());
    end
    
    % Select single frame
    if size(im,4) > 1
        im = im(:,:,:,args.frame);
    end
    
    % Check fov
    if isempty(args.fov)
        warning('FOV is not set, so assuming fov from dimensions (may produce weird images)');
        args.fov = size(im) / max(size(im),[],'all');
    end
    
    % If dimensions are not the same size, inerpolate & extrapolate
    if (size(im,1)~=size(im,2) || size(im,2)~=size(im,3))
        [x,y,z] = ndgrid(...
            linspace(-args.fov(1)/2,args.fov(1)/2,size(im,1)), ...
            linspace(-args.fov(2)/2,args.fov(2)/2,size(im,2)), ...
            linspace(-args.fov(3)/2,args.fov(3)/2,size(im,3)));
        imgrid = griddedInterpolant(x,y,z,im,'nearest','none');
        
        [X,Y,Z] = ndgrid(...
            linspace(-max(args.fov(:))/2,max(args.fov(:))/2,max(size(im))), ...
            linspace(-max(args.fov(:))/2,max(args.fov(:))/2,max(size(im))), ...
            linspace(-max(args.fov(:))/2,max(args.fov(:))/2,max(size(im))));
        
        im = imgrid(X,Y,Z);
    end
    
    % Get cuts
    args.offset = args.offset ./ abs(args.offset) .* ...
        mod(abs(args.offset), round(size(im)/2));
    im_Sag = im(round(size(im,1)/2) + args.offset(1), :, :);
    im_Cor = im(:, round(size(im,2)/2) + args.offset(2), :);
    im_Ax = im(:, :, round(size(im,3)/2) + args.offset(3));
    
    % Concatenate cuts
    im_ortho = [squeeze(im_Sag)', squeeze(im_Cor)', squeeze(im_Ax)'];
    
    if nargout < 1
        
        % Set auto caxis
        if strcmpi(args.caxis,'auto')
            args.caxis = [min(im_ortho(:)),max(im_ortho(:))];
        end
        
        % Get underlay image
        if ~isempty(args.ax)
            uax = args.ax;
            underlay = getimage(uax);
            ucaxis = uax.CLim;
            ucmap = uax.Colormap;
        elseif ~isempty(get(gcf,'CurrentAxes'))
            uax = get(gcf,'CurrentAxes');
            underlay = getimage(uax);
            ucaxis = uax.CLim;
            ucmap = uax.Colormap;
        else
            underlay = zeros(size(im_ortho));
            ucaxis = [0 1];
            ucmap = zeros(64,3);
        end
        
        % Check that underlay image size matches overlay
        if ~isequal(size(underlay),size(im_ortho))
            warning('imcompatible image sizes, cannot overlay on axis');
            underlay = zeros(size(im_ortho));
            ucaxis = [0 1];
            ucmap = zeros(64,3);
        end
        
        % Normalize and clip images
        underlay = (underlay - ucaxis(1)) / diff(ucaxis);
        im_ortho = (im_ortho - args.caxis(1)) / diff(args.caxis);
        underlay(underlay < 0) = 0;
        underlay(underlay > 1) = 1;
        underlay(im_ortho > 0) = 1 + eps;
        im_ortho(im_ortho < 0) = 0;
        im_ortho(im_ortho > 1) = 1;
        
        % Display ortho using imagesc
        imagesc(underlay + im_ortho);
        set(gca,'Ydir','normal');
        grid off
        axis off
        caxis([0 2+eps])
        colormap([ucmap; colormap(args.colormap)]);
        if args.colorbar
            cb = colorbar;
            set(cb, ...
                'XLim', [1,2] + eps, ...
                'XTick', [1,2] + eps, ...
                'XTickLabel', args.caxis);
        end
        
        % Clear im_lb if not returned so it won't be printed to console
        clear im_ortho
        
    end
    
end