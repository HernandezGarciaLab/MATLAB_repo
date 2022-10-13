function im_lb = lbview(im, varargin)
% function im_lb = lbview(im, varargin)
%
% Part of umasl project by Luis Hernandez-Garcia and David Frey
% @ University of Michigan 2022
%
% Description: Function to display 3d images in lightbox view
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
%   - 'frame':
%       - frame of timeseries image to display
%       - integer describing desired frame index
%       - only applicable if size of 4th dimension > 1
%       - default is 1
%   - 'slices':
%       - slices of image to show
%       - either an integer array describing desired slice indices or 'all'
%       - if 'all' is passed, all slices will be shown
%       - default is 'all'
%   - 'logscale'
%       - option to display image in logarithmic scale
%       - boolean integer (0 or 1) describing whether or not to use
%       - default is 0
%   - 'nrows':
%       - number of rows to display slices in
%       - integer describing number of rows
%       - if 'all' is passed, nrows will shape lightbox into a square
%       - default is 'all'
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
%   - im_lb
%       - lightbox concatenated image
%       - 2D array of slicewise images (size is dependent on input specs)
%       - no output will be returned unless specified (don't have to use ;
%           to avoid answer being printed to console) 
%       - if returned, lightbox image will not be shown
%

    % Define default arguments
    defaults = struct(...
        'frame',    1, ...
        'slices',   'all', ...
        'logscale', 0, ...
        'nrows',    'auto', ...
        'ax',       [], ...
        'caxis',    'auto', ...
        'colormap', 'gray', ...
        'colorbar', 1 ...
        );
    
    % Parse through variable inputs using matlab's built-in input parser
    args = vararginparser(defaults,varargin{:});
    
    % If im is a nii file name, read in from file
    if ischar(im)
        im = readnii(im);
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
    
    % Select slices
    if ~strcmp(args.slices,'all')
        im = im(:,:,args.slices);
    end
    
    % Determine number of rows and columns
    if strcmp(args.nrows, 'auto')
        args.nrows = floor(sqrt(size(im,3)));
    end
    ncols = ceil(size(im,3) / args.nrows);
    
    % Reshape image into concatenated slices
    im_lb = [];
    for rown = 1:args.nrows
        im_lb_row = [];
        for coln = 1:ncols
            slicen = (rown-1)*ncols + coln;
            if slicen <=  size(im,3)
                im_lb_row = [im_lb_row, im(:,:,slicen)'];
            else
                im_lb_row = [im_lb_row, min(im(:))*ones(size(im,1),size(im,2))];
            end
        end
        im_lb = [im_lb; im_lb_row];
    end
    
    if nargout < 1
        
        % Set auto caxis
        if strcmpi(args.caxis,'auto')
            args.caxis = [min(im_lb(:)),max(im_lb(:))];
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
            underlay = zeros(size(im_lb));
            ucaxis = [0 1];
            ucmap = zeros(64,3);
        end
        
        % Check that underlay image size matches overlay
        if ~isequal(size(underlay),size(im_lb))
            warning('imcompatible image sizes, cannot overlay on axis');
            underlay = zeros(size(im_lb));
            ucaxis = [0 1];
            ucmap = zeros(64,3);
        end
            
        % Normalize and clip images
        underlay = (underlay - ucaxis(1)) / diff(ucaxis);
        im_lb = (im_lb - args.caxis(1)) / diff(args.caxis);
        underlay(underlay < 0) = 0;
        underlay(underlay > 1) = 1;
        underlay(im_lb > 0) = 1 + eps;
        im_lb(im_lb < 0) = 0;
        im_lb(im_lb > 1) = 1;

        % Display lightbox using imagesc
        imagesc(underlay + im_lb);
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
        clear im_lb
        
    end
    
end

