function lbview(im, varargin)
% function lbview(im, varargin)
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

    % Define default arguments
    defaults = struct(...
        'frame',    1, ...
        'slices',   'all', ...
        'logscale', 0, ...
        'nrows',    'auto', ...
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
                im_lb_row = [im_lb_row, NaN(size(im,1),size(im,2))];
            end
        end
        im_lb = [im_lb; im_lb_row];
    end
    
    % Display lightbox using imagesc
    imagesc(im_lb);
    set(gca,'Ydir','normal');
    grid off
    axis off
    if args.colorbar
        colorbar;
    end
    colormap(args.colormap);
    caxis(args.caxis);
    
end
