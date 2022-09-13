function lbview(im, varargin)

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

