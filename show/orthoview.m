function orthoview(im,varargin)

    % Define default arguments
    defaults = struct(...
        'fov',      [], ...
        'frame',    1, ...
        'offset',   [0,0,0], ...
        'logscale', 0, ...
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
    im_Sag = circshift(im(round(size(im,1)/2), :, :),args.offset(1),1);
    im_Cor = circshift(im(:, round(size(im,2)/2), :),args.offset(2),2);
    im_Ax = circshift(im(:, :, round(size(im,3)/2)),args.offset(3),3);
    
    % Concatenate cuts
    im_ortho = [squeeze(im_Sag)', squeeze(im_Cor)', squeeze(im_Ax)'];
    
    % Display orthoviews using imagesc
    imagesc(im_ortho);
    set(gca,'Ydir','normal');
    grid off
    axis off
    if args.colorbar
        colorbar;
    end
    colormap(args.colormap);
    caxis(args.caxis);
    
end