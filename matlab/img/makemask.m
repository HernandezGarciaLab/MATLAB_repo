function mask = makemask(im,varargin)
% function im_sub = aslsub(im,varargin)
%
% Part of umasl project by Luis Hernandez-Garcia and David Frey
% @ University of Michigan 2022
%
% Description: Function to perform ASL subtraction on timeseries images
%
%
% Notes:
%   - if output is returned, nii files will not be saved to conserve space
%       and prevent overwriting, and vice verse for when output is not
%       returned (see 'im_sub' under 'Function output')
%   - default values in help message may not be up to date - check defaults
%       structure under the function header
%
% Dependencies:
%   - matlab default path
%       - can be restored by typing 'restoredefaultpath'
%   - umasl
%       - github: fmrifrey/umasl
%       - umasl/matlab/ and subdirectories must be in current path
%   - spm12
%       - github: spm/spm12
%       - spm12/ and subdirectories must be in current path
%
% Static input arguments:
%   - im:
%       - timeseries image to perform subtraction on
%       - either a float/double 3D image array or name of a .nii file
%       - if passing in a 3D image array, must also specify fov and TR
%       - default is 'timeseries_mag'
%
% Variable input arguments (type 'help varargin' for usage info):
%   - 'fov'
%       - field of view of image
%       - double/float array of size 1x3 describing FOV (cm)
%       - not necessary if reading timeseries from file
%       - must specify (no default) if im is passed as image array
%   - 'frame'
%       - frame to use if image is 4D|
%       - int describing frame index
%       - default is 1
%   - 'fwhm'
%       - full width half max of gaussian smoothing kernel
%       - double/float describing fwhm as fraction of FOV
%       - default is 1/32
%   - 'thresh'
%       - threshold as a fraction of standard deviation above mean
%       - double/float describing fraction of std to use for thresholding
%       - default is 3/4
%   - 'show'
%       - option to show image masking
%       - boolean integer (0 or 1) describing whether or not to use
%       - default is 0
%
% Function output:
%   - mask
%       - binary mask image
%       - array of image dimensions
%       - if mask is returned, mask will not be saved to file
%       - if mask is not returned, mask will be saved to mask.nii
%

    % Define default arguments
    defaults = struct(...
        'fov',          [], ... % FOV (cm) (if im is not read from file)
        'frame',        1, ... % Frame to use if im is 4D
        'fwhm',         1/32, ... % FWHM of gaussian kernel as fraction of FOV
        'thresh',       3/4, ... % Mask threshold as frac. of std
        'show',         0 ... % Option to show mask
        );
    
    % Start timer
    t = tic;
    
    % Parse through variable inputs using matlab's built-in input parser
    args = vararginparser(defaults,varargin{:});
    
    % Define default for im
    if nargin < 1 || isempty(im)
        im = 'timeseries_mag';
    elseif iscomplex(im)
        warning('Complex images are not supported, using absolute value');
        im = abs(im);
    end
    
    % If im is a nii file name, read in from file
    if ischar(im)
        [im,h] = readnii(im);
        args.fov = h.dim(2:4) .* h.pixdim(2:4);
    elseif isempty(args.fov)
        error('Must specify fov if image is not read from file');
    end

    % Use only first frame if im is 4D
    if size(im,4) > 1
        im = im(:,:,:,args.frame);
    end
    
    % Smooth image and read in smoothed image
    fprintf('\nSmoothing with a FWHM of %.4f * FOV', args.fwhm);
    im_smooth = zeros(size(im));
    spm_smooth(im,im_smooth,args.fwhm*size(im),4);
    
    % Begin masking
    fprintf('\nMasking with a threshold of mean + %.4f * std', args.thresh);
    
    % Loop through slices along each dimension
    mask = ones(size(im));
    for slicex = 1:h.dim(2)
        imsl = im_smooth(slicex,:,:);
        % Remove pixels with intensity less than threshold
        mask(slicex,:,:) = ...
            (imsl > mean(imsl,'all') + args.thresh*std(imsl,[],'all'));
        % Fill in holes
        mask(slicex,:,:) = imfill(squeeze(mask(slicex,:,:)),'holes');
    end
    % Repeat for y
    for slicey = 1:h.dim(3)
        imsl = im_smooth(:,slicey,:).*mask(:,slicey,:);
        mask(:,slicey,:) = ...
            (imsl > mean(imsl,'all') + args.thresh*std(imsl,[],'all'));
        mask(:,slicey,:) = imfill(squeeze(mask(:,slicey,:)),'holes');
    end
    % Repeat for z
    for slicez = 1:h.dim(4)
        imsl = im_smooth(:,:,slicez).*mask(:,:,slicez);
        mask(:,:,slicez) = ...
            (imsl > mean(imsl,'all') + args.thresh*std(imsl,[],'all'));
        mask(:,:,slicez) = imfill(squeeze(mask(:,:,slicez)),'holes');
    end
    
    % Reject all volumes other than greatest volume
    CC = bwconncomp(mask);
    nvols = length(CC.PixelIdxList);
    volsizes = zeros(nvols,1);
    for r = 1:nvols
        volsizes(r) = length(CC.PixelIdxList{r});
    end
    [~,imaxvol] = max(volsizes);
    mask = zeros(size(im));
    mask(CC.PixelIdxList{imaxvol}) = 1;
    
    % Show masking results
    if args.show
        cfigopen('Image masking');
        lbview(im), hold on
        lbview(im.*mask,'colormap', ...
            [linspace(0,1,64)'.*ones(64,2), ...
            0.5+0.5*(linspace(0,1,64).^0.5)'])
        title(sprintf('Masking results\nfwhm: %.3f, thresh: %.3f', ...
            args.fwhm,args.thresh));
        hold off
    end
    
    % Save results that aren't returned to nifti:
    if nargout < 1
        % Save mask
        writenii('./mask.nii', 1.*mask, ...
            args.fov, 1, 0);
        fprintf('\nMask saved to mask.nii');
        
        % Clear im so it won't be returned
        clear mask;
    else
        fprintf('\nMask will not be saved to file since it is returned');
    end
    
    % Save and print elapsed time
    t = toc(t);
    fprintf('\nMasking complete. Total elapsed time: %.2fs\n',t);
    
end

