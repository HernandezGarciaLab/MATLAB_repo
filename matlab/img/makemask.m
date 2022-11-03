function mask = makemask(im,varargin)
% function mask = makemask(im,varargin)
%
% Part of umasl project by Luis Hernandez-Garcia and David Frey
% @ University of Michigan 2022
%
% Description: Function to create mask from image
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
%       - image to make mask from
%       - either a float/double 3D image array or name of a .nii file
%       - if passing in a 3D image array without a return val, must also
%           specify fov so image can be saved to file
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
%       - default is 0.05
%   - 'thresh'
%       - threshold as a fraction of standard deviation above mean
%       - double/float describing fraction of std to use for thresholding
%       - default is 0.75
%   - 'silent'
%       - option for no fprint messages
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
        'fwhm',         0.05, ... % FWHM of gaussian kernel as fraction of FOV
        'thresh',       0.75, ... % Mask threshold as frac. of std
        'silent',       0 ... % Option for no fprints
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
    elseif isempty(args.fov) && nargout < 1
        error('Must specify fov if image is not read from file');
    end

    % Use only first frame if im is 4D
    if size(im,4) > 1
        im = im(:,:,:,args.frame);
    end
    
    % Smooth image and read in smoothed image
    if args.fwhm > 0
        if ~args.silent
            fprintf('\nSmoothing with a FWHM of %.4f * FOV', args.fwhm);
        end
        im = smoothim(im,'fwhm',args.fwhm,'silent',1);
    end
    
    % Begin masking
    if ~args.silent
        fprintf('\nMasking with a threshold of mean + %.4f * std', args.thresh);
    end
    
    % Loop through slices along each dimension
    mask = ones(size(im));
    for slicex = 1:size(im,1)
        imsl = im(slicex,:,:);
        % Remove pixels with intensity less than threshold
        mask(slicex,:,:) = ...
            (imsl > mean(im,'all') + args.thresh*std(im,[],'all'));
        % Fill in holes
        mask(slicex,:,:) = imfill(squeeze(mask(slicex,:,:)),'holes');
    end
    % Repeat for y
    for slicey = 1:size(im,2)
        imsl = im(:,slicey,:).*mask(:,slicey,:);
        mask(:,slicey,:) = ...
            (imsl > mean(im,'all') + args.thresh*std(im,[],'all'));
        mask(:,slicey,:) = imfill(squeeze(mask(:,slicey,:)),'holes');
    end
    % Repeat for z
    for slicez = 1:size(im,3)
        imsl = im(:,:,slicez).*mask(:,:,slicez);
        mask(:,:,slicez) = ...
            (imsl > mean(im,'all') + args.thresh*std(im,[],'all'));
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
    
    % Save results that aren't returned to nifti:
    if nargout < 1
        
        % Save mask
        writenii('./mask.nii', 1.*mask, ...
            'fov', args.fov, 'tr', 1, 'doscl', 0);
        if ~args.silent
            fprintf('\nMask saved to mask.nii');
        end
        
        % Clear im so it won't be returned
        clear mask;
    else
        if ~args.silent
            fprintf('\nMask will not be saved to file since it is returned');
        end
    end
    
    % Save and print elapsed time
    t = toc(t);
    if ~args.silent
        fprintf('\nMasking complete. Total elapsed time: %.2fs\n',t);
    end
    
end

