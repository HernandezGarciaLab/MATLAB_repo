function im_smooth = smoothim(im,varargin)
% function im_smooth = smoothim(im,varargin)
%
% Part of umasl project by Luis Hernandez-Garcia and David Frey
% @ University of Michigan 2022
%
% Description: Function to smooth an image using a gaussian convolution
%   kernel
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
%
% Static input arguments:
%   - im:
%       - image to smooth
%       - either a float/double 3D image array or name of a .nii file
%       - if passing in a 3D image array without a return val, must also
%           specify fov and tr so image can be saved to file
%       - default is 'timeseries_mag'
%
% Variable input arguments (type 'help varargin' for usage info):
%   - 'tr'
%       - temporal frame repetition time of timeseries
%       - double/float describing tr (ms)
%       - not necessary if reading timeseries from file
%       - must specify (no default) if im is passed as image array
%   - 'fov'
%       - field of view of image
%       - double/float array of size 1x3 describing FOV (cm)
%       - not necessary if reading timeseries from file
%       - must specify (no default) if im is passed as image array
%   - 'fwhm'
%       - full width half max of gaussian smoothing kernel
%       - double/float describing fwhm as fraction of FOV
%       - default is 0.05
%   - 'scaleoutput'
%       - option to scale nii files to full dynamic range
%       - boolean integer (0 or 1) to use or not
%       - type 'help writenii' for more information
%       - default is 1
%   - 'silent'
%       - option for no fprint messages
%       - boolean integer (0 or 1) describing whether or not to use
%       - default is 0
%
% Function output:
%   - im_smooth
%       - smoothed image
%       - array of image dimensions
%       - if image is returned, it will not be saved to file
%       - if image is not returned, it will be saved to smoothed.nii
%

    % Define default arguments
    defaults = struct(...
        'fov',          [], ... % FOV (cm) (if im is not read from file)
        'tr',           [], ... % TR (ms) (if im is not read from file)
        'fwhm',         0.05, ... % FWHM of gaussian kernel as fraction of FOV
        'scaleoutput',  1, ... % Option to scale nii files to full dynamic range
        'silent',       0 ... % Option for no fprints
        );
    
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
        args.tr = h.pixdim(5);
    elseif isempty(args.fov) && nargout < 1
        error('Must specify fov if image is not read from file');
    end
    
    % Get dimensions
    dim = [size(im,1),size(im,2),size(im,3)];
    
    % Make convolution kernel
    [X,Y,Z] = ndgrid( ...
        linspace(-1,1,dim(1)), ...
        linspace(-1,1,dim(2)), ...
        linspace(-1,1,dim(3)));
    K = reshape(exp(-(vecnorm([X(:),Y(:),Z(:)],2,2)/args.fwhm).^2),dim);
    K = K / sum(K(:));
    
    % Convolve
    im_smooth = fftshift( ifftn( fftn(im) .* fftn(K) ) );

    % Save results that aren't returned to nifti:
    if nargout < 1
        
        % Save image
        writenii('./smoothed.nii', im_smooth, ...
            'fov', args.fov, 'tr', 1, 'doscl', args.scaleoutput);
        if ~args.silent
            fprintf('\nSmoothed image saved to smoothed.nii');
        end
        
        % Clear im so it won't be returned
        clear im_smooth;
    else
        if ~args.silent
            fprintf('\nSmoothed image will not be saved to file since it is returned');
        end
    end
    
end

