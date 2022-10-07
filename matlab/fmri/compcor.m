function A_noise = compcor(im,varargin)
% function A_noise = compcor(im,varargin)
%
% Part of umasl project by Luis Hernandez-Garcia and David Frey
% @ University of Michigan 2022
%
% Description: Function to derive significant principal compenents for
%   noise regions of interest in fmri data based on temporal variance, as
%   presented in Behzadi et. al. (2007)
%
%
% Notes:
%   - default values in help message may not be up to date - check defaults
%       structure under the function header
%
% Path dependencies:
%   - matlab default path
%       - can be restored by typing 'restoredefaultpath'
%   - umasl
%       - github: fmrifrey/umasl
%       - umasl/matlab/ and subdirectories must be in current path
%
% Static input arguments:
%   - im:
%       - timeseries image to perform analysis on
%       - either a float/double 3D image array or name of a .nii file
%       - if passing in a 3D image array, must also specify fov and TR
%       - default is 'timeseries_mag'
%
% Variable input arguments (type 'help varargin' for usage info):
%   - 'stdthresh':
%       - threshold of std for noisy pixels
%       - double/float describing percentile (fraction) of std
%       - default is 0.25
%   - 'N':
%       - number of principal component regressors to return
%       - integer describing number of components
%       - default is 10
%   - 'mask'
%       - image mask
%       - str describing nii file or binary array of image dimensions
%       - if left empty, image will not be masked
%       - default is empty, which results in a full mask
%   - 'A'
%       - design matrix for fmri experiment (before noise regressors)
%       - 2D float/double array containing regressors of size (number of
%           temporal frames x number of regressors)
%       - if passed, function will determine correlation between regressors
%           in A_noise to regressors in A; if the regressors are
%           correlated enough, they will be removed from A_noise
%       - A_noise will not contain regressors from A
%       - default is empty
%
% Function output:
%   - A_noise:
%       - constructed design matrix of noise regressors
%       - double/float column vector of (nframes x N)
%

    % Define default arguments
    defaults = struct(...
        'stdthresh',    0.25, ... % std threshold for noise ROI
        'N',            10, ... % Number of components to analyze
        'mask',         [], ... % Mask image
        'A',            [] ... % Design matrix for fmri experiment
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
        im = readnii(im);
    end
    
    % Get dimensions
    dim = [size(im,1),size(im,2),size(im,3)];
    nframes = size(im,4);
    
    % Mask image
    if ischar(args.mask) % if mask points to a file name
        args.mask = readnii(args.mask);
    elseif ~isequal(size(args.mask),dim) % check dimensions
        warning('Mask size does not match image size, proceeding with full mask');
        args.mask = ones(dim);
    elseif isempty(args.mask)
        args.mask = ones(dim);
    end
    im = args.mask.*im;
    
    % Determine im with noisy ROI masked
    tstd = std(im,[],4); % temporal standard deviation
    im_noise = im .* (tstd >= args.stdthresh * max(tstd(:)));
    
    % Decompose and store singular values of noise
    [~,~,v] = svd(reshape(im_noise,prod(dim),nframes),0);
    A_noise = v(:,1:args.N) - mean(v(:,1:args.N),1);
    
    % Remove correlated regressors if design matrix is passed
    if ~isempty(args.A)
        for n = 1:args.N
            r = A_noise(:,n) - args.A*pinv(args.A)*A_noise(:,n);
            A_noise(:,n) = r - mean(r);
        end
    end
    
end

