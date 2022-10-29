function [A_noise,im_clean] = compcor(im,varargin)
% function [A_noise,im_clean] = compcor(im,varargin)
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
%   - if 2nd output is returned, nii files will not be saved to conserve
%       space and prevent overwriting, and vice verse for when output is
%       not returned (see 'im_clean' under 'Function output')
%   - default values in help message may not be up to date - check defaults
%       structure under the function headerr
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
%       - default is 0.9
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
%   - 'show'
%       - option to show figures for each step of CompCor
%       - boolean integer (0 or 1) describing whether or not to use
%       - default is 0
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
%   - 'tol'
%       - image zero detection tolerance
%       - double/float scalar value describing 0 value roundoff
%       - default is 1e-3
%   - 'scaleoutput'
%       - option to scale nii files to full dynamic range
%       - boolean integer (0 or 1) to use or not
%       - type 'help writenii' for more information
%       - default is 1
%
% Function output:
%   - A_noise:
%       - constructed design matrix of noise regressors
%       - double/float column vector of (nframes x N)
%   - im_clean:
%       - image after removing principal noise components
%       - float/double array of same dimensions as im
%       - if not returned, im_clean will be saved to compcorclean.nii
%       - if returned, im_clean will not be saved to file
%

    % Define default arguments
    defaults = struct(...
        'stdthresh',    0.9, ... % std threshold for noise ROI
        'N',            10, ... % Number of components to analyze
        'mask',         [], ... % Mask image
        'tol',          1e-3, ... % Zero value tolerance        
        'A',            [], ... % Design matrix for fmri experiment
        'show',         0,  ... % Option to show figures
        'tr',           [], ... % TR (ms) (if im is not read from file)
        'fov',          [], ... % FOV (cm) (if im is not read from file)
        'scaleoutput',  1 ... % Option to scale output to full dynamic range
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
        args.tr = h.pixdim(5);
        args.fov = h.dim(2:4) .* h.pixdim(2:4);
    elseif (isempty(args.tr) || isempty(args.fov)) && nargout < 2
        error('Must specify tr and fov if image is not read from file');
    end
    
    % Get dimensions
    dim = [size(im,1),size(im,2),size(im,3)];
    nframes = size(im,4);
    
    % If mask is a nii file name, read in from file
    if ischar(args.mask) % if mask points to a file name
        args.mask = readnii(args.mask);
    end
    
    % Check dimensions of mask
    if ~isempty(args.mask) && ~isequal(size(args.mask),dim)
        warning('Mask size does not match image size, proceeding without mask');
        args.mask = ones(dim);
    end
    
    % Normalize and round mask
    args.mask = (args.mask - min(args.mask(:))) / ...
        (max(args.mask(:)) - min(args.mask(:)));
    args.mask = round(args.mask);
    
    % Determine im with noisy ROI masked
    tstd = std(im,[],4); % temporal standard deviation
    [~,p] = sort(tstd(:)); % sort the std and store indicies
    if ~isempty(args.mask)
        p(args.mask(p) == 0) = []; % ignore inidices outside mask
    else
        p(tstd(p) < args.tol) = []; % ignore indicies where std = 0
    end
    im_noise = im .* (tstd >= tstd(p(round(args.stdthresh*length(p)))));
    
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
    
    % Estimate noise and remove from clean timeseries
    beta_noise = pinv(A_noise) * reshape(permute(im,[4 1:3]),[],prod(dim));
    im_clean = im - ...
        permute(reshape(A_noise*beta_noise, [nframes,dim]),[2:4,1]);
    
    % Save data to files
    if nargout < 2
        % Save cleaned image      
        writenii('./compcorclean.nii',im_clean, ...
            'fov', args.fov, 'tr', args.tr, 'doscl', args.scaleoutput);
        fprintf('\nCleaned timeseries image saved to compcorclean.nii');
    else
        fprintf('\nImages will not be saved to file since sub image is returned');
    end
    
    % Save and print elapsed time
    t = toc(t);
    fprintf('\nCompCor completed. Elapsed time: %.2fs\n',t);
    
end

