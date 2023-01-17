function [im_clean,A_noise] = compcor(im,varargin)
% function [im_clean,A_noise] = compcor(im,varargin)
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
%   - 'thresh_std':
%       - noise masking threshold
%       - float/double in range [0 100] describing top percentage of noise
%           to reject based on temporal standard deviation
%       - default is 30
%   - 'thresh_rho':
%       - regressor correlation rejection threshold
%       - float/double in range [0 1] describing maximum allowed
%           correlation between passed in design matrix and noise
%           regressors
%       - default is 0.5
%   - 'N':
%       - number of principal component regressors to return
%       - integer describing number of components
%       - default is 10
%   - 'A':
%       - design matrix for fmri experiment (before noise regressors)
%       - 2D float/double array containing regressors of size (number of
%           temporal frames x number of regressors)
%       - if passed, function will determine correlation between regressors
%           in A_noise to regressors in A; if the regressors are
%           correlated enough, they will be removed from A_noise
%       - A_noise will not contain regressors from A
%       - default is empty
%   - 'show':
%       - option to show images
%       - boolean integer (0 or 1) describing whether or not to show
%       - default is 0
%
% Function output:
%   - im_clean:
%       - image after removing principal noise components
%       - float/double array of same dimensions as im
%       - if not returned, im_clean will be saved to compcorclean.nii
%       - if returned, im_clean will not be saved to file
%   - A_noise:
%       - constructed design matrix of noise regressors
%       - double/float column vector of (nframes x N)
%

    % Define default arguments
    defaults = struct(...
        'thresh_std',   30, ... % Top percent of std to threshold
        'thresh_rho',   0.5, ... % Highest allowed corellation of noise to regressor
        'N',            10, ... % Number of components to analyze
        'A',            [], ... % Design matrix for fmri experiment
        'show',         0 ...   % option to show images
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
        im = readnii(im);
    end
    
    % Get std map & mask
    sigma = std(im,[],4);
    sigma_ordered = sort(sigma(sigma(:)>0));
    Nintgrl = cumsum(sigma_ordered)/sum(sigma_ordered) * 100;
    thval = sigma_ordered(find(Nintgrl>100-args.thresh_std,1));
    msk = 1*(sigma > thval);
    
    if args.show
        % Show mask and std map
        cfigopen('CompCor std map');
        subplot(2,1,1)
        lbview(sigma)
        title('Standard deviation');
        subplot(2,1,2);
        lbview(msk)
        title(sprintf('Mask (top %.1f%% of std)',args.thresh_std));
    end
    
    % Get eigenvalues and noise components from masked image
    im_msk = msk.*im;
    [~,s,v] = svd(reshape(im_msk,[],size(im,4)),0);
    A_noise = v(:,1:args.N);
    A_noise = A_noise - mean(A_noise,1); % mean center
    
    if args.show
        % Show eigenvalues and noise components
        cfigopen('CompCor noise components');
        subplot(2,1,1)
        plot(diag(s))
        title('Eigenvalues')
        subplot(2,1,2)
        plot(A_noise)
        title(sprintf('First %d noise components',args.N));
    end
    
    if ~isempty(args.A)
        % Decorrelate noise components from the design matrix
        badinds = [];
        for n = 1:args.N
            design = args.A*pinv(args.A) * A_noise(:,n);
            rho = corrcoef(A_noise(:,n),design);
            if abs(rho(1,2)) > args.thresh_rho
                badinds = [badinds; n];
            end
        end
        A_noise(:,badinds) = [];
    end
    
    % Clean up timeseries by removing the noise
    bhat = pinv(A_noise) * reshape(im,[],size(im,4))';
    im_clean = im - reshape((A_noise*bhat)',size(im));
    
    if args.show
        % Show std before and after
        cfigopen('CompCor noise reduction');
        subplot(2,2,1)
        lbview(std(im,[],4));
        title('Standard deviation (before)')
        subplot(2,2,2)
        lbview(std(im_clean,[],4));
        title('Standard deviation (after)');
        
        % Calculate and show BIC for noise model
        RSS = sum(im_clean.^2,4);
        n = size(im_clean,4);
        k = size(A_noise,2);
        BIC = n * log(RSS/n) + k * log(n);
        BIC = BIC .* msk;
        BIC(BIC==0) = NaN;
        subplot(2,2,3);
        histogram(BIC(:),50);
        title(sprintf('BIC histogram, mean = %.2f', mean(BIC(~isnan(BIC(:))))));
        subplot(2,2,4);
        lbview(BIC);
        title('voxel-wise BIC');
    end
    
    % Save and print elapsed time
    t = toc(t);
    fprintf('\nCompCor completed. Elapsed time: %.2fs\n',t);
    
end

