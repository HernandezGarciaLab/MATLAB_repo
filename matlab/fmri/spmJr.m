function zscore = spmJr(im,A,varargin)
% function spmJr(im,A,varargin)
%
% Part of umasl project by Luis Hernandez-Garcia and David Frey
% @ University of Michigan 2022
%
% Description: Function to solve GLM specified in design matrix A for the
%   parameters by ordinary least squares
%
%
% Notes:
%   - if output (zscore map) is returned, nii files will not be saved to
%       conserve space and prevent overwriting, and vice verse for when
%       output is not returned (see 'im' under 'Function output')
%   - default values in help message may not be up to date - check defaults
%       structure under the function header
%
% Path dependencies:
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
%       - timeseries image to perform analysis on
%       - either a float/double 3D image array or name of a .nii file
%       - if passing in a 3D image array, must also specify fov and TR
%       - default is 'timeseries_mag'
%   - A:
%       - design matrix for GLM
%       - 2D float/double array containing regressors of size (number of
%           temporal frames x number of regressors)
%       - no default, must be specified
%
% Variable input arguments (type 'help varargin' for usage info):
%   - 'C':
%       - contrast matrix for ordinary least squares
%       - 2D double/float array of size (number of contrasts x number of
%           contrasts)
%       - if left empty, contrast matrix will be an identity
%       - default empty, which results in an identity matrix
%   - 'fov'
%       - field of view of image
%       - double/float array of size 1x3 describing FOV (cm)
%       - not necessary if reading timeseries from file
%       - no default, must specify if im is passed as image array
%   - 'mask'
%       - image mask
%       - str describing nii file or binary array of image dimensions
%       - if left empty, image will not be masked
%       - default is empty, which results in a full mask
%   - 'scaleoutput'
%       - option to scale nii files to full dynamic range
%       - boolean integer (0 or 1) to use or not
%       - type 'help writenii' for more information
%       - default is 1
%
% Function output:
%   - zscore:
%       - output timeseries image
%       - complex array of image dimension
%       - if zscore is not returned, images (betas, tscores, zscores) will
%           be saved to nii files
%       - if zscore is returned, will not be saved to nii files
%

    % Define default arguments
    defaults = struct(...
        'C',            [], ... % Contrast matrix
        'fov',          [], ... % FOV (cm) (if im is not read from file)
        'mask',         [], ... % Mask image
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
        args.fov = h.dim(2:4) .* h.pixdim(2:4);
    elseif isempty(args.fov)
        error('Must specify fov if image is not read from file');
    end
    
    % Get and check all dimensions
    dim = [size(im,1),size(im,2),size(im,3)];
    nframes = size(im,4);
    if size(A,1) ~= nframes
        error('First dimension of design matrix must be equal to number of temporal frames');
    end
    ncon = size(A,2);
    if nargin < 3 || isempty(args.C)
        args.C = eye(ncon);
    elseif size(args.C,2) ~= ncon
        warning('2nd dimension of contrast must be equal to second dimension of design matrix, proceeding with an identity matrix');
        args.C = eye(ncon);
    end
    
    % Mask image
    if ischar(args.mask) % if mask points to a file name
        args.mask = readnii(args.mask);
    elseif ~isequal(size(args.mask),dim) % check dimensions
        warning('Mask size does not match image size, proceeding with full mask -- NOT RECOMMENDED!');
        args.mask = ones(dim);
    elseif isempty(args.mask)
        args.mask = ones(dim);
    end
    im = args.mask.*im;
    
    % Estimate beta map by using ordinary least squares
    beta = pinv(A) * reshape(permute(im,[4 1:3]),[],prod(dim));
    beta = permute(reshape(beta,[ncon,dim]),[2:4,1]);
    
    % Calculate temporal variance of data
    V = std(im,[],4).^2;
    
    % Calculate tscores for each contrast
    ncon = size(args.C,1); % Only estimate tscores for desired contrasts
    df = nframes - ncon + 1; % Degrees of freedom
    tscore = zeros([dim,ncon]); % tscore map
    variance = zeros([dim,ncon]); % constrast variance map
    for conn = 1:ncon
        rc = args.C(conn,:) * pinv(A); % Index current contrast in A
        variance(:,:,:,conn) = reshape(rc .* V(:) * rc', dim);
        tscore(:,:,:,conn) = ...
            beta(:,:,:,conn) ./ sqrt(variance(:,:,:,conn) + eps()) .* args.mask;
        if isempty(args.mask)
            % if no mask, try to correct for div0 issues
            tscore(im(:,:,:,1) < min(im(:,:,:,1),[],'all') + ...
                std(im(:,:,:,1),[],'all')) = 0;
        end
    end
    
    % Convert tscores to zscores using spm
    zscore = spm_t2z(tscore,df);
    
    % Save data to files
    if nargout < 1
        % Save beta map   
        writenii('./beta.nii',beta,args.fov,1,args.scaleoutput);
        fprintf('\nBeta estimate map saved to beta.nii');
        
        % Save tscore map     
        writenii('./tscore.nii',tscore,args.fov,1,args.scaleoutput);
        fprintf('\ntscore map saved to tscore.nii');
        
        % Save zscore map      
        writenii('./zscore.nii',zscore,args.fov,1,args.scaleoutput);
        fprintf('\nzscore map saved to zscore.nii');
        
        % Save contrast variance map    
        writenii('./variance.nii',variance,args.fov,1,args.scaleoutput);
        fprintf('\nTemporal contrast variance map saved to variance.nii');
        
        % Clear zscore so it won't be returned as ans
        clear zscore;
    else
        fprintf('\nImages will not be saved to file since zscore is returned');
    end
    
    % Save and print elapsed time
    t = toc(t);
    fprintf('\nspmJr completed. Elapsed time: %.2fs\n',t);
    
end

