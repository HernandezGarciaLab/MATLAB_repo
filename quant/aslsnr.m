function [tSNR,sSNR] = aslsnr(im,mask)
% function [tSNR,sSNR] = aslsnr(im,mask)
%
% Part of umasl project by Luis Hernandez-Garcia and David Frey
% @ University of Michigan 2022
%
% Description: Function to calculate temporal and spatial SNR for
%   ASL subtraction timeseries
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
%       - image should be a subtraction timeseries (so there is no
%           oscillation between M0/control/label)
%       - either a float/double 3D image array or name of a .nii file
%       - default is 'sub'
%   - mask:
%       - image mask
%       - str describing nii file or binary array of image dimensions
%       - default is 'mask'
%
% Function output:
%   - tSNR:
%       - temporal signal to noise ratio
%       - double/float describing tSNR
%   - sSNR:
%       - spatial signal to noise ratio
%       - double/float describing sSNR
%

    % Define default for im
    if nargin < 1 || isempty(im)
        im = 'sub';
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
    
    % Define default for mask
    if nargin < 2 || isempty(mask)
        mask = 'mask';
    elseif iscomplex(mask)
        warning('Complex images are not supported, using absolute value');
        mask = abs(mask);
    end

    % If mask is a nii file name, read in from file
    if ischar(mask)
        mask = readnii(mask);
    end
    
    % Check dimensions of mask
    if ~isempty(mask) && ~isequal(size(mask),dim)
        error('Mask size does not match image size');
    end    
    
    % Normalize and round mask
    mask = (mask - min(mask(:))) / (max(mask(:)) - min(mask(:)));
    mask = round(mask);
    
    % Calculate noise
    tmean = mean(im,4); % Temporal mean
    tstd = std(im,[],4); % Temporal std
    signal = mean(tmean(mask == 1), 'all');
    tNoise = mean(tstd(mask == 0), 'all');
    sNoise = std(tmean(mask == 0), [], 'all');
    
    % Calculate SNR
    tSNR = signal / tNoise;
    sSNR = signal / sNoise;

end

