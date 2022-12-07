function out = fftc(in,dim)
% function out = fftc(in,dim)
%
% Part of umasl project by Luis Hernandez-Garcia and David Frey
% @ University of Michigan 2022
%
% Description: Function to calculate Nd fourier transform with
%   fft shifts and scaling accounted for
%
%
% Notes:
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
%   - in:
%       - data to perform fft on
%       - float/double matrix containing data at uniformly spaced points
%       - no default, required argument
%   - dim:
%       - dimension(s) to perform fft along
%       - integer array describing desired dimensions, or 'all' for all
%           dimensions
%       - default is 'all'
%
% Function output:
%   - out:
%       - fourier transformed data
%       - float/double matrix containing data at uniformly spaced
%           frequencies
%

    % Set default dimensions
    if nargin<2 || isempty(dim)
        dim = 1:ndims(in);
    elseif any(dim(:) > ndims(in)) || any(dim(:) < 1)
        error('dimensions out of range');
    end
    
    % Define fourier transform with scaling and shifts
    fftc1d = @(x,d) 1/sqrt(size(x,d))*fftshift(fft(fftshift(x,d),[],d),d);
    
    % Fourier transform along each requested dimension
    out = in;
    for n = dim
        out = fftc1d(out,n);
    end
    
end