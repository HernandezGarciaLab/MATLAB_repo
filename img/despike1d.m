function [data_corr,nspikes] = despike1d(data,dim,varargin)
% data_corr = despike1d(data,dim,varargin)
%
% Part of umasl project by Luis Hernandez-Garcia and David Frey
% @ University of Michigan 2023
%
% Description: function to despike data along one specified dimension
%   using isoutlier() routine for detection and linear regression for correction
%
% Dependencies:
%   - matlab default path
%       - can be restored by typing 'restoredefaultpath'
%   - umasl
%       - github: fmrifrey/umasl
%       - umasl/matlab/ and subdirectories must be in current path
%
% Static input arguments:
%   - data:
%       - data to despike
%       - array of data of any size
%       - no default, required argument
%   - dim:
%       - dimension to despike along
%       - integer < ndims(data) describing dimension to despike
%       - default is 1 (despike along first dimension)
%
% Variable input arguments (type 'help varargin' for usage info):
%   - 'linked':
%       - groups of indicies to despike
%       - array of indicies or cell array containing arrays of indicies to despike
%       - examples:
%		% despike the 5th-7th indicies of a dataset along dim 1:
%		data_corr = despike1d(data,1,5:7);
%		% despike odd and even indicies seperately along dim 4:
%		data_corr = despike1d(data,4,{1:2:size(data,4),2:2:size(data,4)});
%       - default is all indicies along dim (1:size(data,dim))
%   - 'outliermethod':
%       - 'method' option for isoutlier() function
%       - string describing method to use for outlier detection
%       - default is 'median'
%
% Function output:
%   - data_corr:
%       - corrected (despiked) data
%

% Define function static arg defaults
if nargin < 2 || isempty(dim)
    dim = 1;
end

if iscomplex(data)
    data_corr = despike1d(real(data),dim,varargin{:}) + ...
        1i*despike1d(imag(data),dim,varargin{:});
    return
end

% Define varargin defaults
defaults = struct( ...
    'linked',           1:size(data,dim), ...
    'outliermethod',    'median' ...
    );

% Parse arguments using matlab's built-in variable input parser
args = vararginparser(defaults,varargin{:});

% Determine permutation by swapping first dimension with desired
permutation = 1:ndims(data);
permutation(dim) = 1;
permutation(1) = dim;

% And permute/reshape the data
data = permute(data,permutation);
datashape = size(data);
data = reshape(data,size(data,1),[]);
data_corr = data;
nspikes = 0;

% Loop through sets of indicies
for s = 1:length(args.linked)
    % Determine outliers for set of frames
    data_s = data(args.linked{s},:);
    [badf,badp] = find(isoutlier(abs(data_s),args.outliermethod,1));
    nspikes = nspikes + length(badf);
    
    % Only index unique points that contain outliers
    data_s = data(args.linked{s},unique(badp));

    % Loop through all unique points that contain at least 1 outlier
    parfor bpi = 1:length(unique(badp))

        % Get data for current point
        datap = data_s(:,bpi);

        % Determine & seperate bad/good frames for point
        f1 = badf(badp == bpi);
        f0 = 1:size(datap,1);
        f0(f1) = [];

        % Replace with nearest neighbor
        datap(f1) = interp1(f0,datap(f0),f1,'nearest','extrap');
        data_s(:,bpi) = datap;

    end

    % Apply 
    data_corr(args.linked{s},unique(badp)) = data_s;

end

% Fix permutation and reshape
data_corr = reshape(data_corr,datashape);
data_corr = permute(data_corr,permutation);

end

