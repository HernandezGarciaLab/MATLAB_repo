function varargout = imgrid(fov,dim)
% function varargout = imgrid(fov,dim)
%
% Part of umasl project by Luis Hernandez-Garcia and David Frey
% @ University of Michigan 2022
%
% Description: Function to quickly create Nd cartesian grids based on field
%   of view and dimension
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
%   - fov:
%       - field of view (maximum coordinate value x 2) along each dimension
%       - vector of same number of elements as dimensions in cartesian grid
%       - must be same size as dim
%       - no default
%   - dim:
%       - grid size along each dimension
%       - vector of same number of elements as dimensions in cartesian grid
%       - must be same size as fov
%       - no default
%
% Variable function output (type 'help varargout' for usage info):
%   - Grid arrays
%       - example: [X,Y,Z] = imgrid([100 100 100], [100 100 100])
%           will produce X, Y, and Z grids for ranging from -50 to 50 in
%           100 steps along each dimension
%

    % Check that fov and dim are same size
    if numel(fov) ~= numel(dim)
        error('fov and dim must be same size');
    end

    % Check that number of grid dimensions in output is equal to number of
    % input dimensions
    if nargout ~= numel(fov)
        error('number of arguments out must be equal to number of dimenstions in fov');
    end

    % Loop through dimensions until each is accounted for in grid
    outcode = 'varargout{1}';
    incode = 'linspace(-fov(1)/2,fov(1)/2,dim(1))';
    varargout = cell(1,nargout);
    if numel(fov) > 1
        for n = 2:numel(fov)
            outcode = [outcode, sprintf(',varargout{%d}',n)];
            incode = [incode, sprintf(',linspace(-fov(%d)/2,fov(%d)/2,dim(%d))',n,n,n)];
        end
    end
    eval(sprintf('[%s] = ndgrid(%s);',outcode,incode));

end

