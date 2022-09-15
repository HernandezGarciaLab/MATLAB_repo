function args = vararginparser(defaults,varargin)
% function args = vararginparser(name)
%
% Part of umasl project by Luis Hernandez-Garcia and David Frey
% @ University of Michigan 2022
%
% Description: Function that translates a comma seperated list (varargin)
%   into argument structure
%
%
% Notes:
%   - This function is meant to be called by other scripts
%
% Usage example:
%   - if given the function:
%       function myargs = myfun(varargin)
%           
%           % Set defaults structure
%           defaults = struct(...
%               'a', 0, ...
%               'b', 0, ...
%               'c', 0 ...
%           );
%           
%           % Parse variable inputs using vararginparser:
%           myargs = vararginparser(defaults,varargin{:});
%           
%       end
%   - then the following would be returned by myfun('a', 1, 'c', 2):
%       ans = 
%       
%           struct with fields:
% 
%               a: 1
%               b: 0
%               c: 3  
%
% Dependencies:
%   - matlab default path
%       - can be restored by typing 'restoredefaultpath'
%
% Static input arguments:
%   - defaults:
%       - default structure array
%       - must contain all field names and default values
%       - no default, argument is required
%   - varargin (technically static):
%       - comma seperated list of variable input arguments
%       - used by passing in another function's varargin as varargin{:}
%       - no default, argument is required
%

    % Make matlab inputParser object
    p = inputParser;
    
    % Loop through fields and assign to parser
    parmnames = fieldnames(defaults);
    for i = 1:size(parmnames,1)
        parmname = char(parmnames{i});
        p.addParameter(parmname,defaults.(parmname),@(x)1);
    end
    
    % Parse and assign to args
    p.parse(varargin{:});
    args = p.Results;
    
end

