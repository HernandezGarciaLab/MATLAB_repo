function TF = iscomplex(val)
% function TF = iscomplex(val)
%
% Part of umasl project by Luis Hernandez-Garcia and David Frey
% @ University of Michigan 2022
%
% Description: Quick macro to determine if value is complex
%
%
% Dependencies:
%   - matlab default path
%       - can be restored by typing 'restoredefaultpath'
%   - umasl
%       - github: fmrifrey/umasl
%       - umasl/matlab/ and subdirectories must be in current path
%
% Static input arguments:
%   - val:
%       - value (or array of values) to determine if it is complex
%       - double/float scalar or array
%       - no default, required argument
%
% Function output:
%   - TF:
%       - logical value containing answer
%       - logical scalar or array of same size as val
%

    % Initialize TF
    TF = false(size(val));
    
    % Determine if any value in val has an imaginary value
    for i = 1:length(val(:))
        TF(i) = (imag(val(i)) > 0);
    end
    
end