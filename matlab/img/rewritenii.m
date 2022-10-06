function rewritenii(niiread_name,niiwrite_name)
% function rewritenii(niiread_name,niiwrite_name)
%
% Part of umasl project by Luis Hernandez-Garcia and David Frey
% @ University of Michigan 2022
%
% Description: Function to rewrite nii file
%
%
% Notes:
%   - this function can be used to reformat nifti file headers so that is
%       consistent for umasl functions
%           - i.e. spm nifti output files may not have the same header
%               formatting as nifti files written with writenii(), this
%               function will rewrite the file so that header formats match
%
% Dependencies:
%   - matlab default path
%       - can be restored by typing 'restoredefaultpath'
%   - umasl
%       - github: fmrifrey/umasl
%       - umasl/matlab/ and subdirectories must be in current path
%
% Static input arguments:
%   - niiread_name:
%       - name of nii file to read from
%       - string describing file path/name
%       - if string does not include '.nii', it will be automatically
%           appended
%       - no default, necessary argument
%   - niiwrite_name:
%       - name of nii file to write to
%       - string describing file path/name
%       - if string does not include '.nii', it will be automatically
%           appended
%       - if left empty, it will be set to niifile_name
%       - default is empty
%

    % Set output file to input file if not specified
    if nargin < 2 || isempty(niiwrite_name)
        niiwrite_name = niiread_name;
    end
    
    % Read in input file
    [im,h] = readnii(niiread_name);
    
    % Get fov and tr
    fov = h.dim(2:4).*h.pixdim(2:4); % in cm
    tr = h.pixdim(5); % in ms
    
    % Write out output file
    writenii(niiwrite_name,im,fov,tr,1);
    
end

