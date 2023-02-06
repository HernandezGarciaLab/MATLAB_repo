function mat2txt(fname,A)
% function mat2txt(fname,A)
%
% Part of umasl project by Luis Hernandez-Garcia and David Frey
% @ University of Michigan 2023
%
% Description: function to quickly write a matrix to text file
%
% Dependencies:
%   - matlab default path
%       - can be restored by typing 'restoredefaultpath'
%   - umasl
%       - github: fmrifrey/umasl
%       - umasl/matlab/ and subdirectories must be in current path
%
% Static input arguments:
%   - fname:
%       - name of text file to write matrix to
%       - string describing name of text file
%       - if '.txt' extension is not included, function will append it
%       - no default, argument is required
%   - A:
%       - matrix to write to file
%       - 2D matrix containing float/double data
%       - no default, argument is required
%

    % Add .txt extension if user left it out
    if ~contains(fname,'txt')
        fname = [fname '.txt'];
    end

    % Open file for writing
    fID = fopen(fname,'w');

    % Write matrix as float table
    for row = 1:size(A,1)
        for col = 1:size(A,2)
            fprintf(fID,'%f \t',A(row,col));
        end
        fprintf(fID,'\n');
    end

    % Close the file
    fclose(fID);
    
end

