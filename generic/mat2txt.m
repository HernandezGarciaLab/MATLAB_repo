function mat2txt(fname,A)

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

