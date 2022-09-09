function [im,h] = readnii(niifile_name)
    
    % Add .nii extension if user left it out
    if ~contains(niifile_name,'.nii')
        niifile_name = [niifile_name '.nii'];
    end

    % Open nifti file for reading
    [niifile,msg_fopen] = fopen(niifile_name,'r','ieee-le');
    if ~isempty(msg_fopen), error(msg_fopen); end
    
    % Read in header
    h = makeniihdr(...
        'sizeof_hdr',       fread(niifile, 1,'int32')', ...	% should be 348!
        'data_type',        fread(niifile,10,'*char')', ...
        'db_name',          fread(niifile,18,'*char')',...
        'extents',          fread(niifile, 1,'int32')', ...
        'session_error',    fread(niifile, 1,'int16')', ...
        'regular',          fread(niifile, 1,'*char')', ...
        'dim_info',         fread(niifile, 1,'*char')', ...
        'dim',              fread(niifile,8,'int16')', ...
        'intent_p1',        fread(niifile,1,'float32')', ... 
        'intent_p2',        fread(niifile,1,'float32')', ...
        'intent_p3',        fread(niifile,1,'float32')', ...
        'intent_code',      fread(niifile,1,'int16')', ...
        'datatype',         fread(niifile,1,'int16')', ...
        'bitpix',           fread(niifile,1,'int16')', ...
        'slice_start',      fread(niifile,1,'int16')', ...
        'pixdim',           fread(niifile,8,'float32')', ...
        'vox_offset',       fread(niifile,1,'float32')', ...
        'scl_slope',        fread(niifile,1,'float32')', ...
        'scl_inter',        fread(niifile,1,'float32')', ...
        'slice_end',        fread(niifile,1,'int16')', ...
        'slice_code',       fread(niifile,1,'*char')', ...
        'xyzt_units',       fread(niifile,1,'*char')', ...
        'cal_max',          fread(niifile,1,'float32')', ...
        'cal_min',          fread(niifile,1,'float32')', ...
        'slice_duration',   fread(niifile,1,'float32')', ...
        'toffset',          fread(niifile,1,'float32')', ...
        'glmax',            fread(niifile,1,'int32')', ...
        'glmin',            fread(niifile,1,'int32')', ...
        'descrip',          fread(niifile,80,'*char')', ...
        'aux_file',         fread(niifile,24,'*char')', ...
        'qform_code',       fread(niifile,1,'int16')', ...
        'sform_code',       fread(niifile,1,'int16')', ...
        'quatern_b',        fread(niifile,1,'float32')', ...
        'quatern_c',        fread(niifile,1,'float32')', ...
        'quatern_d',        fread(niifile,1,'float32')', ...
        'qoffset_x',        fread(niifile,1,'float32')', ...
        'qoffset_y',        fread(niifile,1,'float32')', ...
        'qoffset_z',        fread(niifile,1,'float32')', ...
        'srow_x',           fread(niifile,4,'float32')', ...
        'srow_y',           fread(niifile,4,'float32')', ...
        'srow_z',           fread(niifile,4,'float32')', ...
        'intent_name',      fread(niifile,16,'*char')', ...
        'magic',            fread(niifile,4,'*char')', ...
        'originator',       fread(niifile, 5,'int16'),...
        'esize',            0, ...
        'ecode',            0, ...
        'edata',            '' ...
    );

    % Read in data
    im = fread(niifile, prod(h.dim(2:4)), 'short')';
    fclose(niifile);
    
    % Reshape and rescale data
    im = reshape(im,h.dim(2:4));
    im = h.scl_slope*im + h.scl_inter;

end

