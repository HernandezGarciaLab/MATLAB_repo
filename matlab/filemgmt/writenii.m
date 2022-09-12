function writenii(niifile_name,im,fov,tr,doscl)

    % Add .nii extension if user left it out
    if ~contains(niifile_name,'.nii')
        niifile_name = [niifile_name '.nii'];
    end
    
    % Open nifti file for writing
    [niifile,msg_fopen] = fopen(niifile_name,'wb','ieee-le');
    if ~isempty(msg_fopen), error(msg_fopen); end

    % Determine scaling factors for saving full dynamic range
    if doscl
        DACFactor = 2^15-1;
        y = im;
        y_min = min(y,[],'all'); y_max = max(y,[],'all');
        m = 2*DACFactor/(y_max - y_min);
        x = m*y - (m*y_min + DACFactor);
        scl_inter = (m*y_min + DACFactor)/m;
        scl_slope = 1/m;
        im = x;
    else
        scl_inter = 0;
        scl_slope = 0;
    end
    
    % Get dim
    dim = [size(im,1),size(im,2),size(im,3)];
    
    % Define header
    h = makeniihdr(...
        'dim',          [4, dim, size(im,4), 0, 0, 0], ...
        'pixdim',       [4, fov./dim, tr, 0, 0, 0], ...
        'datatype',     4, ...
        'bitpix',       16, ...
        'scl_inter',    scl_inter, ...
        'scl_slope',    scl_slope);
    
    % Write header info
    fwrite(niifile, h.sizeof_hdr,       'int32');
    fwrite(niifile, h.data_type,        'ubit8');
    fwrite(niifile, h.db_name,          'ubit8');
    fwrite(niifile, h.extents,          'int32');
    fwrite(niifile, h.session_error,    'int16');
    fwrite(niifile, h.regular,          'ubit8');
    fwrite(niifile, h. dim_info,        'ubit8');
    fwrite(niifile, h.dim,              'int16');
    fwrite(niifile, h.intent_p1,        'float32');
    fwrite(niifile, h.intent_p2,        'float32');
    fwrite(niifile, h.intent_p3,        'float32');
    fwrite(niifile, h.intent_code,      'int16');
    fwrite(niifile, h.datatype,         'int16');
    fwrite(niifile, h.bitpix,           'int16');
    fwrite(niifile, h.slice_start,      'int16');
    fwrite(niifile, h.pixdim,           'float32');
    fwrite(niifile, h.vox_offset,       'float32');
    fwrite(niifile, h.scl_slope,        'float32');
    fwrite(niifile, h.scl_inter,        'float32');
    fwrite(niifile, h.slice_end,        'int16');
    fwrite(niifile, h.slice_code,       'ubit8');
    fwrite(niifile, h.xyzt_units,       'ubit8');
    fwrite(niifile, h.cal_max,          'float32');
    fwrite(niifile, h.cal_min,          'float32'); 
    fwrite(niifile, h.slice_duration,   'float32');
    fwrite(niifile, h.toffset,          'float32');
    fwrite(niifile, h.glmax,            'int32');
    fwrite(niifile, h.glmin,            'int32');
    fwrite(niifile, h.descrip,          'ubit8');
    fwrite(niifile, h.aux_file,         'ubit8');
    fwrite(niifile, h.qform_code,       'int16');
    fwrite(niifile, h.sform_code,       'int16');
    fwrite(niifile, h.quatern_b,        'float32');
    fwrite(niifile, h.quatern_c,        'float32');
    fwrite(niifile, h.quatern_d,        'float32');
    fwrite(niifile, h.qoffset_x,        'float32');
    fwrite(niifile, h.qoffset_y,        'float32');
    fwrite(niifile, h.qoffset_z,        'float32');
    fwrite(niifile, h.srow_x,           'float32');
    fwrite(niifile, h.srow_y,           'float32');
    fwrite(niifile, h.srow_z,           'float32');
    fwrite(niifile, h.intent_name,      'ubit8');
    fwrite(niifile, h.magic,            'ubit8');
    if length(h.magic)==3
        fwrite(niifile, 0,              'ubit8');
    end
    fwrite(niifile, 0.0,                'float32');
    fwrite(niifile, repmat(' ',1,13),   'ubit8');
    
    % Write data
    fseek(niifile, h.vox_offset, 'bof');
    if h.dim(5) >1
        fwrite(niifile, im(:)', 'short');
    else
        fwrite(niifile, im(:), 'short');
    end
    
    fclose(niifile);

end