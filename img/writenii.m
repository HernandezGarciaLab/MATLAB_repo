function writenii(niifile_name,im,varargin)
% function writenii(niifile_name,im,varargin)
%
% Part of umasl project by Luis Hernandez-Garcia and David Frey
% @ University of Michigan 2022
%
% Description: Function to write nii image file from Nd image array
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
%   - niifile_name:
%       - name of nii file to save to
%       - string describing file path/name
%       - if string does not include '.nii', it will be automatically
%           appended
%       - no default, necessary argument
%   - im:
%       - image array
%       - Nd array with image data
%       - no default, necessary argument
%
% Variable input arguments (type 'help varargin' for usage info):
%   - 'hdr':
%       - manual nifti image header
%       - nifti header structure as made with makeniihdr()
%       - if header is passed with other variable inputs, variable inputs
%           will override specified values in header
%       - default is empty
%   - 'fov':
%       - image field of view
%       - array of 1ximage size describing image fov (standard: cm)
%       - if empty, fov from header will be used; if header is also empty,
%           default of [1 1 1] will be used
%       - default is empty
%   - 'tr':
%       - temporal frame repetition time
%       - double/float describing tr (standard: ms)
%       - if empty, tr from header will be used; if header is also empty,
%           default of 1000 ms will be used
%       - default is empty
%   - 'precision':
%       - precision of data type
%       - C/Fortran string format (i.e. 'float32' or 'int16')
%       - default is 'int16' (signed 16 bit integer)
%   - 'savecomplex':
%       - option to save complex data as 2 seperate files
%       - boolean integer (0 or 1) describing whether or not to use
%       - data will be saved as two seperate files: '*_mag.nii' and
%           '*_ang.nii' for magnitude and phase, respectively
%       - 'readcomplex' argument in readnii() will search for these files
%       - default is 0
%   - 'doscl':
%       - option to scale output to full dynamic range in save file
%       - boolean integer (0 or 1) describing whether or not to use
%       - operation makes use of scl_* nifti header fields, which is not
%           supported by some outside functions (all umasl functions
%           support this)
%       - default is 1
%

    % Define default arguments
    defaults = struct(...
        'hdr',            [], ... % header
        'fov',          [], ... % fov (cm)
        'tr',           [], ... % TR (s)
        'precision',    'int16', ... % datatype precision
        'savecomplex',  0, ... % option to save 2 files for complex data
        'doscl',        1 ... % option to scale
        );
    
    % Parse through variable inputs using matlab's built-in input parser
    args = vararginparser(defaults,varargin{:});

    % Remove .nii extension if included
    if contains(niifile_name,'.nii')
        niifile_name = strsplit(niifile_name,'.nii');
        niifile_name = niifile_name{1};
    end
    
    % Get dim
    dim = [size(im,1),size(im,2),size(im,3)];
    
    % Define header
    if isempty(args.hdr)
        [datatype,bitpix] = niidatatype(args.precision);
        h = makeniihdr('datatype', datatype, 'bitpix', bitpix);
        h.dim(1) = ndims(im);
        h.pixdim(1) = ndims(im);
        if isempty(args.fov)
            args.fov = [1 1 1];
        end
        if isempty(args.tr)
            args.tr = 1000;
        end
    elseif ~isequal(args.hdr.dim(2:5), [dim,size(im,4)])
        error('header dimensions do not match array size');
    else % if header is passed
        h = args.hdr;
        if isempty(args.fov)
            args.fov = h.dim(2:4).*h.pixdim(2:4);
        end
        if isempty(args.tr)
            args.tr = h.pixdim(5);
        end
        args.precision = niidatatype(h.datatype);
    end
    
    % Set fov and tr
    h.dim(2:5) = [dim size(im,4)];
    h.pixdim(2:5) = [args.fov./dim args.tr];

    % Determine scaling factors for saving full dynamic range
    if args.doscl
        maxval = 2^(h.bitpix-1)-1;
        y = im;
        y_min = min(y,[],'all'); y_max = max(y,[],'all');
        m = 2*maxval/(y_max - y_min);
        x = m*y - (m*y_min + maxval);
        h.scl_inter = (m*y_min + maxval)/m;
        h.scl_slope = 1/m;
        im = x;
    else
        h.scl_inter = 0;
        h.scl_slope = 1;
    end

    % Recurse if saving complex data
    if iscomplex(im) && args.savecomplex
        writenii([niifile_name,'_mag'],abs(im),'h',h);
        writenii([niifile_name,'_ang'],angle(im),'h',h);
    elseif iscomplex(im)
        warning('Use savecomplex = 1 to save data as mag and phase files');
        im = abs(im);
    end

    % Open nifti file for writing
    [niifile,msg_fopen] = fopen([niifile_name,'.nii'],'wb','ieee-le');
    if ~isempty(msg_fopen), error(msg_fopen); end
    
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
        fwrite(niifile, im(:)', args.precision);
    else
        fwrite(niifile, im(:), args.precision);
    end
    
    fclose(niifile);
end