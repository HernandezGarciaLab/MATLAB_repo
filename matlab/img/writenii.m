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
%   - 'h':
%       - manual nifti image header
%       - nifti header structure as made with makeniihdr()
%       - if header is passed with other variable inputs, variable inputs
%           will override specified values in header
%       - default is empty
%   - 'fov':
%       - image field of view
%       - array of 1ximage size describing image fov (standard: cm)
%       - default is [2 2 2]
%   - 'tr':
%       - temporal frame repetition time
%       - double/float describing tr (standard: ms)
%       - default is 1
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
        'h',        [], ... % Raw data
        'fov',      [2 2 2], ... % Info structure
        'tr',       1, ... % Search string for Pfile
        'doscl',    1 ... % Kspace distance tolerance
        );
    
    % Parse through variable inputs using matlab's built-in input parser
    args = vararginparser(defaults,varargin{:});

    % Add .nii extension if user left it out
    if ~contains(niifile_name,'.nii')
        niifile_name = [niifile_name '.nii'];
    end
    
    % Open nifti file for writing
    [niifile,msg_fopen] = fopen(niifile_name,'wb','ieee-le');
    if ~isempty(msg_fopen), error(msg_fopen); end
    
    % Check for complex image
    if iscomplex(im)
        warning('Complex images are not supported, using absolute value');
        im = abs(im);
    end
    
    % Get dim
    dim = [size(im,1),size(im,2),size(im,3)];
    
    % Define header
    if isempty(args.h)
        h = makeniihdr('datatype', 4, 'bitpix', 16);
        h.dim(1) = ndims(im);
        h.pixdim(1) = ndims(im);
    elseif ~isequal(args.h.dim(2:5), size(im))
        error('header dimensions do not match array size');
    else
        h = args.h;
    end
    
    % Set fov if passed
    if ~isempty(args.fov)
        h.dim(2:4) = dim;
        h.pixdim(2:4) = args.fov./dim;
    end
    
    % Set tr if passed
    if ~isempty(args.tr)
        h.dim(5) = size(im,4);
        h.pixdim(5) = args.tr;
    end

    % Determine scaling factors for saving full dynamic range
    if args.doscl
        DACFactor = 2^15-1;
        y = im;
        y_min = min(y,[],'all'); y_max = max(y,[],'all');
        m = 2*DACFactor/(y_max - y_min);
        x = m*y - (m*y_min + DACFactor);
        h.scl_inter = (m*y_min + DACFactor)/m;
        h.scl_slope = 1/m;
        im = x;
    else
        h.scl_inter = 0;
        h.scl_slope = 1;
    end
    
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