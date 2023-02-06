function im_sub = aslsub(im,varargin)
% function im_sub = aslsub(im,varargin)
%
% Part of umasl project by Luis Hernandez-Garcia and David Frey
% @ University of Michigan 2022
%
% Description: Function to perform ASL subtraction on timeseries images
%
%
% Notes:
%   - if output is returned, nii files will not be saved to conserve space
%       and prevent overwriting, and vice verse for when output is not
%       returned (see 'im_sub' under 'Function output')
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
%   - im:
%       - timeseries image to perform subtraction on
%       - either a float/double 3D image array or name of a .nii file
%       - if passing in a 3D image array without a return val, must also
%           specify fov and tr so image can be saved to file
%       - default is 'timeseries_mag'
%
% Variable input arguments (type 'help varargin' for usage info):
%   - 'fstart':
%       - first frame to use in subtraction
%       - integer describing index of first frame
%       - default is 3 (assuming 2 M0 frames)
%   - 'fend':
%       - last frame to use in subtraction
%       - integer describing index of first frame, or 'auto'
%       - if 'auto' is passed, the last frame will be the last frame in
%           timeseries
%       - default is 'auto'
%   - 'order':
%       - order of subtraction
%       - boolean integer (0 or 1) describing order
%       - 0 is control, label, control, label
%       - 1 is label, contrl, label, control
%       - default is 0
%   - 'tr'
%       - temporal frame repetition time of timeseries
%       - double/float describing tr (ms)
%       - not necessary if reading timeseries from file
%       - must specify (no default) if im is passed as image array
%   - 'fov'
%       - field of view of image
%       - double/float array of size 1x3 describing FOV (cm)
%       - not necessary if reading timeseries from file
%       - must specify (no default) if im is passed as image array
%   - 'sur'
%       - option to use "surround" algorithm to preserve temporal
%           resolution
%       - boolean integer (0 or 1) decribing whether or not to use
%       - default is 1
%   - 'rel'
%       - option output subtraction as relative percent signal change to M0
%       - boolean integer (0 or 1) decribing whether or not to use
%       - will return error if fstart <= 1 since it has no M0 frames to use
%       - default is 0
%   - 'scaleoutput'
%       - option to scale nii files to full dynamic range
%       - boolean integer (0 or 1) to use or not
%       - type 'help writenii' for more information
%       - default is 1
%
% Function output:
%   - im_sub
%       - subtraction timeseries image
%       - array of image dimensions
%       - if im_sub is returned, images will not be saved to file
%       - if im_sub is not returned, images (mean, standard deviation, &
%           timeseries of control, label, and subtraction images) will be
%           saved to nii files
%

    % Define default arguments
    defaults = struct(...
        'fstart',       3, ... % First frame to use in subtraction
        'fend',         'auto', ... % Last frame to use in subtraction
        'order',        0, ... % Order of subtraction
        'tr',           [], ... % TR (ms) (if im is not read from file)
        'fov',          [], ... % FOV (cm) (if im is not read from file)
        'sur',          1, ... % Use 'sur' algorithm for subtractions
        'rel',          0, ... % Output subtraction as relative signal change to M0
        'scaleoutput',  1 ... % Option to scale output to full dynamic range
        );
    
    % Start timer
    t = tic;
    
    % Parse through variable inputs using matlab's built-in input parser
    args = vararginparser(defaults,varargin{:});
    
    % Define default for im
    if nargin < 1 || isempty(im)
        im = 'timeseries_mag';
    elseif iscomplex(im)
        warning('Complex images are not supported, using absolute value');
        im = abs(im);
    end
    
    % If im is a nii file name, read in from file
    if ischar(im)
        [im,h] = readnii(im);
        args.tr = h.pixdim(5);
        args.fov = h.dim(2:4) .* h.pixdim(2:4);
    elseif (isempty(args.tr) || isempty(args.fov)) && nargout < 1
        error('Must specify tr and fov if image is not read from file');
    end
    
    % Auto-determine last frame for subtraction
    if strcmpi(args.fend,'auto')
        args.fend = size(im,4);
    end
       
    % Seperate control and label frames
    % Order 0: control, label, control, label
    % Order 1: label, control, label, control
    im_con = im(:,:,:,args.fstart+args.order:2:args.fend);
    im_tag = im(:,:,:,args.fstart+1*~args.order:2:args.fend);
    
    % Perform subtraction
    if ~args.sur % For pairwise algorithm
        fprintf('\nSubtracting timeseries using pairwise algorithm');
        
        % Pairwise subtract label - control
        im_sub = im_tag - im_con;
        
        % TR is not preserved
        args.tr = 2*args.tr;
        
    else % For sur algorithm
        fprintf('\nSubtracting timeseries using sur algorithm');
        
        % Initialize subtraction ts
        im_sub = zeros(size(im(:,:,:,args.fstart:args.fend)));
        
        % Loop through frames
        for framen = args.fstart:args.fend
            
            % Determine sign for current frame
            fsign = (-1)^(framen - args.fstart + 1*~args.order);
            
            % Perform sur subtraction with specific cases for first/last
            % frames
            switch framen
                case args.fstart
                    im_sub(:,:,:,framen - args.fstart + 1) = fsign * ...
                        (im(:,:,:,framen) - im(:,:,:,framen+1));
                case args.fend
                    im_sub(:,:,:,framen - args.fstart + 1) = fsign * ...
                        (im(:,:,:,framen-1) - im(:,:,:,framen));
                otherwise
                    im_sub(:,:,:,framen - args.fstart + 1) = fsign * ...
                        (2*im(:,:,:,framen) - im(:,:,:,framen-1) - im(:,:,:,framen+1));
            end
            
        end
    end
    
    % Determine relative signal change w.r.t M0 frames
    if args.rel && args.fstart > 1
        fprintf('\nReturning subtractions as relative (%%) signal change');
        M0 = mean(im(:,:,:,1:args.fstart-1),4);
        im_sub = div0(100 * im_sub, M0);
    elseif args.rel
        warning('Cannot compute relative signal change without M0 frames');
    end
    
    % Save data to files
    if nargout < 1
        % Save subtraction timeseries/mean/std
        mean_sub = mean(im_sub,4); std_sub = std(im_sub,[],4);        
        writenii('./sub.nii',im_sub, ...
            'fov', args.fov, 'tr', args.tr, 'doscl', args.scaleoutput);
        fprintf('\nSubtraction timeseries saved to sub.nii');
        writenii('./mean_sub.nii',mean_sub, ...
            'fov', args.fov, 'tr', args.tr, 'doscl', args.scaleoutput);
        fprintf('\nTemporal mean subtraction saved to mean_sub.nii');
        writenii('./std_sub.nii',std_sub, ...
            'fov', args.fov, 'tr', args.tr, 'doscl', args.scaleoutput);
        fprintf('\nTemporal standard deviation of subtraction saved to std_sub.nii');
        
        % Save tag timeseries/mean/std
        mean_tag = mean(im_tag,4); std_tag = std(im_tag,[],4);        
        writenii('./tag.nii',im_tag, ...
            'fov', args.fov, 'tr', args.tr, 'doscl', args.scaleoutput);
        fprintf('\nTag timeseries saved to tag.nii');
        writenii('./mean_tag.nii',mean_tag, ...
            'fov', args.fov, 'tr', args.tr, 'doscl', args.scaleoutput);
        fprintf('\nTemporal mean tag saved to mean_tag.nii');
        writenii('./std_tag.nii',std_tag, ...
            'fov', args.fov, 'tr', args.tr, 'doscl', args.scaleoutput);
        fprintf('\nTemporal standard deviation of tag saved to std_tag.nii');
        
        % Save control timeseries/mean/std
        mean_con = mean(im_con,4); std_con = std(im_con,[],4);        
        writenii('./con.nii',im_con, ...
            'fov', args.fov, 'tr', args.tr, 'doscl', args.scaleoutput);
        fprintf('\nControl timeseries saved to con.nii');
        writenii('./mean_con.nii',mean_con, ...
            'fov', args.fov, 'tr', args.tr, 'doscl', args.scaleoutput);
        fprintf('\nTemporal mean control saved to mean_con.nii');
        writenii('./std_con.nii',std_con, ...
            'fov', args.fov, 'tr', args.tr, 'doscl', args.scaleoutput);
        fprintf('\nTemporal standard deviation of control saved to std_con.nii');
        
        % Clear im_sub so it won't be returned as ans
        clear im_sub;
    else
        fprintf('\nImages will not be saved to file since sub image is returned');
    end
    
    % Save and print elapsed time
    t = toc(t);
    fprintf('\nASL Subtraction completed. Elapsed time: %.2fs\n',t);
    
end
