function im_sub = aslsub(im,varargin)

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
        im = 'timeseries';
    elseif iscomplex(im)
        warning('Complex images are not supported, using absolute value');
        im = abs(im);
    end
    
    % If im is a nii file name, read in from file
    if ischar(im)
        [im,h] = readnii(im);
        args.tr = h.pixdim(5);
        args.fov = h.dim(2:4) .* h.pixdim(2:4);
    elseif (isempty(args.tr) || isempty(args.fov))
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
        im_sub = 100 * im_sub ./ (M0 + eps());
    elseif args.rel
        warning('Cannot compute relative signal change without M0 frames');
    end
    
    % Save data to files
    if nargout < 1
        % Save subtraction timeseries/mean/std
        mean_sub = mean(im_sub,4); std_sub = std(im_sub,[],4);        
        writenii('./sub.nii',im_sub,args.fov,args.tr,args.scaleoutput);
        fprintf('\nSubtraction timeseries saved to sub.nii');
        writenii('./mean_sub.nii',mean_sub,args.fov,args.tr,args.scaleoutput);
        fprintf('\nTemporal mean subtraction saved to mean_sub.nii');
        writenii('./std_sub.nii',std_sub,args.fov,args.tr,args.scaleoutput);
        fprintf('\nTemporal standard deviation of subtraction saved to std_sub.nii');
        
        % Save tag timeseries/mean/std
        mean_tag = mean(im_tag,4); std_tag = std(im_tag,[],4);        
        writenii('./tag.nii',im_tag,args.fov,args.tr,args.scaleoutput);
        fprintf('\nTag timeseries saved to tag.nii');
        writenii('./mean_tag.nii',mean_tag,args.fov,args.tr,args.scaleoutput);
        fprintf('\nTemporal mean tag saved to mean_tag.nii');
        writenii('./std_tag.nii',std_tag,args.fov,args.tr,args.scaleoutput);
        fprintf('\nTemporal standard deviation of tag saved to std_tag.nii');
        
        % Save control timeseries/mean/std
        mean_con = mean(im_con,4); std_con = std(im_con,[],4);        
        writenii('./con.nii',im_con,args.fov,args.tr,args.scaleoutput);
        fprintf('\nControl timeseries saved to con.nii');
        writenii('./mean_con.nii',mean_con,args.fov,args.tr,args.scaleoutput);
        fprintf('\nTemporal mean control saved to mean_con.nii');
        writenii('./std_con.nii',std_con,args.fov,args.tr,args.scaleoutput);
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