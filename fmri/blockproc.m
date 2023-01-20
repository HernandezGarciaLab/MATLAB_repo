function blockproc(dataname, varargin)
% function blockproc(dataname, varargin)
%
% Part of umasl project by Luis Hernandez-Garcia and David Frey
% @ University of Michigan 2022
%
% Description: Function to quickly perform fmri analysis on block stimulus
%   experiments
%
%
% Notes:
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
%   - dataname:
%       - path to image data (folder that includes timeseries)
%       - string describing the path
%       - no default; required argument
%
% Variable input arguments (type 'help varargin' for usage info):
%   - 'fname':
%       - name of timeseries file to use
%       - string describing name of file
%       - default is 'timeseries_mag'
%   - 'discard':
%       - number of frames to discard from beginning of timeseries
%       - integer describing number of frames
%       - default is 0
%   - 'delay':
%       - delay from first readout (frame) to beginning of first rest time
%       - float/double describing time in seconds
%       - default is 0
%   - 'stimtime':
%       - duration of stimulus
%       - float/double describing time in seconds
%       - default is 30
%   - 'resttime':
%       - duration of rest time
%       - float/double describing time in seconds
%       - default is 30
%   - 'dosub':
%       - option to perform asl sur subtraction on timeseries
%       - boolean integer (0/1) decribing whether or not to use
%       - default is 0
%   - 'realignmode':
%       - option to realign images using spm_realign and spm_reslice
%       - integer in range [0,3] describing which mode
%           - mode 0: no realignment
%           - mode 1: realign images
%           - mode 2: append realignment parameters to design matrix
%           - mode 3: realign images and append parameters
%       - default is 0
%   - 'ncompcor':
%       - number of noise regressors to append to design matrix using
%           compcor
%       - integer describing number of regressors
%       - default is 10
%   - 'zthresh':
%       - range of zscores to show
%       - 1x2 float/double array describing range of zscores
%       - default is [1,5]
%   - 'viewtype':
%       - type of view to display results in
%       - string - either 'lbview' or 'orthoview'
%       - default is 'lbview'
%   - 'viewargs':
%       - arguments to use in viewer
%       - cell array containing variable input parameters for view function
%       - default is empty
%   - 'savepng':
%       - option to save png image
%       - boolean integer (0/1) describing whether or not to save
%       - default is 0
%

    % Define default arguments
    defaults = struct( ...
        'fname',        'timeseries_mag', ...
        'discard',      0, ...
        'delay',        0, ...
        'stimtime',     30, ...
        'resttime',     30, ...
        'dosub',        0, ...
        'realignmode',  0, ...
        'ncompcor',     10, ...
        'zthresh',      [1,5], ...
        'viewtype',     'lbview', ...
        'viewargs',     {{}}, ...
        'savepng',      0 ...
        );
    
    % Parse through variable inputs using matlab's built-in input parser
    args = vararginparser(defaults,varargin{:});

    % save current directory
    curdir = pwd;

    % cd to data directory
    cd(dataname)
    savename = strsplit(dataname,'/');
    savename = savename{end};

    % Add .nii extension if user left it out
    if ~contains(args.fname,'.nii')
        args.fname = [args.fname,'.nii'];
    end

    % read in timeseries
    if ~isfile(args.fname)
        error('%s does not exist',args.fname);
    else
        [im,h] = readnii(args.fname);
    end
    
    % realignment
    rp = [];
    if args.realignmode > 0
        spm_realign(args.fname);
        spm_reslice(args.fname);
        switch args.realignmode
            case 1
                im = readnii(['r',args.fname]);
            case 2
                tmp = strsplit(args.fname,'.');
                rp = load(['rp_',tmp{1},'.txt']);
            case 3
                im = readnii(['r',args.fname]);
                tmp = strsplit(args.fname,'.');
                rp = load(['rp_',tmp{1},'.txt']);
            otherwise
                error('invalid realign mode: %d',args.realignmode);
        end
        savename = [savename,sprintf('_rm%d',args.realignmode)];
    end

    % read in or make mask
    if isfile('mask.nii')
        mask = readnii('mask.nii');
    else
        mask = makemask(im,'thresh',0.1,'silent',1);
    end

    % determine base image
    base = im(:,:,:,1);

    % get TR and nframes from header
    TR = h.pixdim(5)*1e-3;
    fov = h.pixdim(2:4).*h.dim(2:4);
    nframes = h.dim(5);

    % discard frames
    im = im(:,:,:,args.discard+1:end)*1e5;
    rp = rp(args.discard+1:end,:);
    nframes = nframes-args.discard;

    % create stimulation waveform
    x = blockstim(nframes, args.delay, args.resttime, args.stimtime, TR, 0);

    % create design matrix
    A = [ones(size(x(:))), x(:), rp];

    % perform subtraction if specified
    order = 0;
    while args.dosub
        sub = aslsub(im,'fstart',1,'order',order);
        ms = mean(sub,4);
        if mean(ms(mask(:)>0)) > 0
            im = sub;
            args.dosub = 0;
        else
            args.order = 1;
        end
    end

    % use compcor to append noise regressors
    if args.ncompcor > 0
        [~,junk]= compcor(im, 'N', args.ncompcor, 'A', A);
        A = [A,junk];
        savename = [savename,sprintf('_cc%d',args.ncompcor)];
    end

    % smooth data
    for n=1:nframes
        im(:,:,:,n) = smooth3(squeeze(im(:,:,:,n)),'gaussian');
    end

    % create contrast matrix
    C = zeros(3,size(A,2));
    for n=1:3, C(n,n) = 1; end

    % use spmJr to perform regression
    zmap = spmJr(im,A,'C',C,'mask',mask,'fov',fov);

    % use only second regressor (stim. activation)
    zmap = squeeze(zmap(:,:,:,2));

    % show
    cfigopen(savename)
    overlayimages(base,[],zmap,args.zthresh,args.viewtype,args.viewargs{:})
    title([savename,' activation zscores'], 'Interpreter', 'none')

    % cd back
    cd(curdir)

    % save image
    if args.savepng
        curP = get(gcf,'Position');
        set(gcf,'Position',[0,0,500,800]);
        F = getframe(gcf);
        im = frame2im(F);
        [imind,cm] = rgb2ind(im,256);
        imwrite(imind,cm,[savename,'.png']);
        set(gcf,'Position',curP);
    end

end

% overlay images function definition
function overlayimages(im1,cax1,im2,cax2,showtype,varargin)

    % normalize underlay to cax1
    if isempty(cax1)
        lbview(im1)
        cax1 = caxis;
    end
    im1 = (im1 - cax1(1)) / diff(cax1);
    im1(im1 < 0) = 0;
    im1(im1 >= 1) = 1 - eps;
    
    % normalize negative overlay to cax2
    im2n = (-im2 - cax2(1)) / diff(cax2);
    im2n(im2n < 0.1) = 0;
    im2n(im2n >= 1) = 1 - eps;
    im1(im2n >= 0.1) = 1;
    
    % normalize positive overlay to cax2
    im2p = (im2 - cax2(1)) / diff(cax2);
    im2p(im2p < 0.1) = 0;
    im2p(im2p >= 1) = 1 - eps;
    im2n(im2p >= 0.1) = 1;
    im1(im2p >= 0.1) = 1;
    
    % make positive and negative colormap
    myhot = hot;
    myblue = myhot(:,[3 2 1]);
    
    % display
    if nargin < 5 || isempty(showtype) || strcmp(showtype,'lbview')
        lbview(im1 + im2n + im2p, varargin{:}, ...
            'colormap', [gray; myblue; myhot], 'caxis', [0 3]);
    elseif strcmp(showtype,'orthoview')
        orthoview(im1 + im2n + im2p, varargin{:}, ...
            'colormap', [gray; myblue; myhot], 'caxis', [0 3]);
    else
        error('invalid type: %s',showtype);
    end
    
    % set colorbars
    c2n = colorbar;
    set(c2n, ...
        'Location', 'WestOutside', ...
        'Ylim', 1+[0 1-eps], ...
        'Ytick', 1+linspace(0,1-eps,3), ...
        'Yticklabel', -linspace(cax2(1),cax2(2),3));
    c2p = colorbar;
    set(c2p, ...
        'Location', 'EastOutside', ...
        'Ylim', 2+[0 1-eps], ...
        'Ytick', 2+linspace(0,1-eps,3), ...
        'Yticklabel', linspace(cax2(1),cax2(2),3));
    
end