%% Read in timeseries, get dimensions
% (this section is NOT skippable)

% **** SET PARAMETERS ****
im_name = 'timeseries_mag.nii'; % set timeseries name
% ************************

% Perform operations
[~,h] = readnii(im_name); % read header
fov = h.dim(2:4) .* h.pixdim(2:4); % get fov in cm
tr = h.pixdim(5); % get TR in ms
nframes = h.dim(5); % get number of frames

%% Realign and reslice image using spm
% (this section is skippable if you don't want to realign/reslice)

% Perform operations
spm_realign(im_name); % realign using spm
spm_reslice(im_name); % reslice using spm
im_name = ['r' im_name]; % set im name to new realigned image
rewritenii(im_name); % rewrite nii file to maintain consistent header fmt

%% Make mask
% (this section is skippable if you already have a mask in mask.nii or
%   don't want to mask at all)

% **** SET PARAMETERS ****
mask_im = '../human_SOS_cbf_4shot_itr1/timeseries_mag.nii'; % name of img to use for mask generation
mask_thresh = 0.55; % intensity threshold for masking
mask_fwhm = 0.1; % fwhm for gaussian smoothing kernel (as fraction of fov)
% ************************

% Perform operations
makemask(mask_im,'fwhm',mask_fwhm,'thresh',mask_thresh) % make mask
mask_name = 'mask.nii';

% Generate figure:
cfigopen('Mask image');
subplot(3,2,1), lbview(im_name); title('Timeseries image');
subplot(3,2,2), lbview('mask'); title('Mask');
subplot(3,2,3:6), lbview(readnii(im_name).*readnii('mask')); title('Masked image');

%% Perform ASL subtraction
% (this section is skippable if timeseries is already a subtraction
%   timeseries)

% **** SET PARAMETERS ****
sub_order = 1; % order of asl subtraction
sub_sur = 1; % option to subtract using 'sur' algorithm (preserves temporal res)
sub_fstart = 5; % first frame to use in subtraction
sub_fend = nframes; % last frame to use in subtraction
% ************************

% Perform operations
aslsub(im_name,'order',sub_order,'sur',sub_sur,'fstart',sub_fstart,...
    'fend',sub_fend); % perform subtraction
im_name = 'sub.nii'; % set im name to new subtraction timeseries image
[~,h] = readnii(im_name); % read header in case tr/nframes changed
tr = h.pixdim(5); % get TR in ms
nframes = h.dim(5); % get number of frames

%% Smooth timeseries using spm
% (this section is skippable if you don't want to smooth)

% **** SET PARAMETERS ****
smooth_fwhm = 2/64; % fwhm for gaussian smoothing kernel (as fraction of fov)
% ************************

% Perform operations
spm_smooth(im_name, ['s' im_name], ...
    smooth_fwhm * max(fov(:)) * ones(1,3), 4); % perform smoothing
im_name = ['s' im_name]; % set im name to new smoothed image
rewritenii(im_name); % rewrite nii file to maintain consistent header fmt

%% Make basic design matrix
% (this section is NOT skippable)

% **** SET PARAMETERS ****
stim_tstart = -18;
stim_toff = 30;
stim_ton = 30;
% ************************

% Perform operations
[x,x_hires] = blockstim(nframes,stim_tstart,stim_toff,stim_ton,...
    tr*1e-3,0); % create activation regressor
baseline = ones(nframes,1); % create baseline regressor
A = [baseline,x]; % make design matrix
C = eye(2); % make contrast matrix

%% Make noise component regressors using CompCor
% (this section is skippable if you don't want to do CompCor)

% **** SET PARAMETERS ****
compcor_stdthresh = 0.9; % std noise threshold (as fraction of max std)
compcor_N = 20; % number of noise principal components to generate
% ************************

% Perform operations
A_noise = compcor(im_name,'stdthresh',compcor_stdthresh,'N',compcor_N,...
    'A',A,'show',1,'mask',mask_name); % perform compcor
A = [A, A_noise]; % add A_noise into design matrix
C = [C,zeros(2,size(A_noise,2))]; % append contrast matrix

%% Estimate parameter maps
% (this section is NOT skippable since it's like the entire point of the
%   script)

% Perform operations
if exist('mask_name','var')
    spmJr(im_name,A,'C',C,'mask',mask_name); % perform spm with mask
else
    spmJr(im_name,A,'C',C); % perform spm without mask
end