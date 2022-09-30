%% Set experimental parameters
workFile = 'timeseries_mag.nii'; % Name of timeseries file to use
dosub = 1; % option to do subtraction in this script
frame1 = 5; % if asl, first frame to use in subtraction
tstart = -9; % since ignoring 1st 2 asl frames
toff = 30; % 
ton = 30;
gFWHM = 0.005; % Gaussian smoothing filter full width half max (as fraction of fov)
domask = 1;

%% Read in nii and get important dimensions
[im,h] = readnii(workFile);
fov = h.dim(2:4).*h.pixdim(2:4); % in cm
tr = h.pixdim(5); % in ms
if domask
    if ~(exist('mask','var') == 1)
        mask = readnii('mask');
    end
    im = mask.*im;
    workFile = ['masked_' workFile];
    writenii(workFile,im,fov,tr,1);
end

%% Realignment/reslicing using spm
spm_realign(workFile);
spm_reslice(workFile);
workFile = ['r' workFile];

%% Perform subtraction
if dosub
    aslsub(workFile,'fstart',frame1,'order',1);
    workFile = 'sub.nii';
end

%% Smoothing using spm
spm_smooth(workFile,'tmpsmooth.nii',gFWHM*max(fov(:))*ones(1,3),4);
s_im = readnii('tmpsmooth.nii'); !rm tmpsmooth.nii

%% Make design matrix
[x,x_hires] = blockstim(size(s_im,4),tstart,toff,ton,tr*1e-3,1);
A = [ones(size(s_im,4),1),x];

%% Get beta and tscore values
beta = pinv(A) * reshape(permute(s_im,[4 1:3]),[],prod(h.dim(2:4)));
beta = permute(reshape(beta,[2,h.dim(2:4)]),[2:4,1]);
tscore = beta ./ std(beta,[],1:3);

writenii('beta.nii',beta,fov,tr,1);
writenii('tscore.nii',tscore,fov,tr,1);