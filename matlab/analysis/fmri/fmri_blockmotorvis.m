%% Set experimental parameters
workFile = 'timeseries.nii'; % Name of timeseries file to use
dosub = 1; % option to do subtraction in this script
frame1 = 5; % if asl, first frame to use in subtraction
tstart = 9; % since ignoring 1st 2 asl frames
toff = 30; % 
ton = 30;
gFWHM = 0.05; % Gaussian smoothing filter full width half max (as fraction of fov)

%% Format nifti for spm operations
% Remove scaling parameters since it does not work with spm and get fov/tr
[im,h] = readnii(workFile);
fov = h.dim(2:4).*h.pixdim(2:4)*10; % in mm
tr = h.pixdim(5)*1e-3; % in s
writenii(workFile,im,fov,tr*1e3,0);

%% Realignment/reslicing using spm
spm_realign(workFile);
spm_reslice(workFile);

%% Perform subtraction
if dosub
    aslsub(['r' workFile],'fstart',frame1,'order',1);
    workFile = 'sub.nii';
    !mv sub.nii rsub.nii
end

%% Smoothing using spm
spm_smooth(['r' workFile],['s' workFile],gFWHM*max(fov(:))*ones(1,3),4);

%% Make design matrix
[ts,h] = readnii(['r' workFile]);
[x,x_hires] = blockstim(h.dim(5),tstart,toff,ton,h.pixdim(5));