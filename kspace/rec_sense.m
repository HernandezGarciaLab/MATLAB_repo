girfFile = '~/matlab/MATLAB_repo/kspace/20230206UHP3T_rbw125.girf';
girfFile = [];

fprintf('\nFirst pass: generating images without sense correction ...');
fprintf('\n but those will be used to create coil sensitivity maps');
recon3dflex3('girf',girfFile,'pdorder',0);

fprintf('\nUsing BART to generate the sensitivity maps ...');
imc = readnii('coils_mag').*exp(1i*readnii('coils_ang'));
smap = bart('ecalib -d0 -m1',fftc(imc,1:3));
writenii('smap_mag',abs(smap));
writenii('smap_ang',angle(smap));

fprintf('\nSecond pass: reconstructing with sensitivity maps from BART ...')
recon3dflex3('smap',smap,'girf',girfFile, 'pdorder',0);

mask = makemask([],'thresh',0.1);

