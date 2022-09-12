function hdr = ge_pfilehdr(strFile,rdbm_rev)
% ge_pfilehdr.m - quick wrapper to GE-supplied header reader functions.
% 
% EXAMPLE
% strFile = 'P01536.7';
% hdr = ge_pfilehdr(strFile)

% $Id: ge_pfilehdr.m 1333 2014-04-22 19:10:45Z klitinas $
% $HeadURL: svn+ssh://klitinas@anger.engin.umich.edu/svn/matlab/img/trunk/ge_pfilehdr.m $

fid = fopen(strFile,'r','l');

if ~exist('rdbm_rev','var')
    ver = fread(fid,1,'float32');
    str = num2str(ver);
    rdbm_rev = str2double(str);
end

% rdb
fseek(fid,0,'bof');
hdr.rdb = read_rdb_hdr(fid,rdbm_rev);

% ps
fseek(fid,hdr.rdb.off_ps,'bof');
hdr.ps = read_psc_header(fid, rdbm_rev);

% exam
fseek(fid,hdr.rdb.off_exam,'bof');
hdr.exam = read_exam_header(fid, rdbm_rev);

% series
fseek(fid,hdr.rdb.off_series,'bof');
hdr.series = read_series_header(fid,rdbm_rev);

% image
fseek(fid,hdr.rdb.off_image,'bof');
hdr.image = read_image_header(fid,rdbm_rev);

% grad (if exists)
if isfield(hdr.rdb,'off_grad_data')
   fseek(fid,hdr.rdb.off_grad_data,'bof');
   hdr.grad = read_grad_header(fid,rdbm_rev);
end

fclose(fid);
