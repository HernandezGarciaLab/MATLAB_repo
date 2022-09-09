function [raw,info] = readraw_v3df(searchstr)

    % Set default for search string
    if nargin < 1 || isempty(searchstr)
        searchstr = 'P*.7';
    end

    % Find Pfile based on search string
    dirp = dir(searchstr);
    if size(dirp,1) > 1
        fprintf('\nMultiple Pfiles found for search string %s:',searchstr);
        for i = 1:size(dirp,1)
            fprintf('\n\t%s',dirp(i).name);
        end
        fprintf('\n--> Only continuing with first Pfile...');
    elseif size(dirp,1) < 1
        fprintf('\nNo Pfiles found for search string %s:',searchstr);
        raw = NaN(1); info = NaN(1);
        return;
    end
    pfile_name = dirp(1).name;
    fprintf('\nReading Pfile: %s', pfile_name);
    
    % Read header using GE wrapper
    h = ge_pfilehdr(pfile_name);
    
    % Create info struct based on header info
    info = struct(...
        'ndat',     h.rdb.da_xres, ... % Number of points / echo
        'nleaves',  h.image.user0, ... % Number of interleaves
        'nframes',  h.image.user1, ... % Number of temporal frames
        'nslices',  h.image.slquant, ... % Number of slices
        'ncoils',   h.rdb.dab(2) - h.rdb.dab(1) + 1, ... % Number of coils
        'tr',       h.image.tr, ... % TR (usec)
        'te',       h.image.te, ... % TE (usec)
        'xydim',    h.image.dim_X, ... % Image x/y dimension
        'xyfov',    h.image.dfov/10, ... % FOV (cm)
        'slthick',  h.image.slthick/10 ... % Slice Thickness (cm)
        );
    
    % Open pfile
    [pfile,msg_fopen] = fopen(pfile_name,'r','ieee-le');
    if ~isempty(msg_fopen), error(msg_fopen); end
    if fseek(pfile, h.rdb.off_data, 'bof'), error('BOF not found\n'); end

    % Read in data
    raw = zeros(info.nframes,info.ndat,info.nleaves,info.nslices,info.ncoils);
    for coiln = 1:info.ncoils
        for slicen = 1:info.nslices
            % Read in baseline
            fread(pfile, 2*info.ndat, 'short');

            % Read in data
            dat = fread(pfile, [2*info.ndat info.nframes*info.nleaves], 'short');
            dat = reshape(dat(1:2:end,:) + 1i * dat(2:2:end,:), ...
                info.ndat, info.nleaves, info.nframes);

            % Store in raw
            raw(:,:,:,slicen,coiln) = permute(dat,[3 1 2]);

        end
    end
    
    % Close pfile
    fclose(pfile);

return;