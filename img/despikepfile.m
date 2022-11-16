function despikepfile(searchstr,linkedframes,thresh)

    % Set default for search string
    if nargin < 1 || isempty(searchstr)
        searchstr = 'P*.7';
    end
    
    % Set default for linkedframes
    if nargin < 2 || isempty(linkedframes)
        linkedframes = 'all';
    end
    
    % Set default for thresh
    if nargin < 3 || isempty(thresh)
        thresh = 2;
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
        error('No Pfiles found for search string %s:',searchstr);
    end
    pfile = [dirp(1).folder '/' dirp(1).name];
    
    fid = fopen([pfile '_despiked'],'w'); 
    [data,hdr] = readpfile(pfile);
    if ~copyfile(pfile, [pfile '_despiked'])
        error('Could not copy file to %s_despiked',pfile);
    end
    
    % Reformat data
    nframes = hdr.image.user1;
    ncoils = hdr.rdb.dab(2) - hdr.rdb.dab(1) + 1;
    data = reshape(data, nframes, [], ncoils);
    
    % Determine linked frames if set to all
    if strcmpi(linkedframes,'all')
        linkedframes = {1:nframes};
    end
    
    totalrmpts = 0;
    msg_cprog = [];
    
    % Loop through coils
    for coiln = 1:ncoils
        
        % Print progress message
        if ~isempty(msg_cprog)
            fprintf(repmat('\b',1,length(msg_cprog)));
        end
        msg_cprog = sprintf('Despiking data... (coil %d/%d: %d/%d points flattened)', ...
            coiln,ncoils,totalrmpts,numel(data));
        fprintf(msg_cprog);
        
        % Loop through groups of linked frames
        for ilf = 1:size(linkedframes,1)
            
            % Determine linked frames and get series
            lframes = linkedframes{ilf};
            ts = data(lframes,:,coiln);
            
            % Find points that exceed threshold
            bad = find(abs(ts - mean(ts,1)) > thresh.*std(ts,[],1));
            [r_bad,c_bad] = ind2sub(size(ts),bad);
            uc_bad = unique(c_bad);
            
            % Loop through bad points
            for ibc = 1:length(uc_bad)
                
                % Get bad frames for each bad k-space point
                k_bad = uc_bad(ibc);
                f_bad = lframes(r_bad(c_bad == k_bad));
                data(f_bad,k_bad,coiln) = mean(data(lframes,k_bad,coiln),1);
                
            end
            
            % Update total number of removed points
            totalrmpts = totalrmpts + length(bad);
            
        end
        
    end
    
    % Write data
    fseek(fid,hdr.rdb.off_data,'bof');
    data = reshape([real(data(:))'; imag(data(:))'],1,[]);
    fwrite(fid,data,'short');

    % Close file
    fprintf('\nDone.\n');
    fclose(fid);

end

