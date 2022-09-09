function cfigopen(name)

    % Check if figure with specified name is open
    if isempty(findobj('type','figure','name',name))
        % If not, open it
        figure('name',name);
    else
        % If so, make it the current figure
        figure(findobj('type','figure','name',name))
    end
    
end

