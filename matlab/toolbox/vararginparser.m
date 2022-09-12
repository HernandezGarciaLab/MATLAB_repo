function args = vararginparser(defaults,varargin)

    % Make matlab inputParser object
    p = inputParser;
    
    % Loop through fields and assign to parser
    parmnames = fieldnames(defaults);
    for i = 1:size(parmnames,1)
        parmname = char(parmnames{i});
        p.addParameter(parmname,defaults.(parmname),@(x)1);
    end
    
    % Parse and assign to args
    p.parse(varargin{:});
    args = p.Results;
    
end

