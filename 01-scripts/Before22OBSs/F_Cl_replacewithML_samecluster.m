% Load data
load('/Users/mczhang/Documents/GitHub/FM4/02-data/Before22OBSs/F_Cl/F_Cl_All.mat')   % Po_Clu
load('/Users/mczhang/Documents/GitHub/FM4/02-data/Before22OBSs/E_Po/E_All_PoReplaced.mat') % Felix

% List of polarity fields
po_fields = {'Po_AS1','Po_AS2','Po_CC1','Po_EC1','Po_EC2','Po_EC3','Po_ID1'};

% Build ID index for Felix for fast lookup
felix_ids = [Felix.ID];
felix_map = containers.Map(felix_ids, 1:numel(Felix));

% Loop through Po_Clu and append Felix polarity
for i = 1:numel(Po_Clu)
    id = Po_Clu(i).ID;
    
    % Only process if ID exists in Felix
    if isKey(felix_map, id)
        j = felix_map(id);
        
        for f = 1:numel(po_fields)
            field = po_fields{f};
            
            if isfield(Po_Clu, field) && isfield(Felix, field)
                
                old_val = Po_Clu(i).(field);
                new_val = Felix(j).(field);
                
                % Append only if Felix polarity exists and is not empty
                if ~isempty(new_val)
                    Po_Clu(i).(field) = [old_val, new_val];
                end
                
            end
        end
    end
end
