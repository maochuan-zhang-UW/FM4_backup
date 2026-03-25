clc; clear;

%load('/Users/mczhang/Documents/GitHub/FM4/02-data/F_Cl/F_Cl_SKHASH3D.mat');
load('/Users/mczhang/Documents/GitHub/FM4/02-data/E_Po/Felix_combined.mat');
Po_Clu=Felix_combined;clear Felix_combined;
fields = fieldnames(Po_Clu);
for j = 1:numel(fields)
    field = fields{j};
    tokens = regexp(field, '^(Po|NSP|W|DDt|DDSt)_(\d{2})A$', 'tokens');

    if ~isempty(tokens)
        prefix = tokens{1}{1};
        index = tokens{1}{2};
        newField = sprintf('%s_A%02d', prefix, str2double(index));

        % Copy and delete field in all structs
        for i = 1:numel(Po_Clu)
            Po_Clu(i).(newField) = Po_Clu(i).(field);
        end
        Po_Clu = rmfield(Po_Clu, field);
    end
end
