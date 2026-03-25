clear;close all;
path='/Users/mczhang/Documents/GitHub/FM4/02-data/';
path_old='/Users/mczhang/Documents/GitHub/FM3/02-data/';
fields={'AS1','AS2','CC1','EC1','EC2','EC3','ID1'};
groups = {'E1','E2','E3','E4','S1','W1','W2'};
nz=1;
resultman=[];
for gp=1:length(groups)
    load([path,'/E_Po/E_',groups{gp},'.mat']);
    resultman=[resultman,Felix];
end
clear Felix;
Felix=resultman;
% Example vector
% Get unique values
% Find unique IDs and their indices, keeping the first occurrence
[~, ia, ~] = unique([Felix.ID], 'first');

% Keep only the first occurrence of each ID
Felix = Felix(ia);

    for i = 1:length(Felix)
        for kz = 1:length(fields)
            if ~isfield(Felix(i), ['SP_', fields{kz}])
                Felix(i).(['SP_', fields{kz}]) = 0;
            end
        end
    end
    for i=1:length(Felix)
        for kz=1:length(fields)
            if ~isempty(Felix(i).(['NSP_', fields{kz}]))
                Felix(i).(['SP_', fields{kz}])=log10(Felix(i).(['NSP_', fields{kz}])(2)/Felix(i).(['NSP_', fields{kz}])(3));
            else
                Felix(i).(['SP_', fields{kz}])=1;
            end
        end
    end
Felix = rmfield(Felix, 'W_AS1');
Felix = rmfield(Felix, 'W_AS2');
Felix = rmfield(Felix, 'W_CC1');
Felix = rmfield(Felix, 'W_EC1');
Felix = rmfield(Felix, 'W_EC2');
Felix = rmfield(Felix, 'W_EC3');
Felix = rmfield(Felix, 'W_ID1');

   for i=1:length(Felix)
        Felix(i).SP_All=sum([abs(sign(Felix(i).SP_AS1)),abs(sign(Felix(i).SP_AS2)),abs(sign(Felix(i).SP_CC1)),abs(sign(Felix(i).SP_EC3)),abs(sign(Felix(i).SP_EC2)),abs(sign(Felix(i).SP_EC1)),abs(sign(Felix(i).SP_ID1))],'omitnan');
   end


save(['/Users/mczhang/Documents/GitHub/FM4/02-data/E_Po/E_All.mat'], 'Felix');
% % Verify duplicates are removed
% [uniqueA, ~, ~] = unique([Felix.ID]);
% counts = hist([Felix.ID], uniqueA);
% duplicates = uniqueA(counts > 1);
% % Display results
% if isempty(duplicates)
%     fprintf('No duplicated numbers found.\n');
% else
%     fprintf('Duplicated numbers: %s\n', mat2str(duplicates));
% end
% 
% 
% % Optionally save the deduplicated Felix structure
% save([path, '/A_ID/Felix_deduplicated.mat'], 'Felix');