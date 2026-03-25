% Define the file name
clear;
filename = '/Users/mczhang/Documents/GitHub/SKHASH/SKHASH7/examples/hash3/OUT/out.txt';
% Read the table from the file
data = readtable(filename);
events = table2struct(data);


% % Initialize the structure array
% events = struct('id', [], 'time', [], 'lat', [], 'lon', [], 'depth', [], ...
%                 'avmech', [], 'avfnorm_uncert', [], 'avslip_uncert', [], ...
%                 'polnum', [], 'polmisfit', [], 'mechqual', [], 'mechprob', [], ...
%                 'stdr', [], 'namp', [], 'mavg', [], 'nmult', [], 'max_azimgap', [], ...
%                 'max_takeoff', [], 'avfnorm', [], 'avslip', []);
% 
% % Populate the structure array with data from the table
% for i = 1:height(data)
%     events(i).id = data.event_id(i);
%     %events(i).time = datetime(data.time{i}, 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
%     events(i).lat = data.origin_lat(i);
%     events(i).lon = data.origin_lon(i);
%     events(i).depth = data.origin_depth_km(i);
%     events(i).avmech = [data.strike(i), data.dip(i), data.rake(i)];
%     events(i).avfnorm_uncert = data.fault_plane_uncertainty(i);
%     events(i).avslip_uncert = data.aux_plane_uncertainty(i);
%     events(i).polnum = data.num_p_pol(i);
%     events(i).polmisfit = data.polarity_misfit(i);
%     events(i).mechqual = data.quality{i};
%     events(i).mechprob = data.prob_mech(i);
%     events(i).stdr = data.sta_distribution_ratio(i);
%     events(i).namp = data.num_sp_ratios(i);
%     events(i).mavg = data.sp_misfit(i);
%     events(i).nmult = data.mult_solution_flag(i);
%     % Placeholder values for max_azimgap and max_takeoff, update as needed
%     events(i).max_azimgap = NaN; 
%     events(i).max_takeoff = NaN;
%     events(i).avfnorm = NaN; % Placeholder
%     events(i).avslip = NaN; % Placeholder
% end

event1=events;
for i=1:length(event1)
        % Change coordinate systems
        [event1(i).avfnorm,event1(i).avslip] =...
            fp2fnorm(event1(i).strike,event1(i).dip,event1(i).rake);
        event1(i).b_axis=cross(event1(i).avfnorm,event1(i).avslip);
        event1(i).t_axis=(event1(i).avfnorm+event1(i).avslip)/sqrt(2);
        event1(i).p_axis=(event1(i).avfnorm-event1(i).avslip)/sqrt(2);
        event1(i).p_l=asin(norm(event1(i).p_axis(3))/norm(event1(i).p_axis))/pi*180;
        event1(i).t_l=asin(norm(event1(i).t_axis(3))/norm(event1(i).t_axis))/pi*180;
        event1(i).b_l=asin(norm(event1(i).b_axis(3))/norm(event1(i).b_axis))/pi*180;
        event1(i).check=sin(event1(i).b_l).^2+sin(event1(i).t_l).^2+sin(event1(i).p_l).^2;
        if event1(i).p_l>= 52 && event1(i).t_l<= 35
            event1(i).faultType = 'N';event1(i).color='b';event1(i).color2=[0,0,1];
        elseif event1(i).p_l<=35 && event1(i).t_l>= 52
            event1(i).faultType = 'R';event1(i).color='r';event1(i).color2=[1,0,0];
        elseif event1(i).p_l<= 40 && event1(i).b_l>= 45 && event1(i).t_l<=40
            event1(i).faultType = 'S';event1(i).color='g';event1(i).color2=[0,1,0];
        elseif event1(i).p_l<= 20 && event1(i).t_l>=40 && event1(i).t_l<= 52
            % olique reverse
            event1(i).faultType = 'R';event1(i).color='r';event1(i).color2=[1,0,0];
        elseif event1(i).p_l>= 40 && event1(i).p_l<= 52 && event1(i).t_l<=20
            % Oblique normal
            event1(i).faultType = 'N';event1(i).color='b';event1(i).color2=[0,0,1];
        else
            event1(i).faultType = 'U';event1(i).color='k';event1(i).color2=[0,0,0];
        end
end

 [event1.id] = event1.event_id; event1 = orderfields(event1,[1:0,34,1:33]); event1 = rmfield(event1,'event_id');
[event1.lat] = event1.origin_lat; event1 = orderfields(event1,[1:16,34,17:33]); event1 = rmfield(event1,'origin_lat');
[event1.lon] = event1.origin_lon; event1 = orderfields(event1,[1:17,34,18:33]); event1 = rmfield(event1,'origin_lon');
[event1.depth] = event1.origin_depth_km; event1 = orderfields(event1,[1:18,34,19:33]); event1 = rmfield(event1,'origin_depth_km');
[event1.mechqual] = event1.quality; event1 = orderfields(event1,[1:4,34,5:33]); event1 = rmfield(event1,'quality');
[B,I]=sort([event1.id]);
event1=event1(I);
% Display the field names of the structure array
ind=unique([event1.id]);
n=1;
for i=1:length(ind)
    ind_fm=find([event1.id]==ind(i));
    event(n)=event1(ind_fm(1));
    n=n+1;
end
clear event1;
event1=event;
fieldnames(event1)

% Preallocate new event1D struct
%event1D_new = repmat(event1(1), 1, length(event1));  % clone structure

for i = 1:length(event1)
    % Copy over matching fields
    event1D_new(i).id    = event1(i).id;
    %event1D_new(i).time  = event1(i).time;
    event1D_new(i).time  = datenum(event1(i).time);

    event1D_new(i).lat   = event1(i).lat;
    event1D_new(i).lon   = event1(i).lon;
    event1D_new(i).depth = event1(i).depth;
    
    % avmech = [strike, dip, rake]
    event1D_new(i).avmech = [event1(i).strike, event1(i).dip, event1(i).rake];
    
    % Fields with direct equivalents
    event1D_new(i).mechqual  = event1(i).mechqual;
    event1D_new(i).avfnorm   = event1(i).avfnorm;
    event1D_new(i).avslip    = event1(i).avslip;
    event1D_new(i).b_axis    = event1(i).b_axis;
    event1D_new(i).t_axis    = event1(i).t_axis;
    event1D_new(i).p_axis    = event1(i).p_axis;
    event1D_new(i).p_l       = event1(i).p_l;
    event1D_new(i).t_l       = event1(i).t_l;
    event1D_new(i).b_l       = event1(i).b_l;
    event1D_new(i).check     = event1(i).check;
    event1D_new(i).faultType = event1(i).faultType;
    event1D_new(i).color     = event1(i).color;
    event1D_new(i).color2    = event1(i).color2;
    
    % Additional fields (not in original event1)
    event1D_new(i).avfnorm_uncert  = NaN;
    event1D_new(i).avslip_uncert   = NaN;
    event1D_new(i).polnum          = event1(i).num_p_pol;
    event1D_new(i).polmisfit       = event1(i).polarity_misfit;
    event1D_new(i).mechprob        = event1(i).prob_mech;
    event1D_new(i).stdr            = event1(i).sta_distribution_ratio;
    event1D_new(i).namp            = event1(i).num_sp_ratios;
    event1D_new(i).mavg            = NaN;
    event1D_new(i).nmult           = event1(i).mult_solution_flag;
    event1D_new(i).max_azimgap     = NaN;
    event1D_new(i).max_takeoff     = NaN;
    event1D_new(i).Mw              = event1(i).magnitude;
    event1D_new(i).Mw2             = NaN;
    event1D_new(i).Mo              = NaN;
end

% Optional: overwrite event1D
event1 = event1D_new;


save('/Users/mczhang/Documents/GitHub/FM3/02-data/G_FM/GA_FM_SKHASH_all_24060_10percent.mat', 'event1');

