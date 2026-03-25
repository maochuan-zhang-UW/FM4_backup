function event1 = readSKHASH3D_skhash2(saveFilename)
    % READSKHASH3D Reads and processes SKHASH earthquake data from a CSV file, saving to a .mat file.
    %
    % Input:
    %   saveFilename - String specifying the path and name for the output .mat file
    %                  (e.g., 'GA_FM_SKHASH_3D_test.mat').
    %
    % Output:
    %   event1 - Structured array containing processed earthquake data with fields such as
    %            id, time, lat, lon, depth, faultType, etc.

    % Hardcoded input file path
    filename = '/Users/mczhang/Documents/GitHub/SKHASH/SKHASH2/examples/smile/OUT/out.csv';

    % Read the table from the file
    data = readtable(filename);

    % Initialize the structure array
    events = struct('id', [], 'time', [], 'lat', [], 'lon', [], 'depth', [], ...
                    'avmech', [], 'avfnorm_uncert', [], 'avslip_uncert', [], ...
                    'polnum', [], 'polmisfit', [], 'mechqual', [], 'mechprob', [], ...
                    'stdr', [], 'namp', [], 'mavg', [], 'nmult', [], 'max_azimgap', [], ...
                    'max_takeoff', [], 'avfnorm', [], 'avslip', []);

    % Populate the structure array with data from the table
    for i = 1:height(data)
        events(i).id = data.event_id(i);
        events(i).time = datenum(data.time(i));
        events(i).lat = data.origin_lat(i);
        events(i).lon = data.origin_lon(i);
        events(i).depth = data.origin_depth_km(i);
        events(i).avmech = [data.strike(i), data.dip(i), data.rake(i)];
        events(i).avfnorm_uncert = data.fault_plane_uncertainty(i);
        events(i).avslip_uncert = data.aux_plane_uncertainty(i);
        events(i).polnum = data.num_p_pol(i);
        events(i).polmisfit = data.polarity_misfit(i);
        events(i).mechqual = data.quality{i};
        events(i).mechprob = data.prob_mech(i);
        events(i).stdr = data.sta_distribution_ratio(i);
        events(i).namp = data.num_sp_ratios(i);
        events(i).mavg = data.sp_misfit(i);
        events(i).nmult = data.mult_solution_flag(i);
        events(i).max_azimgap = NaN;
        events(i).max_takeoff = NaN;
        events(i).avfnorm = NaN;
        events(i).avslip = NaN;
    end

    % Process each event to compute fault parameters and classify fault type
    event1 = events;
    for i = 1:length(event1)
        % Change coordinate systems
        [event1(i).avfnorm, event1(i).avslip] = ...
            fp2fnorm(event1(i).avmech(1), event1(i).avmech(2), event1(i).avmech(3));
        event1(i).b_axis = cross(event1(i).avfnorm, event1(i).avslip);
        event1(i).t_axis = (event1(i).avfnorm + event1(i).avslip) / sqrt(2);
        event1(i).p_axis = (event1(i).avfnorm - event1(i).avslip) / sqrt(2);
        event1(i).p_l = asin(abs(event1(i).p_axis(3)) / norm(event1(i).p_axis)) / pi * 180;
        event1(i).t_l = asin(abs(event1(i).t_axis(3)) / norm(event1(i).t_axis)) / pi * 180;
        event1(i).b_l = asin(abs(event1(i).b_axis(3)) / norm(event1(i).b_axis)) / pi * 180;
        event1(i).check = sin(event1(i).b_l)^2 + sin(event1(i).t_l)^2 + sin(event1(i).p_l)^2;
        if event1(i).p_l >= 52 && event1(i).t_l <= 35
            event1(i).faultType = 'N';
            event1(i).color = 'b';
            event1(i).color2 = [0, 0, 1];
        elseif event1(i).p_l <= 35 && event1(i).t_l >= 52
            event1(i).faultType = 'R';
            event1(i).color = 'r';
            event1(i).color2 = [1, 0, 0];
        elseif event1(i).p_l <= 40 && event1(i).b_l >= 45 && event1(i).t_l <= 40
            event1(i).faultType = 'S';
            event1(i).color = 'g';
            event1(i).color2 = [0, 1, 0];
        elseif event1(i).p_l <= 20 && event1(i).t_l >= 40 && event1(i).t_l <= 52
            event1(i).faultType = 'R';
            event1(i).color = 'r';
            event1(i).color2 = [1, 0, 0];
        elseif event1(i).p_l >= 40 && event1(i).p_l <= 52 && event1(i).t_l <= 20
            event1(i).faultType = 'N';
            event1(i).color = 'b';
            event1(i).color2 = [0, 0, 1];
        else
            event1(i).faultType = 'U';
            event1(i).color = 'k';
            event1(i).color2 = [0, 0, 0];
        end
    end

    % Sort and select unique events by id
    [~, I] = sort([event1.id]);
    event1 = event1(I);
    ind = unique([event1.id]);
    event = event1(1);
    n = 1;
    for i = 1:length(ind)
        ind_fm = find([event1.id] == ind(i));
        event(n) = event1(ind_fm(1));
        n = n + 1;
    end
    event1 = event;

    % Save the processed data to the specified file
    save(saveFilename, 'event1');
end