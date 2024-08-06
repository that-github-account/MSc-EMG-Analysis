%% Settings

clear all
close all
clc

% Choose which datasets to graph

plot_force = true; %either both or one must be true
plot_emg = true;

plot_right_hand = true; % either both or one must be true 
plot_left_hand = true;

plot_fail = true;
plot_success = true;

plot_passive = false;

all_subjects = true; % specify subject if false 
specific_subject = 1;

which_data = 2; %1 = trials excluded, 2 = trials not excluded

%% Set-Up Variables

if all_subjects == true
    n_subjects = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","19"];
else
    n_subjects = string(specific_subject);
end

namespace_hands = "";
if plot_right_hand == true
    namespace_hands = "con_R_";
end
if plot_left_hand == true
    if namespace_hands == ""
        namespace_hands = "con_L_";
    else
        namespace_hands = [namespace_hands, "con_L_"];
    end
end
if namespace_hands == ""
    disp("Please Specify A Hand!")
    return
end

namespace_type = "";
if plot_fail == true
    for i = 1:length(namespace_hands)
        var = (namespace_hands(i) + "fail_");
        if namespace_type == ""
            namespace_type = var;
        else
            namespace_type = [namespace_type, var];
        end
    end
end

if plot_success == true
    for i = 1:length(namespace_hands)
        var = (namespace_hands(i) + "succ_");
        if namespace_type == ""
            namespace_type = var;
        else
            namespace_type = [namespace_type, var];
        end
    end
end

namespace_force = ["force", "forcevel"];
namespace_measure = "";
if plot_force == true
    for i = 1:length(namespace_type)
        for o = 1:length(namespace_force)
            var = (namespace_type(i) + namespace_force(o));
            if namespace_measure == ""
                namespace_measure = var;
            else
                namespace_measure = [namespace_measure, var];
            end
        end
    end
end

namespace_emg = ["emg", "emgvel"];
if plot_emg == true
    for i = 1:length(namespace_type)
        for o = 1:length(namespace_emg)
            var = (namespace_type(i) + namespace_emg(o));
            if namespace_measure == ""
                namespace_measure = var;
            else
                namespace_measure = [namespace_measure, var];
            end
        end
    end
end

namespace_complete = namespace_measure;
if plot_passive == true
    for i = 1:length(namespace_measure)
        var = namespace_measure(i) + "_pass";
        namespace_complete = [namespace_complete, var];
    end
end

clear namespace_hands
clear namespace_type
clear namespace_measure
clear namespace_force
clear namespace_emg

%% Read in all Data

if which_data == 1
    f_text = 'Analog_all_S';
elseif which_data == 2
    f_text = 'Analog_all_noTriExc_S';
end

filenames = "";
if length(n_subjects) > 1
    for i = 1:length(n_subjects)
        if filenames == ""
            filenames = f_text + n_subjects(i);
        else
            filenames = [filenames, f_text + n_subjects(i)];
        end
    end
else
    filenames = f_text + n_subjects;
end

%% PLOTS

all_force_data = [];
all_emg_data = [];

for p = 1:length(filenames)
    filename = filenames(p);

    load(filename);

    %time variables
    time_resplock = -500:1000;%resp-locked data starts 500ms before response till 1000ms after
    time_stoplock = -1000:1000;%stop-locked data starts 1000ms before stop signal till 1000ms after stop signal

    clear spike_stop_cond_base_mn
    clear spike_stop_cond_base
    clear spike_stop_cond_mn
    clear spike_resp_cond_mn

    spike_resp_cond = spike_resplock_conds;
    clear spike_resplock_conds

    %all data is stop-locked

    %Active hand
    con_R_succ_force = cell2mat(spike_stop_cond(3,5));
    con_R_succ_forcevel = cell2mat(spike_stop_cond(3,6));
    con_R_succ_emg = cell2mat(spike_stop_cond(3,10));
    con_R_succ_emgvel = cell2mat(spike_stop_cond(3,11));
    con_R_fail_force = cell2mat(spike_stop_cond(4,5));
    con_R_fail_forcevel = cell2mat(spike_stop_cond(4,6));
    con_R_fail_emg = cell2mat(spike_stop_cond(4,10));
    con_R_fail_emgvel = cell2mat(spike_stop_cond(4,11));
    con_L_succ_force = cell2mat(spike_stop_cond(5,13));
    con_L_succ_forcevel = cell2mat(spike_stop_cond(5,14));
    con_L_succ_emg = cell2mat(spike_stop_cond(5,18));
    con_L_succ_emgvel = cell2mat(spike_stop_cond(5,19));
    con_L_fail_force = cell2mat(spike_stop_cond(6,13));
    con_L_fail_forcevel = cell2mat(spike_stop_cond(6,14));
    con_L_fail_emg = cell2mat(spike_stop_cond(6,18));
    con_L_fail_emgvel = cell2mat(spike_stop_cond(6,19));

    %passive hand
    con_R_succ_force_pass = cell2mat(spike_stop_cond(3,13));
    con_R_succ_forcevel_pass = cell2mat(spike_stop_cond(3,14));
    con_R_succ_emg_pass = cell2mat(spike_stop_cond(3,18));
    con_R_succ_emgvel_pass = cell2mat(spike_stop_cond(3,19));
    con_R_fail_force_pass = cell2mat(spike_stop_cond(4,13));
    con_R_fail_forcevel_pass = cell2mat(spike_stop_cond(4,14));
    con_R_fail_emg_pass = cell2mat(spike_stop_cond(4,18));
    con_R_fail_emgvel_pass = cell2mat(spike_stop_cond(4,19));
    con_L_succ_force_pass = cell2mat(spike_stop_cond(5,5));
    con_L_succ_forcevel_pass = cell2mat(spike_stop_cond(5,6));
    con_L_succ_emg_pass = cell2mat(spike_stop_cond(5,10));
    con_L_succ_emgvel_pass = cell2mat(spike_stop_cond(5,11));
    con_L_fail_force_pass = cell2mat(spike_stop_cond(6,5));
    con_L_fail_forcevel_pass = cell2mat(spike_stop_cond(6,6));
    con_L_fail_emg_pass = cell2mat(spike_stop_cond(6,10));
    con_L_fail_emgvel_pass = cell2mat(spike_stop_cond(6,11));
    
    linewidth = 1.5;

    figure_placeholder = p*100;
    
%% FORCE PLOTS

    if plot_force == true

        n_figures = namespace_complete(~cellfun(@isempty, strfind(namespace_complete, "force")));

        for k = 1:2:length(n_figures)

            dataset = eval(n_figures(k));
            setvel = eval(n_figures(k+1));

            rmcon = erase(n_figures(k), "con");
            rm_ = erase(rmcon, "_");
            figure_name = char("P" + n_subjects(p) + ": " + rm_);

            for i = 1:size(dataset, 1)
                n_trial = i;
                stopping_amplitude = max(dataset(n_trial, 1001:2001));
                stopping_amplitude_time = find(dataset(n_trial,:) == stopping_amplitude);

                stopping_slope = max(setvel(n_trial, :));

                stopping_end = min(abs(dataset(n_trial, stopping_amplitude_time:2001)));
                stopping_end_time = find(abs(dataset(n_trial,:)) == stopping_end);
                stopping_start_time = 1001;

                avg_volatility_range = 801;
                avg_volatility_range = 1:avg_volatility_range;
                avg_volatility = mean(abs(dataset(n_trial, avg_volatility_range)));
                avg_volatilityvel = mean(abs(setvel(n_trial, avg_volatility_range)));

                time_adjust = -100;
                set_adj = dataset(n_trial, 1:stopping_amplitude_time+time_adjust) - avg_volatility;
                setvel_adj = setvel(n_trial, 1:stopping_amplitude_time+time_adjust) - avg_volatilityvel;

                loc_0_vel = find(abs(setvel_adj(:)) < 0.05);
                if isempty(loc_0_vel)
                    vel_0_closest_to_peak = loc_0_vel;
                else
                    vel_0_closest_to_peak = loc_0_vel(length(loc_0_vel));
                end

                if isempty(vel_0_closest_to_peak)
                    disp("skipped");
                else
                    stopping_start_time = vel_0_closest_to_peak;
                end

                stopping_duration = stopping_end_time - stopping_start_time; %gets recalculated again below in case values change
                if stopping_duration == 0
                    disp("no surface")
                    stopping_surface = 0;
                else
                    stopping_surface = trapz(stopping_start_time:stopping_end_time, dataset(i, stopping_start_time:stopping_end_time));
                end

                participant_id = p;

                if participant_id == 18
                    participant_id = 19;
                end

                go_direction = char(rm_);
                go_direction = go_direction(1);

                if contains(n_figures(k), "succ") == true
                    classification = "success";
                else
                    classification = "fail";
                end

                all_force_data = [all_force_data; stopping_start_time-1001, stopping_amplitude_time-1001, stopping_end_time-1001, stopping_amplitude, stopping_surface, stopping_duration, stopping_slope, avg_volatility, avg_volatilityvel, participant_id, go_direction, classification];
            end
        end
    end

%% EMG PLOTS

    if plot_emg == true

        n_figures = namespace_complete(~cellfun(@isempty, strfind(namespace_complete, "emg")));

        for k = 1:2:length(n_figures)

            dataset = eval(n_figures(k));
            setvel = eval(n_figures(k+1));

            rmcon = erase(n_figures(k), "con");
            rm_ = erase(rmcon, "_");
            figure_name = char("P" + n_subjects(p) + ": " + rm_);

            for i = 1:size(dataset, 1)
                n_trial = i;

                stopping_amplitude = max(dataset(n_trial, 1001:2001));
                stopping_amplitude_time = find(dataset(n_trial,:) == stopping_amplitude);

                stopping_slope = max(setvel(n_trial, :));

                avg_volatility_range = 801;
                avg_volatility_range = 1:avg_volatility_range;
                avg_volatility = mean(abs(dataset(n_trial, avg_volatility_range)));
                avg_volatilityvel = mean(abs(setvel(n_trial, avg_volatility_range)));

                set_adj = dataset(n_trial,:) - avg_volatility;
                intersection_points = find(abs(set_adj(:)) < 0.1);
                intersection_points_past_peak = intersection_points(intersection_points > stopping_amplitude_time);

                if isempty(intersection_points_past_peak)
                    stopping_end = min(abs(dataset(n_trial, stopping_amplitude_time:2001)));
                    stopping_end_time = find(abs(dataset(n_trial,:)) == stopping_end);
                else
                    stopping_end_time = intersection_points_past_peak(1);
                end

                stopping_start_time = 1001;

                time_adjust = -100;
                set_adj = dataset(n_trial, 1:stopping_amplitude_time+time_adjust) - avg_volatility;
                setvel_adj = setvel(n_trial, 1:stopping_amplitude_time+time_adjust) - avg_volatilityvel;

                loc_0_vel = find(abs(setvel_adj(:)) < 0.1);
                if isempty(loc_0_vel)
                    vel_0_closest_to_peak = loc_0_vel;
                else
                    vel_0_closest_to_peak = loc_0_vel(length(loc_0_vel));
                end

                if isempty(vel_0_closest_to_peak)
                    disp("skipped");
                else
                    stopping_start_time = vel_0_closest_to_peak;
                end

                stopping_duration = stopping_end_time - stopping_start_time; %gets recalculated again below in case values change
                if stopping_duration == 0
                    disp("no surface")
                    stopping_surface = 0;
                else
                    stopping_surface = trapz(stopping_start_time:stopping_end_time, dataset(i, stopping_start_time:stopping_end_time));
                end

                participant_id = p;

                if participant_id == 18
                    participant_id = 19;
                end

                go_direction = char(rm_);
                go_direction = go_direction(1);

                if contains(n_figures(k), "succ") == true
                    classification = "success";
                else
                    classification = "fail";
                end

                all_emg_data = [all_emg_data; stopping_start_time-1001, stopping_amplitude_time-1001, stopping_end_time-1001, stopping_amplitude ,stopping_surface, stopping_duration, stopping_slope, avg_volatility, avg_volatilityvel, participant_id, go_direction, classification];                
            end
        end
    end

save("force_data.mat", "all_force_data");
save("emg_data.mat", "all_emg_data");

writematrix(all_force_data, "force_data.csv");
writematrix(all_emg_data, "emg_data.csv");

end

