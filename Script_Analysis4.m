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
    
%% Extract Force Data

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

                end_cutoff = 400;

                stopping_amplitude = max(dataset(n_trial, 1001:2001 - end_cutoff));
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
                    start_value_replace = find(dataset(n_trial, vel_0_closest_to_peak:2001)>0, 1) + vel_0_closest_to_peak-1;
                    if isempty(start_value_replace)
                        disp("no more positives")
                    else
                        stopping_start_time = start_value_replace;
                    end
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

%% Extract EMG Data

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

                end_cutoff = 400;

                stopping_amplitude = max(dataset(n_trial, 1001:2001 - end_cutoff));
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
                    start_value_replace = find(dataset(n_trial, vel_0_closest_to_peak:2001)>0, 1) + vel_0_closest_to_peak-1;
                    if isempty(start_value_replace)
                        disp("no more positives")
                    else
                        stopping_start_time = start_value_replace;
                    end
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
end

%% Correct Force/EMG Data

% for trials where extracted max amplitude over the response window was
% negative there was no response so ampltiude, duration, surface, and slope
% should be 0 

% for trials where extracted surface is negative (can happen due to chances) the area classified is
% predominantly negative and so does not actually reflect and increase
% despite the amplitude still being large (on these trials it will be
% comparateively very small) so theses are also changed to 0

index_of_values_to_replace = str2double(all_force_data(:, 4:5)) < 0;
index_of_values_to_replace = (index_of_values_to_replace(:, 1) + index_of_values_to_replace(:, 2)) > 0;
all_force_data(index_of_values_to_replace, 4:7) = 0;

index_of_values_to_replace = str2double(all_emg_data(:, 4:5)) < 0;
index_of_values_to_replace = (index_of_values_to_replace(:, 1) + index_of_values_to_replace(:, 2)) > 0;
all_emg_data(index_of_values_to_replace, 4:7) = 0;

%% Settings

quantile_size = 0.35;

%

participant_numbers = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","19"];

participant_velocities = [];
participant_velocities_LR = [];

% select success classified stop trials

reclassification_data = all_force_data(all_force_data(:, 12) == "success",:);
reclassification_data = [reclassification_data, str2double(reclassification_data(:, 7)) - str2double(reclassification_data(:, 9))];

% find median/quantile for each hand of each participant 

for i = participant_numbers

    participant_number = str2num(i);

    if participant_number == 19
        participant_number = 18;
    end


    participant_data = reclassification_data(reclassification_data(:,10) == participant_numbers(participant_number), :);
    
    participant_data_R = participant_data(participant_data(:, 11) == "R", :);
    participant_data_L = participant_data(participant_data(:, 11) == "L", :);
    
    velocity_R = str2double(participant_data_R(:, 13));
    velocity_L = str2double(participant_data_L(:, 13));

    median_velocity_R = median(velocity_R);
    median_velocity_L = median(velocity_L);

    lower_quantile_R = quantile(velocity_R, quantile_size);
    lower_quantile_L = quantile(velocity_L, quantile_size);

    participant_velocities_LR = [participant_velocities_LR; i, "R", median_velocity_R, lower_quantile_R; i, "L", median_velocity_L, lower_quantile_L];

end

% merge median/quantile velocity into all data, perform a comparison to
% to reclassify, extract those values

reclassification_identifier = append(reclassification_data(:, 10), reclassification_data(:, 11));
reclassification_data = [reclassification_data, reclassification_identifier];

velocities_identifier = append(participant_velocities_LR(:, 1), participant_velocities_LR(:, 2));
participant_velocities_LR = [participant_velocities_LR, velocities_identifier];

[~, ii] = ismember(reclassification_data(:,14), participant_velocities_LR(:,5));

joined_data_LR = [reclassification_data, participant_velocities_LR(ii,3:4)];

% 15 for median, 16 for specified quantile 

include_LR = str2double(joined_data_LR(:, 13)) > str2double(joined_data_LR(:, 16));

reclassification_data_joined_LR = [reclassification_data, include_LR];

columns_to_include_LR = reclassification_data_joined_LR(:, [10,11,15]);

%% Graph group average and hands for each participant


all_t_rec_fail_force = [];
all_t_rec_succ_force = [];

all_t_rec_fail_emg = [];
all_t_rec_succ_emg = [];

for p = 1:length(filenames)

    p_str = string(p);

    filename = filenames(p);

    load(filename);

    clear spike_stop_cond_base_mn
    clear spike_stop_cond_base
    clear spike_stop_cond_mn
    clear spike_resp_cond_mn

    spike_resp_cond = spike_resplock_conds;
    clear spike_resplock_conds

    %all data is stop-locked

    %Active hand
    R_succ_force = cell2mat(spike_stop_cond(3,5));
    R_succ_forcevel = cell2mat(spike_stop_cond(3,6));
    R_succ_emg = cell2mat(spike_stop_cond(3,10));
    R_succ_emgvel = cell2mat(spike_stop_cond(3,11));

    L_succ_force = cell2mat(spike_stop_cond(5,13));
    L_succ_forcevel = cell2mat(spike_stop_cond(5,14));
    L_succ_emg = cell2mat(spike_stop_cond(5,18));
    L_succ_emgvel = cell2mat(spike_stop_cond(5,19));

    clear spike_resp_cond spike_stop_cond

    if p_str == "18"
        p_str = "19";
    end

    % EMG

    participant_index = columns_to_include_LR(columns_to_include_LR(:, 1) == p_str, :);

    rh_index = strcmp(participant_index(participant_index(:, 2) == "R", :), "true");
    lh_index = strcmp(participant_index(participant_index(:, 2) == "L", :), "true");

    R_succ_emg_reclassified_fail = R_succ_emg(rh_index(:, 3), :);
    R_succ_emg_reclassified_succ = R_succ_emg(~rh_index(:, 3), :);

    L_succ_emg_reclassified_fail = L_succ_emg(lh_index(:, 3), :);
    L_succ_emg_reclassified_succ = L_succ_emg(~lh_index(:, 3), :);

    R_succ_emg_reclassified_fail_grp_mean = mean(R_succ_emg_reclassified_fail);
    R_succ_emg_reclassified_succ_grp_mean = mean(R_succ_emg_reclassified_succ);

    L_succ_emg_reclassified_fail_grp_mean = mean(L_succ_emg_reclassified_fail);
    L_succ_emg_reclassified_succ_grp_mean = mean(L_succ_emg_reclassified_succ);

    participant_index = strcmp(participant_index(:, 3), "true");

    succ_emg = [R_succ_emg; L_succ_emg];
    succ_emg_reclassified_fail = succ_emg(participant_index, :);
    succ_emg_reclassified_succ = succ_emg(~participant_index, :);

    all_t_rec_fail_emg = [all_t_rec_fail_emg; succ_emg_reclassified_fail];
    all_t_rec_succ_emg = [all_t_rec_succ_emg; succ_emg_reclassified_succ];

    succ_emg_reclassified_fail = mean(succ_emg_reclassified_fail);
    succ_emg_reclassified_succ = mean(succ_emg_reclassified_succ);

    figure(p);
    subplot(2, 3, 1);
    plot(time_stoplock, succ_emg_reclassified_fail, 'Color', 'r', 'LineWidth', linewidth); hold on;
    plot(time_stoplock, succ_emg_reclassified_succ, 'Color', 'b', 'LineWidth', linewidth);
    ylim([-5 8])
    title(append(p_str, " avg EMG"))

    subplot(2, 3, 2);
    plot(time_stoplock, R_succ_emg_reclassified_fail_grp_mean, 'Color', 'r', 'LineWidth', linewidth); hold on;
    plot(time_stoplock, R_succ_emg_reclassified_succ_grp_mean, 'Color', 'b', 'LineWidth', linewidth);
    ylim([-5 8])
    title(append(p_str, " R EMG"))

    subplot(2, 3, 3);
    plot(time_stoplock, L_succ_emg_reclassified_fail_grp_mean, 'Color', 'r', 'LineWidth', linewidth); hold on;
    plot(time_stoplock, L_succ_emg_reclassified_succ_grp_mean, 'Color', 'b', 'LineWidth', linewidth);
    ylim([-5 8])
    title(append(p_str, " L EMG"))
    legend('Fail', 'Success')

    % Force

    participant_index = columns_to_include_LR(columns_to_include_LR(:, 1) == p_str, :);

    rh_index = strcmp(participant_index(participant_index(:, 2) == "R", :), "true");
    lh_index = strcmp(participant_index(participant_index(:, 2) == "L", :), "true");

    R_succ_force_reclassified_fail = R_succ_force(rh_index(:, 3), :);
    R_succ_force_reclassified_succ = R_succ_force(~rh_index(:, 3), :);

    L_succ_force_reclassified_fail = L_succ_force(lh_index(:, 3), :);
    L_succ_force_reclassified_succ = L_succ_force(~lh_index(:, 3), :);

    R_succ_force_reclassified_fail_grp_mean = mean(R_succ_force_reclassified_fail);
    R_succ_force_reclassified_succ_grp_mean = mean(R_succ_force_reclassified_succ);

    L_succ_force_reclassified_fail_grp_mean = mean(L_succ_force_reclassified_fail);
    L_succ_force_reclassified_succ_grp_mean = mean(L_succ_force_reclassified_succ);

    participant_index = strcmp(participant_index(:, 3), "true");

    succ_force = [R_succ_force; L_succ_force];
    succ_force_reclassified_fail = succ_force(participant_index, :);
    succ_force_reclassified_succ = succ_force(~participant_index, :);

    all_t_rec_fail_force = [all_t_rec_fail_force; succ_force_reclassified_fail];
    all_t_rec_succ_force = [all_t_rec_succ_force; succ_force_reclassified_succ];

    succ_force_reclassified_fail = mean(succ_force_reclassified_fail);
    succ_force_reclassified_succ = mean(succ_force_reclassified_succ);

    figure(p);
    subplot(2, 3, 4);
    plot(time_stoplock, succ_force_reclassified_fail, 'Color', 'r', 'LineWidth', linewidth); hold on;
    plot(time_stoplock, succ_force_reclassified_succ, 'Color', 'b', 'LineWidth', linewidth);
    ylim([-3 3])
    title(append(p_str, " avg Force"))

    subplot(2, 3, 5);
    plot(time_stoplock, R_succ_force_reclassified_fail_grp_mean, 'Color', 'r', 'LineWidth', linewidth); hold on;
    plot(time_stoplock, R_succ_force_reclassified_succ_grp_mean, 'Color', 'b', 'LineWidth', linewidth);
    ylim([-3 3])
    title(append(p_str, " R Force"))

    subplot(2, 3, 6);
    plot(time_stoplock, L_succ_force_reclassified_fail_grp_mean, 'Color', 'r', 'LineWidth', linewidth); hold on;
    plot(time_stoplock, L_succ_force_reclassified_succ_grp_mean, 'Color', 'b', 'LineWidth', linewidth);
    ylim([-3 3])
    title(append(p_str, " L Force"))

end

all_t_rec_fail_force = mean(all_t_rec_fail_force);
all_t_rec_succ_force = mean(all_t_rec_succ_force);

all_t_rec_fail_emg = mean(all_t_rec_fail_emg);
all_t_rec_succ_emg = mean(all_t_rec_succ_emg);

figure()
subplot(2,1,1)
plot(time_stoplock, all_t_rec_fail_emg, 'Color', 'r', 'LineWidth', linewidth); hold on;
plot(time_stoplock, all_t_rec_succ_emg, 'Color', 'b', 'LineWidth', linewidth);
legend('Fail EMG', 'Succ. EMG');
title('Grand Average EMG');

subplot(2,1,2)
plot(time_stoplock, all_t_rec_fail_force, 'Color', 'r', 'LineWidth', linewidth); hold on;
plot(time_stoplock, all_t_rec_succ_force, 'Color', 'b', 'LineWidth', linewidth);
legend('Fail Force', 'Succ. Force');
title('Grand Average Force');


%% merging all into final data

%adding identifiers to force

all_force_data = [all_force_data, append(all_force_data(:, 10), all_force_data(:, 11), all_force_data(:, 12))];

unique_names = unique(all_force_data(:, 13), 'stable');
count_of_names = cellfun(@(x) sum(ismember(all_force_data(:, 13),x)), unique_names);
count = [unique_names, count_of_names];

trial_numbers = [];
for i = 1:length(count)
    t = [1:str2double(count(i, 2))]';
    trial_numbers = [trial_numbers; t];
end

all_force_data = [all_force_data, trial_numbers];

all_force_data = [all_force_data, append(all_force_data(:, 13), all_force_data(:, 14))];


% adding identifiers to emg

all_emg_data = [all_emg_data, append(all_emg_data(:, 10), all_emg_data(:, 11), all_emg_data(:, 12))];

unique_names = unique(all_emg_data(:, 13), 'stable');
count_of_names = cellfun(@(x) sum(ismember(all_emg_data(:, 13),x)), unique_names);
count = [unique_names, count_of_names];

trial_numbers = [];
for i = 1:length(count)
    t = [1:str2double(count(i, 2))]';
    trial_numbers = [trial_numbers; t];
end

all_emg_data = [all_emg_data, trial_numbers];

all_emg_data = [all_emg_data, append(all_emg_data(:, 13), all_emg_data(:, 14))];


% adding identifiers to reclassification index

columns_to_include_LR = [columns_to_include_LR, append(columns_to_include_LR(:, 1), columns_to_include_LR(:, 2), "success")];

unique_names = unique(columns_to_include_LR(:, 4), 'stable');
count_of_names = cellfun(@(x) sum(ismember(columns_to_include_LR(:, 4),x)), unique_names);
count = [unique_names, count_of_names];

trial_numbers = [];
for i = 1:length(count)
    t = [1:str2double(count(i, 2))]';
    trial_numbers = [trial_numbers; t];
end

columns_to_include_LR = [columns_to_include_LR, trial_numbers];

columns_to_include_LR = [columns_to_include_LR, append(columns_to_include_LR(:, 4), columns_to_include_LR(:, 5))];


% merging data 

[~, ii] = ismember(all_force_data(:, 15), columns_to_include_LR(:, 6));

columns_to_include_LR = [columns_to_include_LR; NaN(1, 6)];

ii = changem(ii, length(columns_to_include_LR));

all_force_data = [all_force_data, columns_to_include_LR(ii, 3)];
all_emg_data = [all_emg_data, columns_to_include_LR(ii, 3)];

all_force_data(:, 13) = [];
all_force_data(:, 14) = [];

all_emg_data(:, 13) = [];
all_emg_data(:, 14) = [];

writematrix(all_force_data, "force_data_reclassified.csv");
writematrix(all_emg_data, "emg_data_reclassified.csv");
