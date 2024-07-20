%% Setup
%script to load individual subject data and plot resp and stop-locked data
clear all
close all
clc

sub = 1;%subject number 
which_data = 2;%1 = trials excluded, 2 = trials not excluded

if which_data == 1
    f_text = 'Analog_all_S';
elseif which_data == 2
    f_text = 'Analog_all_noTriExc_S';
end

%Load data
subnum_text = num2str(sub);
filename = [f_text subnum_text];
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

%% Reading All Data

%stop-locked analysis
%Active hand
con_R_succ_force = cell2mat(spike_stop_cond(3,5));
con_R_succ_vel = cell2mat(spike_stop_cond(3,6));
con_R_succ_emg = cell2mat(spike_stop_cond(3,10));
con_R_succ_emgvel = cell2mat(spike_stop_cond(3,11));
con_R_fail_force = cell2mat(spike_stop_cond(4,5));
con_R_fail_vel = cell2mat(spike_stop_cond(4,6));
con_R_fail_emg = cell2mat(spike_stop_cond(4,10));
con_R_fail_emgvel = cell2mat(spike_stop_cond(4,11));
con_L_succ_force = cell2mat(spike_stop_cond(5,13));
con_L_succ_vel = cell2mat(spike_stop_cond(5,14));
con_L_succ_emg = cell2mat(spike_stop_cond(5,18));
con_L_succ_emgvel = cell2mat(spike_stop_cond(5,19));
con_L_fail_force = cell2mat(spike_stop_cond(6,13));
con_L_fail_vel = cell2mat(spike_stop_cond(6,14));
con_L_fail_emg = cell2mat(spike_stop_cond(6,18));
con_L_fail_emgvel = cell2mat(spike_stop_cond(6,19));

%passive hand
con_R_succ_force_pass = cell2mat(spike_stop_cond(3,13));
con_R_succ_vel_pass = cell2mat(spike_stop_cond(3,14));
con_R_succ_emg_pass = cell2mat(spike_stop_cond(3,18));
con_R_succ_emgvel_pass = cell2mat(spike_stop_cond(3,19));
con_R_fail_force_pass = cell2mat(spike_stop_cond(4,13));
con_R_fail_vel_pass = cell2mat(spike_stop_cond(4,14));
con_R_fail_emg_pass = cell2mat(spike_stop_cond(4,18));
con_R_fail_emgvel_pass = cell2mat(spike_stop_cond(4,19));
con_L_succ_force_pass = cell2mat(spike_stop_cond(5,5));
con_L_succ_vel_pass = cell2mat(spike_stop_cond(5,6));
con_L_succ_emg_pass = cell2mat(spike_stop_cond(5,10));
con_L_succ_emgvel_pass = cell2mat(spike_stop_cond(5,11));
con_L_fail_force_pass = cell2mat(spike_stop_cond(6,5));
con_L_fail_vel_pass = cell2mat(spike_stop_cond(6,6));
con_L_fail_emg_pass = cell2mat(spike_stop_cond(6,10));
con_L_fail_emgvel_pass = cell2mat(spike_stop_cond(6,11));

%% Recreate Plot from Jack

%now create graphs and then extract the averages etc
%choose right set i.e., only stop trials
%get amplitude so max, duration until it goes back to baseline, and
%area/surface
%decide on appropriate threshold to use for the area under the curve and
%also for reclassification of successful trials, i.e., only those under
%e.g., 1% increase in MVC are ok, maybe I can get average volatility of in
%responses and use that

linewidz = 1.5;
figure;
for i = 1:size(con_R_fail_force, 1)
    subplot(sqrt(size(con_R_fail_force, 1)), sqrt(size(con_R_fail_force, 1)), i);
    plot(time_stoplock,con_R_fail_force (i,:), 'Color', [0 0 0], 'LineWidth', linewidz);hold on;
    plot(time_stoplock,con_R_fail_vel (i,:), 'Color', [0.5 0.5 0.5], 'LineWidth', linewidz);
    plot(time_stoplock,con_L_fail_force (i,:), 'Color', [1 0 0], 'LineWidth', linewidz);
    plot(time_stoplock,con_L_fail_vel (i,:), 'Color', [0.7 0 0], 'LineWidth', linewidz);
    yline(0, '--');
    xlim([-1000 1000]);xlabel('Time from stop siganl (ms)');ylabel('Force');title('Failed stops: Force');
end
legend('rel R fail. force','rel R fail. force vel.','rel L fail. force','rel L fail. force vel.');

%% Some General Variables
n_trial = 26;

%base min level on average volatility, i.e., if responses were more around
%1% MVC then use that rather than 0

%get average volatility as the average of points above 0 in the first 500
%ms
%different methods to see which estimation works best
avg_volatility1 = mean(abs(con_R_fail_force(n_trial, 1:501)));
avg_volatility2 = mean(abs(con_R_fail_force(n_trial, 1:801)));

avg_volatility1vel = mean(abs(con_R_fail_vel(n_trial, 1:501)));
avg_volatility2vel = mean(abs(con_R_fail_vel(n_trial, 1:801)));



%get the stopping amplitude
stopping_amplitude = max(con_R_fail_force(n_trial, 1001:2001));
stopping_amplitude_time = find(con_R_fail_force(n_trial,:) == stopping_amplitude);

%get the stopping duration
stopping_end = min(abs(con_R_fail_force(n_trial, stopping_amplitude_time:2001)));
stopping_end_time = find(abs(con_R_fail_force(n_trial,:)) == stopping_end);

stopping_start_time = 1001;

stopping_duration = stopping_end_time - stopping_start_time;

%get the stopping surface
stopping_surface = trapz(stopping_start_time:stopping_end_time, con_R_fail_force(n_trial, stopping_start_time:stopping_end_time));

%improve on finding the start of the response
%find the x-intercept before the peak response
%use xinterceptvalue = min(abs(con_R_fail_force(1, stopping_start_time:stopping_amplitude_time)))
%then xintercept = find(abs(con_R_fail_force(1,
%stopping_start_time:stopping_amplitude_time)) == xinterceptvalue)
%after use this to find the velocity intercept closest to this value to
%find the starting point

force_data = [stopping_start_time, stopping_amplitude_time, stopping_end_time, stopping_amplitude ,stopping_surface, stopping_duration];

%% Better Force Plotting

%new approach, based on velocity and adjusted for fluctuation

%first get average variation in baseline

avg_volatility_range = 1001;
avg_volatility_range = 1:avg_volatility_range;
avg_volatility = mean(abs(con_R_fail_force(n_trial, avg_volatility_range)));
avg_volatilityvel = mean(abs(con_R_fail_vel(n_trial, avg_volatility_range)));

%now find 0 intercepts of velocity prior to force peak adjusted to new
%baseline
%first adjust and cut variables, note that this cuts to the force peak
%edit: cut to 10 ms before force peak to avoid finding velocity 0 at peak
time_adjust = -10;
con_R_fail_force_adj = con_R_fail_force(n_trial, 1:stopping_amplitude_time+time_adjust) - avg_volatility;
con_R_fail_vel_adj = con_R_fail_vel(n_trial, 1:stopping_amplitude_time+time_adjust) - avg_volatilityvel;
time_stoplock_adj = time_stoplock(1:stopping_amplitude_time+time_adjust);

%then get all the 0s and the last value as this will be cloesest to the
%peak
loc_0_force = find(abs(con_R_fail_force_adj(:)) < 0.01);
force_0_closest_to_peak = loc_0_force(length(loc_0_force));

loc_0_vel = find(abs(con_R_fail_vel_adj(:)) < 0.01);
vel_0_closest_to_peak = loc_0_vel(length(loc_0_vel));

%verify graphically
% figure;
% plot(time_stoplock_adj, con_R_fail_force_adj(:), 'Color', [0 0 0], 'LineWidth', linewidz);hold on;
% plot(time_stoplock_adj, con_R_fail_vel_adj(:), 'Color', [0.5 0.5 0.5], 'LineWidth', linewidz);
% yline(0, '--');
% xline(vel_0_closest_to_peak-1001);

%plot as regular figure
stopping_start_time = vel_0_closest_to_peak;

linewidz = 1.5;
figure;
plot(time_stoplock,con_R_fail_force(n_trial,:), 'Color', [0 0 0], 'LineWidth', linewidz);hold on;
area((stopping_start_time-1001):(stopping_end_time-1001), con_R_fail_force(n_trial, stopping_start_time:stopping_end_time));
plot(time_stoplock,con_R_fail_vel(n_trial,:), 'Color', [0.5 0.5 0.5], 'LineWidth', linewidz);
yline(0, '--');
yline(avg_volatilityvel, '-');
xline(stopping_start_time-1001);
xline(stopping_end_time-1001);
xline(stopping_amplitude_time-1001);
xlim([-1000 1000]);xlabel('Time from stop siganl (ms)');ylabel('Force');title('Failed stops: Force');
legend('rel R fail. force','area R fail. force', 'rel R fail. force vel.');

%%

%plot for all
all_force_data = [];
linewidz = 1.5;
figure;
for i = 1:size(con_R_fail_force, 1)
    n_trial = i;
    stopping_amplitude = max(con_R_fail_force(n_trial, 1001:2001));
    stopping_amplitude_time = find(con_R_fail_force(n_trial,:) == stopping_amplitude);
    stopping_end = min(abs(con_R_fail_force(n_trial, stopping_amplitude_time:2001)));
    stopping_end_time = find(abs(con_R_fail_force(n_trial,:)) == stopping_end);
    stopping_start_time = 1001;
    stopping_duration = stopping_end_time - stopping_start_time; %gets recalculated again below in case values change
    stopping_surface = trapz(stopping_start_time:stopping_end_time, con_R_fail_force(n_trial, stopping_start_time:stopping_end_time));
    
    avg_volatility_range = 1001;
    avg_volatility_range = 1:avg_volatility_range;
    avg_volatility = mean(abs(con_R_fail_force(n_trial, avg_volatility_range)));
    avg_volatilityvel = mean(abs(con_R_fail_vel(n_trial, avg_volatility_range)));
    
    time_adjust = -10;
    con_R_fail_force_adj = con_R_fail_force(n_trial, 1:stopping_amplitude_time+time_adjust) - avg_volatility;
    con_R_fail_vel_adj = con_R_fail_vel(n_trial, 1:stopping_amplitude_time+time_adjust) - avg_volatilityvel;

    loc_0_vel = find(abs(con_R_fail_vel_adj(:)) < 0.01);
    vel_0_closest_to_peak = loc_0_vel(length(loc_0_vel));

    if isempty(vel_0_closest_to_peak)
        disp("skipped");
    else
        stopping_start_time = vel_0_closest_to_peak;
    end

    stopping_duration = stopping_end_time - stopping_start_time; %gets recalculated again below in case values change
    stopping_surface = trapz(stopping_start_time:stopping_end_time, con_R_fail_force(i, stopping_start_time:stopping_end_time));
    
    all_force_data = [all_force_data; stopping_start_time-1001, stopping_amplitude_time-1001, stopping_end_time-1001, stopping_amplitude ,stopping_surface, stopping_duration];

    subplot(sqrt(size(con_R_fail_force, 1)), sqrt(size(con_R_fail_force, 1)), n_trial);
    plot(time_stoplock,con_R_fail_force (n_trial,:), 'Color', [0 0 0], 'LineWidth', linewidz);hold on;
    area((stopping_start_time-1001):(stopping_end_time-1001), con_R_fail_force(n_trial, stopping_start_time:stopping_end_time));
    plot(time_stoplock,con_R_fail_vel (n_trial,:), 'Color', [0.5 0.5 0.5], 'LineWidth', linewidz);
    yline(0, '--');
    yline(avg_volatilityvel, '-');
    xline(stopping_start_time-1001);
    xline(stopping_end_time-1001);
    xline(stopping_amplitude_time-1001);
    xlim([-1000 1000]);xlabel('Time from stop siganl (ms)');ylabel('Force');title('Failed stops: Force');
end
legend('rel R fail. force', 'area R fail. force', 'rel R fail. force vel.');



%% Old Force Plotting

%test force plot

loc_0_force = find(abs(con_R_fail_force(n_trial,:)) < 0.01); %adjust to set different precision for intercepts, althought might also adjust with separate baseline
loc_0_vel = find(abs(con_R_fail_vel(n_trial,:)) < 0.01);
comparison_cutoff = -200+1001; %adjust level here

loc_0_force = loc_0_force(loc_0_force > comparison_cutoff);
loc_0_vel = loc_0_vel(loc_0_vel > comparison_cutoff);

matrix_of_diff_values = abs(loc_0_force(1) - loc_0_vel);
min_diff = min(matrix_of_diff_values);
location_in_diff_matrix = find(matrix_of_diff_values == min_diff);
corresponding_vel_value = loc_0_vel(location_in_diff_matrix);

if isnan(corresponding_vel_value)
    disp("skipped");
elseif corresponding_vel_value+2 > (stopping_amplitude_time) %adjust for amplitude timing
    disp("skipped")
else
    stopping_start_time = corresponding_vel_value;
end

linewidz = 1.5;
figure;
plot(time_stoplock,con_R_fail_force (n_trial,:), 'Color', [0 0 0], 'LineWidth', linewidz);hold on;
area((stopping_start_time-1001):(stopping_end_time-1001), con_R_fail_force(n_trial, stopping_start_time:stopping_end_time));
plot(time_stoplock,con_R_fail_vel (n_trial,:), 'Color', [0.5 0.5 0.5], 'LineWidth', linewidz);
yline(0, '--');
yline(avg_volatility1, '-');
xline(stopping_start_time-1001);
xline(stopping_end_time-1001);
xline(stopping_amplitude_time-1001);
for i = 1:size(loc_0_force, 2)
   xline(loc_0_force(1,i)-1001, '--', 'Color', [0 0 1]); %blue color
end
for i = 1:size(loc_0_vel, 2)
    xline(loc_0_vel(1,i)-1001, '--', 'Color', [1 0 0]); %red color
end
xlim([-1000 1000]);xlabel('Time from stop siganl (ms)');ylabel('Force');title('Failed stops: Force');
legend('rel R fail. force','area R fail. force', 'rel R fail. force vel.');

%%

% plot all force
all_force_data = [];
linewidz = 1.5;
figure;
for i = 1:size(con_R_fail_force, 1)
    stopping_amplitude = max(con_R_fail_force(i, 1001:2001));
    stopping_amplitude_time = find(con_R_fail_force(i,:) == stopping_amplitude);
    stopping_end = min(abs(con_R_fail_force(i, stopping_amplitude_time:2001)));
    stopping_end_time = find(abs(con_R_fail_force(i,:)) == stopping_end);
    stopping_start_time = 1001;
    stopping_duration = stopping_end_time - stopping_start_time; %gets recalculated again below in case values change
    stopping_surface = trapz(stopping_start_time:stopping_end_time, con_R_fail_force(i, stopping_start_time:stopping_end_time));
    
    avg_volatility2vel = mean(abs(con_R_fail_vel(i, 1:801)));

    loc_0_force= find(abs(con_R_fail_force(i,:)) < 0.01); %adjust to set different precision for intercepts
    loc_0_vel = find(abs(con_R_fail_vel(i,:)) < 0.01);
    comparison_cutoff = -200+1001; %adjust level here
    loc_0_force = loc_0_force(loc_0_force > comparison_cutoff);
    loc_0_vel = loc_0_vel(loc_0_vel > comparison_cutoff);
    matrix_of_diff_values = abs(loc_0_force(1) - loc_0_vel);
    min_diff = min(matrix_of_diff_values);
    location_in_diff_matrix = matrix_of_diff_values == min_diff;
    corresponding_vel_value = loc_0_vel(location_in_diff_matrix);
    if isempty(corresponding_vel_value)
        disp("skipped");
    elseif corresponding_vel_value+2 > stopping_amplitude_time %adjust for amplitude timing
        disp("skipped")
    else
        stopping_start_time = corresponding_vel_value;
    end

    stopping_duration = stopping_end_time - stopping_start_time; %gets recalculated again below in case values change
    stopping_surface = trapz(stopping_start_time:stopping_end_time, con_R_fail_force(i, stopping_start_time:stopping_end_time));
    
    all_force_data = [all_force_data; stopping_start_time-1001, stopping_amplitude_time-1001, stopping_end_time-1001, stopping_amplitude ,stopping_surface, stopping_duration];

    subplot(sqrt(size(con_R_fail_force, 1)), sqrt(size(con_R_fail_force, 1)), i);
    plot(time_stoplock,con_R_fail_force (i,:), 'Color', [0 0 0], 'LineWidth', linewidz);hold on;
    area((stopping_start_time-1001):(stopping_end_time-1001), con_R_fail_force(i, stopping_start_time:stopping_end_time));
    plot(time_stoplock,con_R_fail_vel (i,:), 'Color', [0.5 0.5 0.5], 'LineWidth', linewidz);
    yline(0, '--');
    yline(avg_volatility2vel, '-');
    xline(stopping_start_time-1001);
    xline(stopping_end_time-1001);
    xline(stopping_amplitude_time-1001);
    xlim([-1000 1000]);xlabel('Time from stop siganl (ms)');ylabel('Force');title('Failed stops: Force');
end
legend('rel R fail. force', 'area R fail. force', 'rel R fail. force vel.');

%% Old EMG Plotting

%test EMG plot
n_trial = 9;

linewidz = 1.5;
figure;

stopping_amplitude = max(con_R_fail_emg(n_trial, 1001:2001));
stopping_amplitude_time = find(con_R_fail_emg(n_trial,:) == stopping_amplitude);
stopping_end = min(abs(con_R_fail_emg(n_trial, stopping_amplitude_time:2001)));
stopping_end_time = find(abs(con_R_fail_emg(n_trial,:)) == stopping_end);
stopping_start_time = 1001;
stopping_duration = stopping_end_time - stopping_start_time;
stopping_surface = trapz(stopping_start_time:stopping_end_time, con_R_fail_emg(n_trial, stopping_start_time:stopping_end_time));

loc_0_emg= find(abs(con_R_fail_emg(n_trial,:)) < 0.05); %adjust to set different precision for intercepts, althought might also adjust with separate baseline
loc_0_emgvel = find(abs(con_R_fail_emgvel(n_trial,:)) < 0.05);
comparison_cutoff = -300+1001; %adjust level here

loc_0_emg = loc_0_emg(loc_0_emg > comparison_cutoff);
loc_0_emgvel = loc_0_emgvel(loc_0_emgvel > comparison_cutoff);

matrix_of_diff_values = abs(loc_0_emg(1) - loc_0_emgvel);
min_diff = min(matrix_of_diff_values);
location_in_diff_matrix = find(matrix_of_diff_values == min_diff);
corresponding_vel_value = loc_0_emgvel(location_in_diff_matrix);

if isnan(corresponding_vel_value)
    disp("skipped");
elseif corresponding_vel_value+2 > (stopping_amplitude_time) %adjust for amplitude timing
    disp("skipped")
elseif corresponding_vel_value < 401
    disp("skipped")
else
    stopping_start_time = corresponding_vel_value;
end

plot(time_stoplock,con_R_fail_emg (n_trial,:), 'Color', [0 0 0], 'LineWidth', linewidz);hold on;
area((stopping_start_time-1001):(stopping_end_time-1001), con_R_fail_emg(n_trial, stopping_start_time:stopping_end_time));
plot(time_stoplock,con_R_fail_emgvel (n_trial,:), 'Color', [0.5 0.5 0.5], 'LineWidth', linewidz);
yline(0, '--');
xline(stopping_start_time-1001);
xline(stopping_end_time-1001);
xline(stopping_amplitude_time-1001);
for i = 1:size(loc_0_emg, 2)
   xline(loc_0_emg(1,i)-1001, '--', 'Color', [0 0 1]); %blue color
end
for i = 1:size(loc_0_emgvel, 2)
    xline(loc_0_emgvel(1,i)-1001, '--', 'Color', [1 0 0]); %red color
end
xlim([-1000 1000]);xlabel('Time from stop siganl (ms)');ylabel('EMG');title('Failed stops: EMG');
legend('rel R fail. EMG', 'area R fail. EMG', 'rel R fail. EMG vel.');

%% 

%plot all EMG

linewidz = 1.5;
figure;

for i = 1:size(con_R_fail_emg, 1)

    n_trial = i;

    stopping_amplitude = max(con_R_fail_emg(n_trial, 1001:2001));
    stopping_amplitude_time = find(con_R_fail_emg(n_trial,:) == stopping_amplitude);
    stopping_end = min(abs(con_R_fail_emg(n_trial, stopping_amplitude_time:2001)));
    %stopping_end = con_R_fail_emg(abs(con_R_fail_emg(n_trial, :)) < 0.01); %this gives y values, find corresponding x and find closest to amplitude time.
    
    stopping_end_time = find(abs(con_R_fail_emg(n_trial,:)) == stopping_end);
    stopping_start_time = 1001;
    stopping_duration = stopping_end_time - stopping_start_time;
    stopping_surface = trapz(stopping_start_time:stopping_end_time, con_R_fail_emg(n_trial, stopping_start_time:stopping_end_time));

    loc_0_emg= find(abs(con_R_fail_emg(n_trial,:)) < 0.05); %adjust to set different precision for intercepts, althought might also adjust with separate baseline
    loc_0_emgvel = find(abs(con_R_fail_emgvel(n_trial,:)) < 0.05);
    comparison_cutoff = -300+1001; %adjust level here

    loc_0_emg = loc_0_emg(loc_0_emg > comparison_cutoff);
    loc_0_emgvel = loc_0_emgvel(loc_0_emgvel > comparison_cutoff);

    matrix_of_diff_values = abs(loc_0_emg(1) - loc_0_emgvel);
    min_diff = min(matrix_of_diff_values);
    location_in_diff_matrix = find(matrix_of_diff_values == min_diff);
    corresponding_vel_value = loc_0_emgvel(location_in_diff_matrix);

    if isnan(corresponding_vel_value)
        disp("skipped");
    elseif corresponding_vel_value+2 > (stopping_amplitude_time) %adjust for amplitude timing
        disp("skipped")
    elseif corresponding_vel_value < 401
        disp("skipped")
    else
        stopping_start_time = corresponding_vel_value;
    end
    
    subplot(sqrt(size(con_R_fail_emg, 1)), sqrt(size(con_R_fail_emg, 1)), n_trial);
    plot(time_stoplock,con_R_fail_emg (n_trial,:), 'Color', [0 0 0], 'LineWidth', linewidz);hold on;
    area((stopping_start_time-1001):(stopping_end_time-1001), con_R_fail_emg(n_trial, stopping_start_time:stopping_end_time));
    plot(time_stoplock,con_R_fail_emgvel (n_trial,:), 'Color', [0.5 0.5 0.5], 'LineWidth', linewidz);
    yline(0, '--');
    xline(stopping_start_time-1001);
    xline(stopping_end_time-1001);
    xline(stopping_amplitude_time-1001);
    xlim([-1000 1000]);xlabel('Time from stop siganl (ms)');ylabel('EMG');title('Failed stops: EMG');
end
legend('rel R fail. EMG', 'area R fail. EMG', 'rel R fail. EMG vel.');

%% Better EMG Plotting

all_EMG_data = [];
linewidz = 1.5;
figure;
for i = 1:size(con_R_fail_emg, 1)
    n_trial = i;
    stopping_amplitude = max(con_R_fail_emg(n_trial, 1001:2001));
    stopping_amplitude_time = find(con_R_fail_emg(n_trial,:) == stopping_amplitude);
    stopping_end = min(abs(con_R_fail_emg(n_trial, stopping_amplitude_time:2001)));

    % alternative approaches below decrease performance

    % stopping_end = con_R_fail_emg(n_trial, abs(con_R_fail_emg(n_trial, stopping_amplitude_time:2001)) < 0.5);

    % avg_volatility_range = 1001;
    % avg_volatility_range = 1:avg_volatility_range;
    % avg_volatility = mean(abs(con_R_fail_emg(n_trial, avg_volatility_range)));
    % con_R_fail_emg_adj = con_R_fail_emg(n_trial,:) - avg_volatility;
    % stopping_end = min(abs(con_R_fail_emg_adj(stopping_amplitude_time:2001)));
    % stopping_end_time = find(abs(con_R_fail_emg_adj(:)) == stopping_end);
    % stopping_end = con_R_fail_emg(n_trial, stopping_end_time);
    
    stopping_end_time = find(abs(con_R_fail_emg(n_trial,:)) == stopping_end);
    stopping_start_time = 1001;
    stopping_duration = stopping_end_time - stopping_start_time; %gets recalculated again below in case values change
    stopping_surface = trapz(stopping_start_time:stopping_end_time, con_R_fail_emg(n_trial, stopping_start_time:stopping_end_time));
    
    avg_volatility_range = 1001;
    avg_volatility_range = 1:avg_volatility_range;
    avg_volatility = mean(abs(con_R_fail_emg(n_trial, avg_volatility_range)));
    avg_volatilityvel = mean(abs(con_R_fail_emgvel(n_trial, avg_volatility_range)));
    
    time_adjust = -20;
    con_R_fail_emg_adj = con_R_fail_emg(n_trial, 1:stopping_amplitude_time+time_adjust) - avg_volatility;
    con_R_fail_emgvel_adj = con_R_fail_emgvel(n_trial, 1:stopping_amplitude_time+time_adjust) - avg_volatilityvel;

    loc_0_vel = find(abs(con_R_fail_emgvel_adj(:)) < 0.5);

    cutoff_for_multiple_peaks = 1001;
    loc_0_vel = loc_0_vel(loc_0_vel < cutoff_for_multiple_peaks);
    vel_0_closest_to_peak = loc_0_vel(length(loc_0_vel));

    if isempty(vel_0_closest_to_peak)
        disp("skipped");
    else
        stopping_start_time = vel_0_closest_to_peak;
    end

    stopping_duration = stopping_end_time - stopping_start_time; %gets recalculated again below in case values change
    stopping_surface = trapz(stopping_start_time:stopping_end_time, con_R_fail_emg(i, stopping_start_time:stopping_end_time));
    
    all_EMG_data = [all_EMG_data; stopping_start_time-1001, stopping_amplitude_time-1001, stopping_end_time-1001, stopping_amplitude ,stopping_surface, stopping_duration];

    subplot(sqrt(size(con_R_fail_emg, 1)), sqrt(size(con_R_fail_emg, 1)), n_trial);
    plot(time_stoplock,con_R_fail_emg (n_trial,:), 'Color', [0 0 0], 'LineWidth', linewidz);hold on;
    area((stopping_start_time-1001):(stopping_end_time-1001), con_R_fail_emg(n_trial, stopping_start_time:stopping_end_time));
    plot(time_stoplock,con_R_fail_emgvel (n_trial,:), 'Color', [0.5 0.5 0.5], 'LineWidth', linewidz);
    yline(0, '--');
    yline(avg_volatilityvel, '-');
    xline(stopping_start_time-1001);
    xline(stopping_end_time-1001);
    xline(stopping_amplitude_time-1001);
    xlim([-1000 1000]);xlabel('Time from stop siganl (ms)');ylabel('EMG');title('Failed stops: EMG');
end
legend('rel R fail. EMG', 'area R fail. EMG', 'rel R fail. EMG vel.');

%% End


