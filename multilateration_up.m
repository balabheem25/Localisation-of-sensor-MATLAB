clc; clear; close all;

clc; clear; close all;

startTime = datetime('now');


%% CONFIGURATION
% Define file paths dynamically
file_paths = {
    "C:\Users\karth\Downloads\Localization mat\Kartheek\grid_T1_R1\data_X10Y-10.csv", ...
    "C:\Users\karth\Downloads\Localization mat\Kartheek\grid_T2_R2\data_X10Y-10.csv", ...
    "C:\Users\karth\Downloads\Localization mat\Kartheek\grid_T3_R3\data_X10Y-10.csv", ...
    "C:\Users\karth\Downloads\Localization mat\Kartheek\grid_T6_R6\data_X10Y-10.csv", ...
    "C:\Users\karth\Downloads\Localization mat\Kartheek\grid_ref\ref_T1_R1_data_X0Y0.csv", ...
    "C:\Users\karth\Downloads\Localization mat\Kartheek\grid_ref\ref_T2_R2_data_X0Y0.csv", ...
    "C:\Users\karth\Downloads\Localization mat\Kartheek\grid_ref\ref_T3_R3_data_X0Y0.csv",...
    "C:\Users\karth\Downloads\Localization mat\Kartheek\grid_ref\ref_T6_R6_data_X0Y0.csv"
};
%% commented file paths for transducer 5
%"C:\Users\karth\Downloads\Localization mat\Kartheek\grid_T5_R5\data_X-20Y30.csv", ...
%%"C:\Users\karth\Downloads\Localization mat\Kartheek\grid_ref\ref_T5_R5_data_X0Y0.csv",...

% Preprocessing parameters
f_sig = 1e6; % Signal frequency
fs_rs_factor = 50; % Resampling factor
lowcut = 0.8e6; % Bandpass lower cutoff frequency
highcut = 1.2e6; % Bandpass upper cutoff frequency
order = 4; % Filter order
envelope_window = 50; % Envelope window size

cylinder_radius = 55; %mm(outer radius of cylinder)
thickness = 4;       % mm (wall thickness)
r = (cylinder_radius - thickness); % usable inner radius = 51 mm
v = 1.48; % speed of sound in water, mm/µs
sensor_radius = 12; % mm (sensor diameter = 24mm)

%% Sensor positions (in mm)
sensor_angles = [ ...
    270, ...            % Transducer 1
    326.25, ...         % Transducer 2
    22.5, ...           % Transducer 3
    90, ...             % Transducer 4
    157.5, ...          % Transducer 5
    213.75 ...          % Transducer 6
    
];

% --- Compute all sensor positions (Nx2 array) ---
num_sensors = length(sensor_angles);
sensor_positions = zeros(num_sensors, 2);
for k = 1:num_sensors
    sensor_positions(k, :) = [r * cosd(sensor_angles(k)), r * sind(sensor_angles(k))];
end

% --- Choose which sensors to use (by index, e.g., [1 2 3 6]) ---
selected_indices = [1 2 3 6]; % Example: use sensors 1, 2, 3, 6
% selected_indices = [1 3 6]; % use sensors 1, 3, 6
selected_positions = sensor_positions(selected_indices, :);

%% Plot the setup 

%plot_setup(r, cylinder_radius, num_sensors, sensor_angles);

%% Loop through files and process signals
global  processed_signals
global proceesed_signal_burst
processed_signals = struct(); % Store processed signals
burst_input_idx = 2;

%% Process burst_input_reference from file , colummn 2
burst_file = file_paths{4};
burst_data = readmatrix(burst_file);

t_burst = burst_data(:, 1);
burst_signal = burst_data(:, 2);

% Resample
fs_rs = f_sig * fs_rs_factor;
t_rs_burst = (t_burst(1):1 / fs_rs:t_burst(end));
burst_signal_rs = interp1(t_burst, burst_signal, t_rs_burst, 'spline');

% Normalize
burst_signal_n = burst_signal_rs ./ max(burst_signal_rs);

% Bandpass filter
Wn = [lowcut highcut] / (fs_rs / 2);
[b, a] = butter(order, Wn, 'bandpass');
burst_signal_f = filtfilt(b, a, burst_signal_n);

% Envelope
burst_signal_env = envelope(burst_signal_f, envelope_window, 'peak');
burst_25_idx = t_rs_burst <= 25e-6; %25microseconds
burst_signal_env_25us = burst_signal_env(burst_25_idx);
t_rs_burst_25us = t_rs_burst(burst_25_idx);


% store this first 25us in processed_signal_burst
burst_idx_psb = 1;
processed_signal_burst(burst_idx_psb).file = burst_file;
processed_signal_burst(burst_idx_psb).time = t_rs_burst_25us;
processed_signal_burst(burst_idx_psb).resampled = burst_signal_rs(burst_25_idx);
processed_signal_burst(burst_idx_psb).normalized = burst_signal_n(burst_25_idx);
processed_signal_burst(burst_idx_psb).filtered = burst_signal_f(burst_25_idx);
processed_signal_burst(burst_idx_psb).envelope = burst_signal_env_25us;


%assign to burst_input_reference for convenience
burst_input_reference = processed_signal_burst(burst_idx_psb);
burst_peak_psb = max(burst_input_reference.envelope);

for i = 1:length(file_paths)
    % Load file
    file = file_paths{i};
    fprintf('Processing file: %s\n', file);
    data = readmatrix(file);

    % Extract time and channels
    t = data(:, 1); % Time column
    dt = t(2) - t(1); % Time step
    fs = 1 / dt; % Original sampling frequency

    signal = data(:, 3);

    %% Resampling
    fs_rs = f_sig * fs_rs_factor; % Resampled frequency
    t_rs = (t(1):1/fs_rs:t(end))'; % Resampled time vector
    signal_rs = interp1(t, signal, t_rs, 'spline'); % Resample signal

    %% Normalization
    signal_n = signal_rs ./ max(signal_rs); % Normalize signal

    %% Bandpass Filtering
    Wn = [lowcut highcut] / (fs_rs / 2); % Normalized cutoff frequencies
    [b, a] = butter(order, Wn, 'bandpass'); % Bandpass filter design
    signal_f = filtfilt(b, a, signal_n); % Apply filter (phase-neutral)

    %% Envelope Detection
    signal_env = envelope(signal_f, envelope_window, 'peak'); % Compute envelope
    %burst_peak = max(burst_input_reference.envelope);
    envelope_max = max(signal_env);
    is_valid_signal = envelope_max >= (0.05 * burst_peak_psb);  % Flag if envelope is strong enough
    

    %% Store processed signals
    processed_signals(i).file = file;
    processed_signals(i).time = t_rs;
    processed_signals(i).resampled = signal_rs;
    processed_signals(i).normalized = signal_n;
    processed_signals(i).filtered = signal_f;
    processed_signals(i).envelope = signal_env;
    processed_signals(i).envelope_max = envelope_max;
    processed_signals(i).is_valid = is_valid_signal;

  
    [category, max_after140, mean_after140] = categorize_signal_after140us(t_rs, signal_env);
    processed_signals(i).after140_category = category;
    processed_signals(i).after140_max = max_after140;
    processed_signals(i).after140_mean = mean_after140;


   
    

    %% Plot Results
    figure(i);
    hold on; grid on;
    plot(t_rs, signal_f, 'LineWidth', 2, 'DisplayName', 'Filtered Signal');
    plot(t_rs, signal_env, 'LineWidth', 2, 'DisplayName', 'Envelope');
    xlabel('Time [s]');
    ylabel('Amplitude [V]');
    legend('Location', 'northeast', 'FontSize', 8);
    title(sprintf('Processed Signal for File %d', i));
    set(gca, 'FontSize', 16, 'FontWeight', 'bold');
    hold off;
end

%-------------------------------------------------------
%% Process burst_input_reference from file , colummn 2
burst_file = file_paths{4};
burst_data = readmatrix(burst_file);

t_burst = burst_data(:, 1);
burst_signal = burst_data(:, 2);

% Resample
fs_rs = f_sig * fs_rs_factor;
t_rs_burst = (t_burst(1):1 / fs_rs:t_burst(end));
burst_signal_rs = interp1(t_burst, burst_signal, t_rs_burst, 'spline');

% Normalize
burst_signal_n = burst_signal_rs ./ max(burst_signal_rs);

% Bandpass filter
Wn = [lowcut highcut] / (fs_rs / 2);
[b, a] = butter(order, Wn, 'bandpass');
burst_signal_f = filtfilt(b, a, burst_signal_n);

% Envelope
burst_signal_env = envelope(burst_signal_f, envelope_window, 'peak');
burst_25_idx = t_rs_burst <= 25e-6; %25microseconds
burst_signal_env_25us = burst_signal_env(burst_25_idx);
t_rs_burst_25us = t_rs_burst(burst_25_idx);

% Store in processed_signals as a new entry
% burst_idx_full = length(processed_signals) + 1;
% processed_signals(burst_idx).file = burst_file;
% processed_signals(burst_idx).time = t_rs_burst;
% processed_signals(burst_idx).resampled = burst_signal_rs;
% processed_signals(burst_idx).normalized = burst_signal_n;
% processed_signals(burst_idx).filtered = burst_signal_f;
% processed_signals(burst_idx).envelope = burst_signal_env;

% store this first 25us in processed signals
burst_idx = length(processed_signals) + 1;
processed_signals(burst_idx).file = burst_file;
processed_signals(burst_idx).time = t_rs_burst_25us;
processed_signals(burst_idx).resampled = burst_signal_rs(burst_25_idx);
processed_signals(burst_idx).normalized = burst_signal_n(burst_25_idx);
processed_signals(burst_idx).filtered = burst_signal_f(burst_25_idx);
processed_signals(burst_idx).envelope = burst_signal_env_25us;


%assign to burst_input_reference for convenience
burst_input_reference = processed_signals(burst_idx);
burst_peak = max(burst_input_reference.envelope);

%% (don't delete this section) plot the input burst_signal which is considered as the reference
% % % % Plot the burst_input_reference envelope
% % % figure;
% % % hold on; grid on;
% % % plot(t_rs_burst, burst_signal_f, 'LineWidth', 2, 'DisplayName', 'Filtered Burst Reference');
% % % plot(t_rs_burst, burst_signal_env, 'LineWidth', 2, 'DisplayName', 'Envelope');
% % % xlabel('Time [s]');
% % % ylabel('Amplitude [V]');
% % % legend('Location', 'northeast', 'FontSize', 8);
% % % title('Processed Burst Input Reference (File 4, Column 2)');
% % % set(gca, 'FontSize', 16, 'FontWeight', 'bold');
% % % hold off;


%---------------------------------------------------------

%% Save Processed Signals to Workspace
assignin('base', 'processed_signals', processed_signals);
fprintf('Processed signals saved to workspace as "processed_signals".\n');


%% Verify Processed Signals
disp('Processed signals structure:');
disp(processed_signals);
fprintf('Number of processed signals: %d\n', length(processed_signals));

if length(processed_signals) < 4
    error('Processed signals structure does not contain enough elements. Ensure all files are processed correctly.');
end
%% Call the function tof_us = extract_tof
%Extract ToF using index values of signals and references in
% "processed signals"
num_selected = length(selected_indices);

% Assume: first N files are signals, next N files are references
signal_indices    = 1:num_selected;
reference_indices = num_selected + (1:num_selected);


tofs = zeros(num_selected, 1);
dists = zeros(num_selected, 1);

% % % % % % % for i = 1:num_selected
% % % % % % %     tofs(i) = extract_tof(signal_indices(i), reference_indices(i), burst_idx, processed_signals);
% % % % % % %     dists(i) = (tofs(i) / 2) * v;
% % % % % % % end

valid_mask = false(num_selected, 1);  % Track which sensors are valid
tofs = NaN(num_selected, 1);  % Default to NaN for invalids
dists = NaN(num_selected, 1);

for i = 1:num_selected
    sig_idx = signal_indices(i);
    if ~processed_signals(sig_idx).is_valid
        fprintf('Skipping Sensor %d: Envelope max %.4f < 0.05\n', ...
                selected_indices(i), processed_signals(sig_idx).envelope_max);
        continue;
    end

    tofs(i) = extract_tof(sig_idx, reference_indices(i), burst_idx, processed_signals);
    % --- NEW: Skip if ToF is NaN ---
    if isnan(tofs(i))
        fprintf('Skipping Sensor %d: Invalid ToF (env after 140 µs >= 1)\n', ...
                selected_indices(i));
        continue;
    end

    dists(i) = (tofs(i) / 2) * v;
    valid_mask(i) = true;
end

% Use only valid signals for position estimation
valid_positions = selected_positions(valid_mask, :);
valid_dists = dists(valid_mask);




% Now dists is a vector of one-way distances for each signal-reference pair
% Select only those used in selected_indices (make sure selected_indices does not exceed num_pairs)
selected_dists = dists;  % dists = [(tof1/2)*v; (tof2/2)*v; ...] are already in the order of selected indices

% Print ToF and distances for each pair
for i = 1:num_selected
    fprintf('ToF%d: %.2f µs, Distance%d: %.2f mm\n', i, tofs(i), i, dists(i));
end




%% Multilateration region approach (use for any number of sensors)
[roi_x, roi_y] = multilateration_region(selected_positions, selected_dists, r);
fprintf('Multilateration region computed. Number of probable points: %d\n', numel(roi_x));

%% Plot
num_selected = length(selected_indices);
locus = cell(num_selected, 1);
locus_truncated = false(num_selected, 1);

% Get after-140us categories for each selected signal
categories = cell(num_selected, 1);
for i = 1:num_selected
    s = selected_positions(i, :);
    d = selected_dists(i);
    %categories{i} = processed_signals(i).after140_category;
    locus{i} = compute_locus(s, d, r);
    locus_truncated(i) = any((locus{i}(:,1).^2 + locus{i}(:,2).^2) > r^2);
end

% fig 7
figure; hold on; axis equal; grid on;
plot(r * cos(linspace(0,2*pi,100)), r * sin(linspace(0,2*pi,100)), '-k', 'DisplayName', 'Cylinder Boundary');
scatter(selected_positions(:,1), selected_positions(:,2), 80, 'b', 'filled', 'DisplayName', 'Transducers/ Receivers');
% % % for k = find(valid_mask)'           %1:num_selected
% % %     color = lines(num_selected);
% % %     x = locus{k}(:,1);
% % %     y = locus{k}(:,2);
% % % 
% % %     % % % angles = atan2(y - s(2), x - s(1));
% % %     % % % [~, sort_idx] = sort(angles);
% % %     % % % plot(x(sort_idx), y(sort_idx), 'Color', color(k,:), 'LineWidth', 2, ...
% % %     % % %     'DisplayName', sprintf('Locus %d (S%d)', k, selected_indices(k)));
% % %     plot(x, y, 'Color', color(k,:), 'LineWidth', 2, ...
% % %         'DisplayName', sprintf('Locus %d (S%d)', k, selected_indices(k)));
% % % end


% After computing locus{i} for all sensors
valid_locus_mask = true(num_selected, 1);
for i = 1:num_selected
    if isempty(locus{i}) || any(isnan(locus{i}(:))) || selected_dists(i) <= 0
        valid_locus_mask(i) = false;
        fprintf('Skipping Sensor %d: invalid locus or distance\n', selected_indices(i));
    end
end

% Combine with distance validity
valid_mask = valid_locus_mask;

% Now filter everything based on valid_mask
filtered_positions = selected_positions(valid_mask, :);
filtered_dists = selected_dists(valid_mask);
filtered_locus = locus(valid_mask);
filtered_indices = selected_indices(valid_mask);

color = lines(num_selected);
for k = 1:length(filtered_locus)
    x = filtered_locus{k}(:,1);
    y = filtered_locus{k}(:,2);
    plot(x, y, 'Color', color(k,:), 'LineWidth', 2, ...
        'DisplayName', sprintf('Locus %d (S%d)', k, filtered_indices(k)));
end


if num_selected == 2
    % --- Two sensors: plot intersection or fallback to mean ---
    [intersection_x, intersection_y] = compute_intersection(locus{1}, locus{2});
    if ~isempty(intersection_x)
        scatter(intersection_x, intersection_y, 150, 'k', 'filled', 'DisplayName', 'Intersection Points');
        for i = 1:length(intersection_x)
            text(intersection_x(i), intersection_y(i), sprintf('P%d', i), 'VerticalAlignment', 'top');
        end
    else
        [Xs, Ys] = most_probable_position(locus{1}, locus{2});
        scatter(Xs, Ys, 100, 'm', 'filled', 'DisplayName', 'Most Probable Mean Position');
        text(Xs, Ys, 'Mean Position', 'VerticalAlignment', 'top');
    end
    legend('show', 'Location', 'best');
    % --- Two sensors: use category-based intersection logic ---
    cat1 = categories{1};
    cat2 = categories{2};
    
    pos_1 = selected_indices(1);
    pos_2 = selected_indices(2);
    fprintf('Category for Sensor {%d}: %s\n', pos_1, cat1);
    fprintf('Category for Sensor {%d}: %s\n', pos_2, cat2);

    if strcmp(cat1, "center_line") && strcmp(cat2, "center_line")
        disp('Both signals are center_line.');
    elseif strcmp(cat1, "outside_cone") && strcmp(cat2, "outside_cone")
        disp('Both signals are outside_cone.');
    elseif strcmp(cat1, "center_line") && (strcmp(cat2, "on_cone_line") || strcmp(cat2, "between_center_and_cone"))
        disp('Sensor %d is center_line, Sensor %d is on/between cone.', pos_1, pos_2);
    elseif strcmp(cat2, "center_line") && (strcmp(cat1, "on_cone_line") || strcmp(cat1, "between_center_and_cone"))
        disp('Sensor %d is center_line, Sensor %d is on/between cone.', pos_2, pos_1);
    elseif strcmp(cat1, "on_cone_line") && strcmp(cat2, "on_cone_line")
        disp('Both signals are on_cone_line.');
    else
        disp('Other category combination.');
    end

else
    % --- More than two sensors: multilateration region ---
    [roi_x, roi_y] = multilateration_region(filtered_positions, filtered_dists, r);
    scatter(roi_x, roi_y, 10, 'm', 'filled', 'DisplayName', 'Most Probable Region');
    % Compute the centroid (mean) of the probable points as a single estimate
    % if ~isempty(roi_x)
    %     centroid_x = mean(roi_x);
    %     centroid_y = mean(roi_y);
    % 
    %     % Plot the centroid with a distinct marker and label
    %     scatter(centroid_x, centroid_y, 150, 'r', 'filled', 'DisplayName', 'Estimated Position');
    %     text(centroid_x, centroid_y, 'Estimated Position', 'VerticalAlignment', 'bottom', 'Color', 'r', 'FontWeight', 'bold');
    % else
    %     disp('No probable multilateration region found to estimate position.');
    % end

    legend('show', 'Location', 'best');
end
hold off;


% --- Generalized locus plotting for all selected sensors ---
num_selected = length(selected_indices);
locus = cell(num_selected, 1);
locus_truncated = false(num_selected, 1);

for i = 1:num_selected
    s = selected_positions(i, :);
    d = selected_dists(i);
    locus{i} = compute_locus(s, d, r);
    locus_truncated(i) = any((locus{i}(:,1).^2 + locus{i}(:,2).^2) > r^2);
end



%% --- Print after-140us signal categories for signals 1 & 2 ---
% % % % % fprintf('\nSummary of signal categories after 140us:\n');
% % % % % for i = 1:num_selected
% % % % %     idx = selected_indices(i);
% % % % %     cat = processed_signals(idx).after140_category;
% % % % %     maxv = processed_signals(idx).after140_max;
% % % % %     meanv = processed_signals(idx).after140_mean;
% % % % %     fprintf('Signal %d (Sensor %d): Category = %s, Max = %.3f, Mean = %.3f\n', i, idx, cat, maxv, meanv);
% % % % % end

%% Time stamp/ processing time
fprintf('Started at: %s\n', datestr(startTime));
endTime = datetime('now');
elapsedTime = endTime - startTime;

fprintf('Ended at: %s\n', datestr(endTime));
fprintf('Elapsed time: %s\n', char(elapsedTime));

%% FUNCTIONS
function plot_setup = plot_setup(r, cylinder_radius, num_sensors, sensor_angles)
    
    % --- Plot inner circle (radius 51 mm) and outer 32-gon (radius 55 mm) with transducers at face centers ---
    
    figure; hold on; axis equal; grid on;
    
    % Inner circle (usable region, radius 51 mm)
    theta = linspace(0, 2*pi, 1000);
    h1 = plot(r * cos(theta), r * sin(theta), '-b', 'LineWidth', 2, 'DisplayName', 'Inner Circle (r=51mm)');
    
    % Outer 32-gon (polygonal cylinder boundary, radius 55 mm)
    n_faces = 32;
    theta_poly = linspace(0, 2*pi, n_faces+1) + pi/n_faces; % rotate so a face is at x=0
    x_poly = cylinder_radius * cos(theta_poly); % cylinder_radius = 55
    y_poly = cylinder_radius * sin(theta_poly);
    h2 = plot(x_poly, y_poly, 'k-', 'LineWidth', 2, 'DisplayName', 'Outer 32-gon (r=55mm)');
    
    % Compute face centers (midpoint between adjacent vertices)
    face_centers_x = zeros(1, n_faces);
    face_centers_y = zeros(1, n_faces);
    for i = 1:n_faces
        x1 = cylinder_radius * cos(theta_poly(i));
        y1 = cylinder_radius * sin(theta_poly(i));
        x2 = cylinder_radius * cos(theta_poly(i+1));
        y2 = cylinder_radius * sin(theta_poly(i+1));
        face_centers_x(i) = (x1 + x2) / 2;
        face_centers_y(i) = (y1 + y2) / 2;
    end
    
    % For each sensor, find the closest face center and plot the transducer there
    h3 = gobjects(num_sensors,1);
    for k = 1:num_sensors
        sensor_angle_rad = wrapTo2Pi(deg2rad(sensor_angles(k)));
        [~, idx] = min(abs(wrapTo2Pi(atan2(face_centers_y, face_centers_x)) - sensor_angle_rad));
        tx = face_centers_x(idx);
        ty = face_centers_y(idx);
        h3(k) = scatter(tx, ty, 100, 'r', 'filled', 'DisplayName', sprintf('Transducer %d', k));
        text(tx, ty, sprintf('T%d', k), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    end
    
    legend([h1 h2 h3(1)], {'Inner Circle (r=51mm)', 'Outer 32-gon (r=55mm)', 'Transducers'}, 'Location', 'best');
    xlabel('X (mm)');
    ylabel('Y (mm)');
    title('Cylinder, Transducers/ Receivers setup');
    hold off;
end

function tof_us = extract_tof(signal_idx, ref_idx, burst_idx, processed_signals)
% Estimate ToF as the average of first rise and score peak in the score curve.
% Plots both points and their average.

    % --- Extract signals ---
    signal_data = processed_signals(signal_idx);
    ref_data = processed_signals(ref_idx);
    burst_data = processed_signals(burst_idx);

    t_us = signal_data.time * 1e6; % microseconds
    signal_env = signal_data.envelope;
    ref_env = ref_data.envelope;
    burst_env = burst_data.envelope;

    % --- Subtract reference envelope from signal envelope ---
    env_sub = signal_env - ref_env;
    %env_sub = signal_env;

    % --- Truncate signals to region of interest ---
    min_time = 45; % us
    max_time = 140; % us
    valid_idx = (t_us >= min_time) & (t_us <= max_time);
    t_us_truncated = t_us(valid_idx);
    env_sub_truncated = env_sub(valid_idx);
    env_140 = env_sub(t_us > max_time);

    % === [NEW] Amplitude rejection after 140 µs ===
    idx_140 = t_us > 140;
    
    if any(idx_140)
        % Extract envelope parts after 140 µs (no subtraction, raw values)
        signal_after140 = signal_env(idx_140);
        ref_after140    = ref_env(idx_140);
    
        % Peak values in this region
        signal_peak140 = max(signal_after140);
        ref_peak140    = max(ref_after140);
    
        % Reject if signal is stronger than or equal to reference
        if signal_peak140 >= ref_peak140
            fprintf('Signal rejected: signal env peak %.3f ≥ ref env peak %.3f after 140 µs\n', ...
                     signal_peak140, ref_peak140);
            tof_us = NaN;
            return;
        end
    end

    % --- Set window size in microseconds ---
    desired_window_us = 25; % Adjust as needed
    dt_us = t_us_truncated(2) - t_us_truncated(1);
    window_len = round(desired_window_us / dt_us);

    % --- Prepare burst reference (normalize and truncate/pad if needed) ---
    burst_env = burst_env(:);
    burst_env = burst_env / max(abs(burst_env));
    if length(burst_env) > window_len
        burst_env = burst_env(1:window_len);
    elseif length(burst_env) < window_len
        burst_env = [burst_env; zeros(window_len - length(burst_env), 1)];
    end

    % --- Sliding window std and coherence ---
    num_windows = length(env_sub_truncated) - window_len + 1;
    std_vals = zeros(num_windows,1);
    coh_vals = zeros(num_windows,1);
    fs = 1/(t_us_truncated(2)-t_us_truncated(1)); % MHz
    amp_thresh = 0.1 * max(env_sub_truncated);
    mean_amp = zeros(num_windows, 1);
    for i = 1:num_windows
        win = env_sub_truncated(i:i+window_len-1);
        mean_amp(i) = mean(abs(win));
        if mean(abs(win)) < amp_thresh
            std_vals(i) = NaN;
            coh_vals(i) = NaN;
            continue;
        end
        win = win / max(abs(win)); % normalize window
        std_vals(i) = std(win - burst_env);
        [cxy, f] = mscohere(win, burst_env, [], [], [], fs);
        freq_band = (f >= 0.7) & (f <= 1.3);
        coh_vals(i) = mean(cxy(freq_band));
    end
    std_norm = (std_vals - nanmin(std_vals)) / (nanmax(std_vals) - nanmin(std_vals) + eps);
    coh_norm = (coh_vals - nanmin(coh_vals)) / (nanmax(coh_vals) - nanmin(coh_vals) + eps);
    score = -std_norm + coh_norm;

    % --- Mask out low amplitude and edges ---
    score(mean_amp < amp_thresh) = -Inf;
    edge_margin = 5;
    center_times = t_us_truncated(1 : num_windows) + (window_len / 2) * (t_us_truncated(2) - t_us_truncated(1));
    valid_edge = (center_times > min_time + edge_margin) & (center_times < max_time - edge_margin);
    score(~valid_edge) = -Inf;

    % --- Find best match (score peak) ---
    if all(~isfinite(score))
        best_idx = NaN;
    else
        [~, best_idx] = max(score);
    end
    % --- Envelope peak ---
    [~, env_peak_idx] = max(env_sub_truncated);
    tof_env_peak = t_us_truncated(env_peak_idx);
    
    %% conisdering tof_us from env_peak no matter the condition (can change later)
    % % % % tof_us = tof_env_peak;
    %%Use score peak instead of envelope peak
    if ~isnan(best_idx)
        tof_us = t_us_truncated(best_idx);
    else
        tof_us = NaN;
    end

    fprintf('ToF (envelope peak): %.2f us\n', tof_us);

    % --- Plot for visualization ---
    figure;
    subplot(4,1,1);
    plot(t_us_truncated, env_sub_truncated, 'b'); hold on;
    plot(t_us_truncated(1:window_len), burst_env, 'r');
    title('Subtracted Envelope and Burst Envelope');
    legend('Subtracted', 'Burst Ref');

    subplot(4,1,2);
    valid_std = ~isnan(std_vals);
    plot(t_us_truncated(1:num_windows), std_vals, 'k'); hold on;
    plot(t_us_truncated(find(valid_std)), std_vals(valid_std), 'b', 'LineWidth', 1); % overlays valid values in blue
    %plot(t_us_truncated(1:num_windows), std_vals, 'k'); hold on;
    title('Sliding Window Std Dev');
    ylabel('Std Dev');
    xlabel('Time (\mus)');
    if ~isnan(best_idx)
        scatter(t_us_truncated(best_idx), std_vals(best_idx), 100, 'm', 'filled');
    end
    legend('std_dev vals', 'least std_dev val');

    subplot(4,1,3);
    valid_coh = ~isnan(coh_vals);
    plot(t_us_truncated(1:num_windows), coh_vals, 'g'); hold on;
    plot(t_us_truncated(find(valid_coh)), coh_vals(valid_coh), 'b', 'LineWidth', 1.5);
    %plot(t_us_truncated(1:num_windows), coh_vals, 'g'); hold on;
    title('Sliding window coherence');
    ylabel('coherence');
    xlabel('Time (\mus)');
    if ~isnan(best_idx)
        scatter(t_us_truncated(best_idx), coh_vals(best_idx), 100, 'r', 'filled');
    end
    legend('coherence vals', 'highest coherence val');

    subplot(4,1,4);
    plot(t_us_truncated, env_sub_truncated, 'b'); hold on;
    if ~isnan(best_idx)
        idx = best_idx:min(best_idx+window_len-1, length(env_sub_truncated));
        plot(t_us_truncated(idx), env_sub_truncated(idx), 'r', 'LineWidth', 2);
        scatter(t_us_truncated(best_idx), env_sub_truncated(best_idx), 100, 'r', 'filled', 'DisplayName', 'Score Peak');
    end
    scatter(tof_env_peak, env_sub_truncated(env_peak_idx), 80, 'g', 'filled', 'DisplayName', 'Envelope Peak');
    title('Envelope with Peaks');
    legend('signal', 'most resembling window', 'score peak', 'envelope peak');
    xlabel('Time (\mus)');

    fprintf('ToF (score peak): %.2f us\n', tof_us);
end


function locus = compute_locus(sensor, distance, radius)
    if distance <= 0
        locus = []; % Invalid locus is distance is non-positive
        %return;
    end
    % Compute the locus (arc) of possible positions
    theta = linspace(0, 2*pi, 361);
    %theta(end) = [];
    x = sensor(1) + distance * cos(theta);
    y = sensor(2) + distance * sin(theta);

    % Keep points inside the cylinder
    inside_idx = (x.^2 + y.^2) <= radius^2;
    locus = [x(inside_idx)', y(inside_idx)'];

  
    %Debugging: print the number of points inside the cylinder
    fprintf('Debug: %d points retained inside the cylinder boundaray.\n', sum(inside_idx));
end

function [x, y] = most_probable_position(locus1, locus2)
    % Check if loci are empty
    if isempty(locus1)
        warning('Locus 1 is empty. Cannot compute the most probable position.');
        x = NaN;
        y = NaN;
        return;
    end
    if isempty(locus2)
        warning('Locus 2 is empty. Cannot compute the most probable position.');
        x = NaN;
        y = NaN;
        return;
    end

    % Initialize variables
    min_dist = inf;
    closest_points = []; % Initialize to avoid undefined variable error

    % Find the closest points between the two loci
    for i = 1:size(locus1, 1)
        for j = 1:size(locus2, 1)
            dist = norm(locus1(i,:) - locus2(j,:));
            if dist < min_dist
                min_dist = dist;
                closest_points = [locus1(i,:); locus2(j,:)];
            end
        end
    end

    % Check if closest points were found
    if isempty(closest_points)
        warning('No closest points found between loci. Returning NaN.');
        x = NaN;
        y = NaN;
        return;
    end

    % Compute the midpoint of the closest points
    x = mean(closest_points(:,1));
    y = mean(closest_points(:,2));

    % Debugging message
    fprintf('Most probable position computed: X = %.2f mm, Y = %.2f mm\n', x, y);
end

function [intersection_x, intersection_y] = compute_intersection(locus1, locus2)
    % Initialize intersection points
    intersection_x = [];
    intersection_y = [];
    
    % Loop through all pairs of line segments from locus1 and locus2
    for i = 1:size(locus1, 1) - 1
        for j = 1:size(locus2, 1) - 1
            % Extract line segment endpoints
            x1 = locus1(i, 1); y1 = locus1(i, 2);
            x2 = locus1(i + 1, 1); y2 = locus1(i + 1, 2);
            x3 = locus2(j, 1); y3 = locus2(j, 2);
            x4 = locus2(j + 1, 1); y4 = locus2(j + 1, 2);
            
            % Compute intersection using line-segment intersection formula
            denom = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);
            if denom == 0
                % Lines are parallel, no intersection
                continue;
            end
            
            % Compute intersection point
            px = ((x1 * y2 - y1 * x2) * (x3 - x4) - (x1 - x2) * (x3 * y4 - y3 * x4)) / denom;
            py = ((x1 * y2 - y1 * x2) * (y3 - y4) - (y1 - y2) * (x3 * y4 - y3 * x4)) / denom;
            
            % Check if intersection point lies within both line segments
            if (px >= min(x1, x2) && px <= max(x1, x2) && py >= min(y1, y2) && py <= max(y1, y2)) && ...
               (px >= min(x3, x4) && px <= max(x3, x4) && py >= min(y3, y4) && py <= max(y3, y4))
                % Add intersection point to the list
                intersection_x = [intersection_x; px];
                intersection_y = [intersection_y; py];
            end
        end
    end
    % Debugging: Print intersection points
    fprintf('Debug: Intersection points computed:\n');
    for k = 1:length(intersection_x)
        fprintf('Point %d: X = %.2f, Y = %.2f\n', k, intersection_x(k), intersection_y(k));
    end
end

% % % % 
% % % % function [x, y] = circle_intersection(center1, radius1, center2, radius2)
% % % %     % Compute the intersection points of two circles
% % % %     x1 = center1(1); y1 = center1(2);
% % % %     x2 = center2(1); y2 = center2(2);
% % % % 
% % % %     d = sqrt((x2 - x1)^2 + (y2 - y1)^2); % Distance between centers
% % % % 
% % % %     % Check if circles intersect
% % % %     if d > (radius1 + radius2) || d < abs(radius1 - radius2) || d == 0
% % % %         % No intersection or circles are coincident
% % % %         x = [];
% % % %         y = [];
% % % %         return;
% % % %     end
% % % % 
% % % %     % Compute intersection points
% % % %     a = (radius1^2 - radius2^2 + d^2) / (2 * d);
% % % %     h = sqrt(radius1^2 - a^2);
% % % % 
% % % %     % Point P2 (midpoint between intersection points)
% % % %     x3 = x1 + a * (x2 - x1) / d;
% % % %     y3 = y1 + a * (y2 - y1) / d;
% % % % 
% % % %     % Offset for intersection points
% % % %     offset_x = h * (y2 - y1) / d;
% % % %     offset_y = h * (x2 - x1) / d;
% % % % 
% % % %     % Intersection points
% % % %     x = [x3 + offset_x, x3 - offset_x];
% % % %     y = [y3 - offset_y, y3 + offset_y];
% % % % end

function inside = is_in_cone(sensor, center, point, cone_angle_deg)
    % Returns true if 'point' is inside the cone from 'sensor' toward 'center'
    v_sensor_to_center = center - sensor;
    v_sensor_to_point = point - sensor;
    angle = acos(dot(v_sensor_to_center, v_sensor_to_point) / ...
        (norm(v_sensor_to_center) * norm(v_sensor_to_point)));
    inside = angle <= deg2rad(cone_angle_deg/2);
end

function [category, max_after140, mean_after140] = categorize_signal_after140us(t_rs, signal_env)
    % Categorize the signal after 140us based on thresholds
    t_after140_idx = t_rs > 140e-6;
    env_after140 = signal_env(t_after140_idx);
    
    max_env = max(signal_env);
    mean_after140 = mean(env_after140);
    max_after140 = max(env_after140);

    % Define thresholds (adjust as needed)
    very_low_thresh = 0.05 * max_env;
    low_thresh      = 0.15 * max_env;
    more_thresh     = 0.4 * max_env;

    if max_after140 < very_low_thresh
        category = "center_line";
    elseif max_after140 < low_thresh
        category = "between_center_and_cone";
    elseif max_after140 < more_thresh
        category = "on_cone_line";
    else
        category = "outside_cone";
    end
end


% Helper: function to pick closest intersection to the original loci intersection
function [cx, cy] = pick_closest(center_x, center_y, ref_x, ref_y)
    if isempty(center_x)
        cx = [];
        cy = [];
        return;
    end
    dists = sqrt((center_x - ref_x).^2 + (center_y - ref_y).^2);
    [~, idx] = min(dists);
    cx = center_x(idx);
    cy = center_y(idx);
end


function [roi_x, roi_y] = multilateration_region(sensor_positions, dists, r)
    % Grid search for region closest to all circles
    grid_x = linspace(-r, r, 300);
    grid_y = linspace(-r, r, 300);
    [X, Y] = meshgrid(grid_x, grid_y);
    mask = (X.^2 + Y.^2) <= r^2;
    total_dist = zeros(size(X));
    for i = 1:length(dists)
        D = sqrt((X - sensor_positions(i,1)).^2 + (Y - sensor_positions(i,2)).^2);
        total_dist = total_dist + abs(D - dists(i));
    end
    min_total = min(total_dist(mask));
    roi_mask = (total_dist <= min_total + 0.05*(max(total_dist(mask))-min_total)) & mask;
    roi_x = X(roi_mask);
    roi_y = Y(roi_mask);
end