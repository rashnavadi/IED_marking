

%% Launch first GUI to manually find the peak in the main channel
% By Tahere Rashnavadi
% Last update: June 11, 2025
% NOTES:
%== 
% by adjusted IED times i mean original IED times, they were named adjusted because the original timings marked by the epileptologists were
%done wrt the fMRI runs, but then i adjusted them wrt the continous iEEG time/run
%== 
% manual IED times: are the output IED timings obtained from the first GUI and refer to the ~25 IEDs that have been marked by the user

%% function 1
function manual_ied_loader_gui()
    % === NEW: Get screen size and center the GUI ===
    screenSize = get(0, 'ScreenSize');  % [left bottom width height]
    figWidth = 1000;
    figHeight = 700;
    left = (screenSize(3) - figWidth) / 2;
    bottom = (screenSize(4) - figHeight) / 2;

    % Create the main GUI window (centered and larger)
    fig = uifigure('Name', 'IED Averaging GUI', ...
        'Position', [left bottom figWidth figHeight]);

    % === Grid layout ===
    % Add layout: 5 rows, 2 columns
    glayout = uigridlayout(fig, [5, 2]);
    glayout.RowHeight = {'fit', 'fit', 'fit', 'fit', '1x'};
    glayout.ColumnWidth = {'fit', '1x'};

    % === Left panel containing vertically stacked controls ===
    leftPanel = uipanel(glayout);
    leftPanel.Layout.Row = [1 5];
    leftPanel.Layout.Column = 1;

    % Vertical layout inside left panel
    leftLayout = uigridlayout(leftPanel, [14, 1]);
    leftLayout.RowHeight = repmat({'fit'}, 1, 12);
    leftLayout.ColumnWidth = {'1x'};

    % === File Selection Buttons ===
    
    % EEG file selection button
    eegButton = uibutton(leftLayout, 'Text', 'Select EEG .bin File', ...
        'ButtonPushedFcn', @(btn,event) select_eeg_file(fig));

    % manaul IED file selection button (output file from the previous GUI)
    iedButton = uibutton(leftLayout, 'Text', 'Select IED Manual Timings .txt File', 'ButtonPushedFcn', @(btn,event) select_manual_ied_file(fig));

    % Select Original Adjusted IED File (all IEDs) marked by epileptologists
    selectAdjustedIEDBtn = uibutton(leftLayout, 'Text', 'Select Adjusted/Orig IEDs .txt File', 'ButtonPushedFcn', @(btn,event) load_adjusted_ieds_file(fig));

    % === Processing Buttons ===

    % Button: Align All IEDs Using Template only on main channel
    alignAllButton = uibutton(leftLayout, 'Text', 'Align All IEDs in main Ch (Cross-corr)', 'ButtonPushedFcn', @(btn,event) align_all_ieds(fig));

    % Button: Average All IEDs Using all aligned IEDs after cross correlaiton with the tmeplate25IEDs
    plotAvgAllIEDsBtn = uibutton(leftLayout, 'Text', 'Plot IED Avg of Main Ch (Using All IEDs)', 'ButtonPushedFcn', @(btn, event) plot_Avg_allIEDs_OfMainCh(fig));

    % === Dropdown to choose which IED times to use for averaging ===
    iedSourceDropdown2 = uidropdown(leftLayout, 'Items', {'25 manually marked IEDs', 'All aligned IEDs'}, 'Value', 'All aligned IEDs', ...
        'Tooltip', 'Choose IED timing source for channel-wise average');
    fig.UserData.iedSourceDropdown2 = iedSourceDropdown2;

    plotAvgBtn = uibutton(leftLayout, 'Text', 'Plot IED Avg of All Channels', 'ButtonPushedFcn', @(btn, event) show_avg_all_channels(fig));

    prevButton = uibutton(leftLayout, 'Text', '◀ Prev. EEG Channel', 'ButtonPushedFcn', @(btn, event) go_to_previous_page(fig));

    % Page buttons for EEG channels: Next and Previous buttons in leftLayout
    nextButton = uibutton(leftLayout, 'Text', 'Next EEG Channel ➤', 'ButtonPushedFcn', @(btn, event) go_to_next_page(fig));
   
    % === Onset & Peak Marking ===
    uibutton(leftLayout, 'Text', 'Mark Onset on Avg', 'ButtonPushedFcn', @(btn, event) enable_onset_marking(fig));

    uibutton(leftLayout, 'Text', 'Mark Peak on Avg', 'ButtonPushedFcn', @(btn, event) enable_peak_marking(fig));

    uibutton(leftLayout, 'Text', 'Save Onset & Peak', 'ButtonPushedFcn', @(btn, event) save_onset_peak_times(fig));

    % === Time Window Zoom Dropdown ===
    % Add a horizontal layout for Time Window
    timeWindowLayout = uigridlayout(leftLayout, [1, 2]);
    timeWindowLayout.Layout.Row = 13;
    timeWindowLayout.ColumnWidth = {'fit', 'fit'};
    timeWindowLayout.RowHeight = {'fit'};
    timeWindowLayout.Padding = [0 0 0 0];
    timeWindowLayout.ColumnSpacing = 5;

    % Label
    % Label
    zoomLabel = uilabel(timeWindowLayout, ...
        'Text', 'Time Window:', ...
        'HorizontalAlignment', 'left', ...
        'FontSize', 11);

    % Dropdown
    zoomDropdown = uidropdown(timeWindowLayout, ...
        'Items', {'±0.5s','±1s','±2s','±5s'}, ...
        'Value', '±0.5s', ...
        'Tag', 'TimeWindowMenu', ...
        'FontSize', 11, ...
        'ValueChangedFcn', @(dd, event) changeTimeWindow(fig));

    % Store default window size
    fig.UserData.win_sec = 0.5;

    % === Right Panel for EEG Visualization ===
    plotPanel = uipanel(glayout);
    plotPanel.Layout.Row = [1 5];
    plotPanel.Layout.Column = 2;
    plotPanel.Title = 'EEG traces';
    fig.UserData.plotPanel = plotPanel;

    leftLayout.RowSpacing = 3;
    leftLayout.Padding = [5 5 5 5];

    % Create 5 stacked axes for plotting
    axesHandles = gobjects(1, 5);
    for i = 1:5
        ax = uiaxes(plotPanel);
        ax.Position = [50, 600 - i*110, 800, 100];
        ax.XLabel.String = 'Time (s)';
        ax.YLabel.String = 'Amplitude';
        axesHandles(i) = ax;
    end
    fig.UserData.axesHandles = axesHandles;

    % === Initialize Shared Data ===
    fig.UserData.avgTemplate25 = [];
    fig.UserData.templateAlignedTimes = [];
    fig.UserData.finalAvgTemplate = [];
    fig.UserData.manualIEDTimes = [];
    fig.UserData.isMarkingOnset = false;
    fig.UserData.isMarkingPeak = false;
    fig.UserData.onsetTimes = [];
    fig.UserData.peakTimes = [];
    fig.UserData.t_ied = [];
    fig.UserData.fs = [];
    fig.UserData.channelNames = {};
    fig.UserData.allEEGData = [];
    fig.UserData.currentPage = 1;
    fig.UserData.tracesPerPage = 5;

    fig.UserData.subjectInfoLabel = uilabel(fig, ...
        'Text', 'Subject Info: N/A', ...
        'FontWeight', 'bold', ...
        'FontSize', 14, ...
        'HorizontalAlignment', 'left', ...
        'Position', [300 figHeight - 65 figWidth - 320 30]);

    % === Main channel map ===
    main_map = load('main_channels_map.mat', 'main_channels_map');
    fig.UserData.main_channels_map = main_map.main_channels_map;
    fig.UserData.main_channels = {};
end


%% function 2
function select_eeg_file(fig)
    
    % Step 1: Prompt user to select the EEG binary file
    [bin_file, bin_path] = uigetfile('*.bin', 'Select the EEG binary file (.bin)');

    if isequal(bin_file, 0)
        disp('❌ EEG file selection cancelled.')
        return;
    end

    full_bin_path = fullfile(bin_path, bin_file);
    fig.UserData.eeg_bin_path = full_bin_path;
    fprintf('[✔] EEG file selected: %s\n', full_bin_path);

    % Load EEG binary data
    [eeg_data, fs, channel_labels, ch_map] = load_eeg_bin_with_labels(fig.UserData.eeg_bin_path);
    fig.UserData.eeg_data = eeg_data;
    fig.UserData.fs = fs;
    disp('EEG file loaded successfully.');

    % Step 2: Automatically find the .txt processing log file in same folder
    % Find all .txt files in folder
    txt_candidates = dir(fullfile(bin_path, '*.txt'));

    % Filter only those containing both substrings
    log_candidates = txt_candidates(contains({txt_candidates.name}, 'ProcessingLog') & ...
        contains({txt_candidates.name}, 'ChannelLabels'));

    % Remove hidden/system files (e.g., .DS_Store)
    log_candidates = log_candidates(~startsWith({log_candidates.name}, '.'));

    if isempty(log_candidates)
        uialert(fig, '❌ No processing log file found in the same folder.', 'Missing File');
        return;
    elseif numel(log_candidates) > 1
        warning('[⚠] Multiple log files detected: %s\n', strjoin({log_candidates.name}, ', '));
        % Pick the longest filename or latest modified one
        [~, idx] = max(cellfun(@length, {log_candidates.name}));
        txt_file = log_candidates(idx).name;
        fprintf('[⚠] Multiple logs — defaulting to: %s\n', txt_file);
    else
        txt_file = log_candidates(1).name;
    end

    full_log_path = fullfile(bin_path, txt_file);
    fprintf('[✔] Automatically loaded processing log: %s\n', full_log_path);

    % Step 3: Read file line-by-line, Extract sampling rate, number of channels, and labels
    fid = fopen(full_log_path, 'r');
    lines = {};
    while ~feof(fid)
        lines{end+1} = fgetl(fid); %#ok<AGROW>
    end
    fclose(fid);

    s_Rate = NaN;
    n_Channels = NaN;
    c_Labels = {};

    for i = 1:length(lines)
        line = strtrim(lines{i});
        if contains(line, 's_Rate:')
            s_Rate = sscanf(line, 's_Rate: %f');
        elseif contains(line, 'n_Channels:')
            n_Channels = sscanf(line, 'n_Channels: %f');
        elseif contains(line, 'c_Labels (Row):')
            % Read all remaining channel labels on that line and possibly the next line(s)
            label_line = strrep(line, 'c_Labels (Row):', '');
            labels_str = strtrim(label_line);
            % If more labels are in the next line(s), keep reading
            j = i + 1;
            while j <= length(lines) && ~startsWith(strtrim(lines{j}), 'c_Labels (Column):')
                labels_str = [labels_str ' ' strtrim(lines{j})];
                j = j + 1;
            end
            c_Labels = strsplit(strtrim(labels_str));
        end
    end

    if isnan(s_Rate) || isnan(n_Channels) || isempty(c_Labels)
        uialert(fig, '❌ Could not extract all required info from log file.', 'Error');
        return;
    end

    % Step 4: Build bipolar channel pairs
    bipolar_labels = cell(1, length(c_Labels) - 1);
    for i = 1:length(c_Labels) - 1
        bipolar_labels{i} = [c_Labels{i} '-' c_Labels{i+1}];
    end

    % Step 5: Save into UserData
    fig.UserData.fs = s_Rate;
    fig.UserData.nchannels = n_Channels;
    fig.UserData.channelnames_mono = c_Labels;
    fig.UserData.channelnames_bipolar = bipolar_labels;
    fig.UserData.totalChannels = length(bipolar_labels);
    fig.UserData.currentPage = 1;

    fig.UserData.eegStruct = struct( ...
        'data', eeg_data, ...
        'path', full_bin_path, ...
        'fs', s_Rate, ...
        'labels', c_Labels ...
        );

    fprintf('Loaded sampling rate: %d Hz\n', s_Rate);
    fprintf('Found %d monopolar channels, %d bipolar pairs.\n', ...
        numel(c_Labels), numel(bipolar_labels));

    % Placeholder: confirmation for now
    disp(['Extracted ', num2str(numel(bipolar_labels)), ' bipolar traces. Sampling Rate: ', num2str(s_Rate), ' Hz']);
end

%% function 3
function select_manual_ied_file(fig)

    disp('Please select the manual IED timing file (.txt)');

    % Step 1: Let user select .txt file
    [ied_file, ied_path] = uigetfile('*_manual.txt', 'Select manual IED timing file');
    if isequal(ied_file, 0)
        disp('[⚠] Manual IED file selection cancelled.');
        return;
    end

    ied_txt_path = fullfile(ied_path, ied_file);
    fig.UserData.ied_txt_path = ied_txt_path;

    % Step 2: Parse subject/run/IED from filename
    [~, ied_filename_noext, ~] = fileparts(ied_file);
    tokens = regexp(ied_filename_noext, '(ICE\d+)_([a-zA-Z0-9]+)_([a-zA-Z0-9]+)_manual', 'tokens');

    if isempty(tokens)
        warning('❌ Could not parse subject/run/IED info from filename: %s', ied_file);
        return;
    end
    tokens = tokens{1};
    txt_subject = tokens{1};
    txt_run     = tokens{2};
    txt_ied     = tokens{3};

    disp(['✔ Loaded manual IED timing file: ', ied_file]);

    % Extract numeric run number from EEG file path
    [~, eeg_filename_noext, ~] = fileparts(fig.UserData.eeg_bin_path);
    eeg_run_match = regexp(eeg_filename_noext, 'Run(\d+)', 'tokens');
    if isempty(eeg_run_match)
        warning('Could not extract run number from EEG file name.');
    else
        eeg_run_num = str2double(eeg_run_match{1}{1});  % e.g., 6 from Run6
        % Now extract numeric part from txt_run (e.g., Run3b → 3)
        txt_run_match = regexp(txt_run, 'Run(\d+)', 'tokens');
        if isempty(txt_run_match)
            uialert(fig, sprintf('Run mismatch: IED file run "%s" is not valid.', txt_run), 'Run Mismatch');
            return;
        end
        txt_run_num = str2double(txt_run_match{1}{1});

        if eeg_run_num ~= txt_run_num
            msg = sprintf(['❌ Run number mismatch:\n' ...
                'EEG file run: Run%d\n' ...
                'IED txt file run: %s\n\n' ...
                'Please load the correct IED file for Run%d.'], ...
                eeg_run_num, txt_run, eeg_run_num);
            uialert(fig, msg, 'Run Mismatch');
            return;
        end
    end

    % Step 3: Read all lines (incl. "skipped by user")
    fid = fopen(ied_txt_path, 'r');
    lines = {};
    marking_channel = 'Unknown';

    line_num = 0;
    while ~feof(fid)
        line = strtrim(fgetl(fid));
        line_num = line_num + 1;

        if line_num == 1 && contains(line, 'Channel used for IED timing adjustment')
            tokens = regexp(line, 'adjustment:\s*(\S+)', 'tokens');
            if ~isempty(tokens)
                marking_channel = tokens{1}{1};
            end
            continue;  % Skip channel header
        end

        if ~startsWith(line, '#') && ~isempty(line)
            lines{end+1} = line; %#ok<AGROW>
        end
    end
    fclose(fid);

    % Store and display marking channel
    fig.UserData.marking_channel_used = marking_channel;
    fig.UserData.markingChannelLabel.Text = sprintf('📍 Marked in: %s', marking_channel);

    % Convert to numeric
    n_ieds = numel(lines);
    manual_ied_times = nan(n_ieds, 1);  % NaN for skipped
    skip_flags = false(n_ieds, 1);      % true for skipped

    for i = 1:n_ieds
        txt = strtrim(lines{i});
        if contains(lower(txt), 'skipped')
            skip_flags(i) = true;
        else
            nums = sscanf(txt, '%f');  % Extract all numeric values in line
            if numel(nums) >= 2
                manual_ied_times(i) = nums(2);  % Use the second number (actual IED time)
            elseif numel(nums) == 1
                manual_ied_times(i) = nums(1);  % If only one number (no line index), use it
            else
                fprintf('[⚠️] Line %d could not be parsed as number: "%s"\n', i, txt);
            end
        end
    end

    fprintf('[DEBUG] Parsed %d valid IED times, %d skipped.\n', ...
        sum(~isnan(manual_ied_times)), sum(isnan(manual_ied_times)));

    % Validate we have at least 1 non-skipped
    if all(skip_flags)
        uialert(fig, '❌ All IEDs marked as skipped. Cannot proceed.', 'Empty File');
        return;
    end

    % Step 4: Save to UserData
    fig.UserData.manualIEDTimes = manual_ied_times(:);  % NaNs for skipped
    fig.UserData.manualSkipFlags = skip_flags(:);
    fig.UserData.manualIED_samples = round(manual_ied_times * fig.UserData.fs);

    % to avoid cross-subject contamination:
    fig.UserData.avgTemplate25 = [];
    fig.UserData.templateAlignedTimes = [];
    fig.UserData.finalAvgTemplate = [];
    fig.UserData.onsetTimes = [];
    fig.UserData.peakTimes = [];

    % Step 5: Update metadata
    fig.UserData.subject_id = txt_subject;
    fig.UserData.run_id     = txt_run;
    fig.UserData.ied_id     = txt_ied;

    % Try to load previous onset/peak if exists
    [ied_folder, ~, ~] = fileparts(fig.UserData.ied_txt_path);
    [parent_dir, ~, ~] = fileparts(ied_folder);  % go up to subject folder
% 
    save_dir = fullfile(parent_dir, 'avg_onset_peak_times');
    basename = sprintf('%s_%s_%s', fig.UserData.subject_id, fig.UserData.run_id, fig.UserData.ied_id);
    mat_path = fullfile(save_dir, [basename '_onset_peak_times.mat']);
    
    % Step 6: Load main channels for optional highlighting
    key = sprintf('%s_%s', txt_subject, txt_ied);
    if isKey(fig.UserData.main_channels_map, key)
        fig.UserData.main_channels = fig.UserData.main_channels_map(key);
        fprintf('Main channels for %s loaded.\n', key);
    else
        fig.UserData.main_channels = {};  % fallback
        warning('No main channels found for %s', key);
    end

    % Step 7: Validate EEG is already loaded
    if isfield(fig.UserData, 'eeg_data') && isfield(fig.UserData, 'channelnames_bipolar')
        fig.UserData.totalChannels = length(fig.UserData.channelnames_bipolar);
        if ~isfield(fig.UserData, 'onsetTimes') || isempty(fig.UserData.onsetTimes)
            fig.UserData.onsetTimes = nan(1, fig.UserData.totalChannels);
        end
        if ~isfield(fig.UserData, 'peakTimes') || isempty(fig.UserData.peakTimes)
            fig.UserData.peakTimes  = nan(1, fig.UserData.totalChannels);
        end
    else
        uialert(fig, 'Missing EEG data or channel labels. Load EEG first.', 'Missing Data');
    end

    % Step 8: Update label
    update_subject_info_label(fig);

    % Final message
    fprintf('[✔] Loaded %d valid IED timings from: %s\n', length(manual_ied_times), ied_file);

    % Step 9: Auto-compute average IED template and plot
    if length(fig.UserData.manualIEDTimes) >= 1
        compute_avg_template(fig);
    end

%     % === Step 10: Attempt to auto-load saved onset/peak if it exists ===
%     % Build expected filename
    parent_dir = fileparts(ied_txt_path);  % This is IED_Cleaned
    fig.UserData.parent_dir = parent_dir;
    save_dir = fullfile(parent_dir, '..', 'avg_onset_peak_times');
    expected_txt = fullfile(save_dir, sprintf('%s_%s_%s_avg_onset_peak_times.txt', txt_subject, txt_run, txt_ied));

    % === Set output_dir for saving aligned results like final aligned IEDs ===
    subject_dir = fileparts(parent_dir);  % go up from IED_Cleaned to subject folder
    output_dir = fullfile(subject_dir, 'avg_onset_peak_times');
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end
    fig.UserData.output_dir = output_dir;
    fprintf('[📁] Output directory set to: %s\n', output_dir);

end

%% function 4
% load_adjusted_ieds_file this is the original/adjusted IED timings 
function load_adjusted_ieds_file(fig)

    [fname, fpath] = uigetfile('*_adjusted.txt', 'Select adjusted IED timing file (all detections)');
    if isequal(fname, 0)
        return;
    end

    full_path = fullfile(fpath, fname);

    % Step 1: Extract subject, run, and IED type from filename
    [~, basename, ~] = fileparts(fname);
    tokens = regexp(basename, '(ICE\d+)_([a-zA-Z0-9]+)_([a-zA-Z0-9]+)_adjusted', 'tokens');

    if isempty(tokens)
        uialert(fig, 'Filename does not match expected format: ICE###_RunX_IEDY_adjusted.txt', 'Invalid Format');
        return;
    end

    tokens = tokens{1};
    adj_subject = tokens{1};
    adj_run     = tokens{2};
    adj_ied     = tokens{3};

    % Step 2: If manual file was already loaded, check match
    if isfield(fig.UserData, 'subject_id')
        man_subject = fig.UserData.subject_id;
        man_run     = fig.UserData.run_id;
        man_ied     = fig.UserData.ied_id;

        if ~strcmp(man_subject, adj_subject) || ...
           ~strcmp(man_run, adj_run) || ...
           ~strcmp(man_ied, adj_ied)

            msg = sprintf(['❌ Subject/run/IED mismatch:\n\n' ...
                'Manual file:   %s %s %s\nAdjusted file: %s %s %s\n\n' ...
                'Please load matching adjusted file.'], ...
                man_subject, man_run, man_ied, ...
                adj_subject, adj_run, adj_ied);
            uialert(fig, msg, 'File Mismatch');
            return;
        end
    else
        % If manual file not loaded yet, store metadata from adjusted file
        fig.UserData.subject_id = adj_subject;
        fig.UserData.run_id     = adj_run;
        fig.UserData.ied_id     = adj_ied;
    end

    % Step 3: Load adjusted IED times
    all_ied_times = load(full_path);
    fig.UserData.adjustedIEDTimes = all_ied_times(:);

    fprintf('[📥] Loaded %d IED times from adjusted file: %s\n', ...
        length(all_ied_times), fname);
end

%% function 5 (correct)
% compute the average IED template from the main channel only over the 25 manually marked IEDs
% only use the main channel (the one noted in the first line of the manual IED .txt file) to compute the average template
% this func gives single-channel template for aligning IEDs across all channels

function compute_avg_template(fig)
    % Compute average IED template only from the manually marked main channel
    fs = fig.UserData.fs;
    eeg = fig.UserData.eeg_data;
    manual_IED_times_sec = fig.UserData.manualIEDTimes;

    % === Step 1: Filter valid IEDs ===
    % 🧠 Filter out "skipped by user" entries (NaNs)
    valid_idx = ~isnan(manual_IED_times_sec);
    manual_IED_times_sec = manual_IED_times_sec(valid_idx);  % Keep only valid manual IEDs
    fprintf('[DEBUG] Using %d of %d manual IEDs (non-skipped).\n', ...
        numel(manual_IED_times_sec), numel(fig.UserData.manualIEDTimes));

    if isempty(manual_IED_times_sec)
        error('No valid IEDs found to compute average template.');
    end

    % === Step 2: Parse main channel ===
    % Get the marking channel used (from first line of manual file)
    target_label = strtrim(fig.UserData.marking_channel_used);    
    channelnames = fig.UserData.channelnames_bipolar;
    mono_map = containers.Map(fig.UserData.channelnames_mono, 1:numel(fig.UserData.channelnames_mono));

    % Validate marking channel exists
    ch_idx = find(strcmp(channelnames, target_label), 1);
    if isempty(ch_idx)
        error('Main marking channel (%s) not found in channel list.', target_label);
    end

    % Get mono indices for that channel
    parts = split(target_label, '-');
    ch1 = parts{1}; ch2 = parts{2};

    if ~isKey(mono_map, ch1) || ~isKey(mono_map, ch2)
        error('Main marking channel parts (%s, %s) not found in mono map.', ch1, ch2);
    end
    idx1 = mono_map(ch1);
    idx2 = mono_map(ch2);

    % === Step 3: Segment extraction parameters ===
    pre_samples = round(0.5 * fs);
    post_samples = round(0.5 * fs);
    win_len = pre_samples + post_samples + 1;
    n_ieds = numel(manual_IED_times_sec);
    % Convert adjusted times to sample indices
    sample_times = round(manual_IED_times_sec * fs);

    % === Step 4: Collect segments and compute average ===
    % Extract IED segments for that channel
    segments = zeros(n_ieds, win_len);
    count = 0;
    for i = 1:n_ieds
        center = sample_times(i);
        start_idx = center - pre_samples;
        end_idx   = center + post_samples;

        if start_idx < 1 || end_idx > size(eeg, 2)
            fprintf('[SKIP] IED %d: window [%d–%d] is out of bounds.\n', i, start_idx, end_idx);
            continue;
        end
        segments(count + 1, :) = eeg(idx1, start_idx:end_idx) - eeg(idx2, start_idx:end_idx);
        count = count + 1;
    end

    if count == 0
        error('No valid IED segments extracted for averaging.');
    end

    avg_template = mean(segments(1:count, :), 1);
    fig.UserData.avgTemplate25 = avg_template;
    fprintf('[✔] Avg template computed from %d IEDs on channel %s.\n', count, target_label);

    fig.UserData.currentTemplateKey = 'avgTemplate25';  % ✅ Track current template for editing

    % === Step 5: Plot the averagtemplat overe 25 IEDs in a separate figure ===
    t = (-pre_samples:post_samples) / fs;
    fig_avg = figure('Name', 'Avg IED Template', 'Color', 'w');
    plot(t, avg_template, 'k', 'LineWidth', 1.5);
    xlabel('Time (s)');
    ylabel('Amplitude');
    title(sprintf('Avg IED Template | Subject: %s | Run: %s | IED: %s | Channel: %s', ...
        fig.UserData.subject_id, fig.UserData.run_id, fig.UserData.ied_id, target_label), ...
        'Interpreter', 'none');
    grid on;
    xlim([-0.5, 0.5]);

    % === Step 6: Save figure and raw data ===
    [ied_folder, ~, ~] = fileparts(fig.UserData.ied_txt_path);
    [parent_dir, ~, ~] = fileparts(ied_folder);
    save_dir = fullfile(parent_dir, 'avg_onset_peak_times');
    if ~exist(save_dir, 'dir'), mkdir(save_dir); end

    base = sprintf('%s_%s_%s', fig.UserData.subject_id, fig.UserData.run_id, fig.UserData.ied_id);

    saveas(fig_avg, fullfile(save_dir, [base '_avg_template25.png']));
    save(fullfile(save_dir, [base '_avg_template25.mat']), 'avg_template', 't');

    % Inform user
    msg = sprintf(['✅ The average IED template was successfully computed from %d manually marked IEDs on channel %s.\n\n' ...
        'The figure and data have been saved to:\n\n%s'], count, target_label, save_dir);

    uialert(fig, msg, 'Avg Template Saved');
end

%% 
function plot_Avg_allIEDs_OfMainCh(fig)
    avg_allIEDs_OfMainCh(fig);  % computes the average IED waveform for the main marking channel only, based on all aligned IEDs

    % Extract
    avg_waveform = fig.UserData.AvgAllIEDs_OfMainCh;
    fs = fig.UserData.fs;
    win_sec = fig.UserData.win_sec;
    t = linspace(-win_sec, win_sec, numel(avg_waveform));

    % Create new figure and make sure it stays open
    f = figure('Name', 'Main Channel Avg (All IEDs)', 'Color', 'white', 'Visible', 'on');  % <== explicitly make it visible
    fig.UserData.mainAvgFig = f;  

    ax = axes(f);
    plot(ax, t, avg_waveform, 'k-', 'LineWidth', 2);
    title(ax, sprintf('Main Channel Avg (All Aligned IEDs)\n%s', strtrim(fig.UserData.marking_channel_used)), 'Interpreter', 'none');
    xlabel(ax, 'Time (s)');
    ylabel(ax, 'Amplitude');
    grid(ax, 'on');
    xlim(ax, [-win_sec win_sec]);

    % === Save figure and waveform ===
    [ied_folder, ~, ~] = fileparts(fig.UserData.ied_txt_path);
    [parent_dir, ~, ~] = fileparts(ied_folder);
    save_dir = fullfile(parent_dir, 'avg_onset_peak_times');
    if ~exist(save_dir, 'dir'), mkdir(save_dir); end

    base = sprintf('%s_%s_%s', fig.UserData.subject_id, fig.UserData.run_id, fig.UserData.ied_id);

    % Save figure
    saveas(f, fullfile(save_dir, [base '_mainChanAvg_allIEDs.png']));

    % Save waveform + time vector
    save(fullfile(save_dir, [base '_mainChanAvg_allIEDs.mat']), 'avg_waveform', 't');

    % Optional alert
    msg = sprintf('✅ The average waveform over ALL aligned IEDs on the main channel was saved to:\n\n%s', save_dir);
    uialert(fig, msg, 'All IEDs Avg Saved');

end


%% computes the average IED waveform for the main marking channel only, based on all aligned IEDs
function avg_allIEDs_OfMainCh(fig)
    eeg = fig.UserData.eeg_data;
    fs = fig.UserData.fs;
    aligned_samples = fig.UserData.templateAlignedTimes;  % [nIEDs x 1], same across channels, This is not in seconds

    main_label = strtrim(fig.UserData.marking_channel_used);  % e.g., 'dRaIN1-dRaIN3'
    channelnames = fig.UserData.channelnames_bipolar;
    mono_map = containers.Map(fig.UserData.channelnames_mono, 1:numel(fig.UserData.channelnames_mono));
    
    % Get index of main channel
    parts = split(main_label, '-');
    ch1 = parts{1}; ch2 = parts{2};

    if ~isKey(mono_map, ch1) || ~isKey(mono_map, ch2)
        warning('[⚠] Main channel %s contains invalid mono names.', main_label);
        return;
    end
    idx1 = mono_map(ch1);
    idx2 = mono_map(ch2);

    pre_samples = round(0.5 * fs);
    post_samples = round(0.5 * fs);
    win_len = pre_samples + post_samples + 1;

    segments = [];
    for i = 1:numel(aligned_samples)
        center = aligned_samples(i);
        start_idx = center - pre_samples;
        end_idx = center + post_samples;
        if start_idx >= 1 && end_idx <= size(eeg, 2)
            segment = eeg(idx1, start_idx:end_idx) - eeg(idx2, start_idx:end_idx);
            segments(end+1, :) = segment; %#ok<AGROW>
        end
    end

    if ~isempty(segments)
        fig.UserData.AvgAllIEDs_OfMainCh = mean(segments, 1);  % 1×T row vector
        fprintf('[✔] Main channel average from all aligned IEDs computed.\n');
    else
        warning('[⚠] No valid IED segments found for main channel.');
        fig.UserData.AvgAllIEDs_OfMainCh = NaN(1, win_len);
    end
end

%% This function computes the average IED waveform for each bipolar channel, using either:
% the 25 manually marked IEDs (manualIED_samples), or
% the IED-aligned times after template matching (templateAlignedTimes) — depending on the dropdown value.
function compute_avg_all_channels(fig)
    fs = fig.UserData.fs;
    eeg = fig.UserData.eeg_data;
    channelnames = fig.UserData.channelnames_bipolar;
    mono_map = containers.Map(fig.UserData.channelnames_mono, 1:numel(fig.UserData.channelnames_mono));

    % Read from dropdown choice
    source = fig.UserData.iedSourceDropdown2.Value;
    if strcmp(source, '25 manually marked IEDs')
        aligned_samples = fig.UserData.manualIED_samples;
    else
        aligned_samples = fig.UserData.templateAlignedTimes; % in samples not sec
    end

    win_sec = fig.UserData.win_sec;
    pre_samples = round(win_sec * fs);
    post_samples = round(win_sec * fs);

    win_len = pre_samples + post_samples + 1;
    n_channels = numel(channelnames);
    avg_matrix = nan(n_channels, win_len);

    for ch = 1:n_channels
        label = channelnames{ch};
        parts = split(label, '-');
        ch1 = parts{1}; ch2 = parts{2};

        if ~isKey(mono_map, ch1) || ~isKey(mono_map, ch2)
            continue;
        end
        idx1 = mono_map(ch1);
        idx2 = mono_map(ch2);

        segments = [];
        for i = 1:numel(aligned_samples)
            center = aligned_samples(i);
            start_idx = center - pre_samples;
            end_idx = center + post_samples;
            if start_idx >= 1 && end_idx <= size(eeg, 2)
                segment = eeg(idx1, start_idx:end_idx) - eeg(idx2, start_idx:end_idx);
                segments(end+1, :) = segment; %#ok<AGROW>
            end
        end

        if ~isempty(segments)
            avg_matrix(ch, :) = mean(segments, 1);
        end
    end
    fig.UserData.avgAllChannelsComputed = true;
    fig.UserData.avgAllChannels = avg_matrix;
    fig.UserData.avgAllChannelsSource = source;
    fig.UserData.avg_t = (-pre_samples:post_samples) / fs;
    plot_avg_template(fig, 'avgAllChannels');
    fig.UserData.currentTemplateKey = 'avgAllChannels';  % NEW LINE

    fig.UserData.finalAvgTemplate = avg_matrix;  % Set final template alias
    save_final_avg_template(fig);  % Reuse the save function


    fprintf('[✔] Computed average EEG for all channels using "%s".\n', source);
end

%%
function show_avg_all_channels(fig)
    compute_avg_all_channels(fig);  % Compute based on dropdown choice
    plot_avg_template(fig, 'avgAllChannels');
    fig.UserData.currentTemplateKey = 'avgAllChannels';
end

%% function 6 
% This is a plotting/visualization function only, not the alignment function
% the main channel
function plot_avg_template(fig, template_key)
    % template_key: 'avgTemplate25' or 'finalAvgTemplate'
    
    fprintf('[DEBUG] Checking for template key: %s\n', template_key);
    if isfield(fig.UserData, template_key)
        fprintf('[DEBUG] Field "%s" exists in fig.UserData.\n', template_key);
    else
        fprintf('[DEBUG] Field "%s" is MISSING in fig.UserData.\n', template_key);
    end
    
    if isfield(fig.UserData, template_key)
        fprintf('[DEBUG] Size of %s: [%s]\n', template_key, mat2str(size(fig.UserData.(template_key))));
    end

     
    if ~isfield(fig.UserData, template_key) || isempty(fig.UserData.(template_key))
        uialert(fig, 'Template not found. Compute it first.', 'Missing Template');
        return;
    end

    avg_data = fig.UserData.(template_key);  % [nChannels x nTimePoints] OR [1 x T]
    fs = fig.UserData.fs;

    % DEBUG info
    fprintf('[DEBUG] Plotting template: %s | Size: [%d x %d]\n', ...
        template_key, size(avg_data, 1), size(avg_data, 2));

    % If this is avgTemplate25 and it's only 1 row, force it to show only for main channel
    if strcmp(template_key, 'avgTemplate25') && size(avg_data, 1) == 1
        main_label = strtrim(fig.UserData.marking_channel_used);
        main_ch_idx = find(strcmp(fig.UserData.channelnames_bipolar, main_label));
        full_data = nan(length(fig.UserData.channelnames_bipolar), size(avg_data, 2));
        if ~isempty(main_ch_idx)
            full_data(main_ch_idx, :) = avg_data;
        end
        avg_data = full_data;
    end

    if strcmp(template_key, 'mainAvgAllIEDs') && size(avg_data, 1) == 1
        main_label = strtrim(fig.UserData.marking_channel_used);
        main_ch_idx = find(strcmp(fig.UserData.channelnames_bipolar, main_label));
        full_data = nan(length(fig.UserData.channelnames_bipolar), size(avg_data, 2));
        if ~isempty(main_ch_idx)
            full_data(main_ch_idx, :) = avg_data;
        end
        avg_data = full_data;
    end

    nChannels = size(avg_data, 1);
    nTime = size(avg_data, 2);
    if isfield(fig.UserData, 'avg_t')
        t = fig.UserData.avg_t;
    else
        win_sec = fig.UserData.win_sec;
        t = linspace(-win_sec, win_sec, nTime);
    end

    % Paging setup
    perPage = fig.UserData.tracesPerPage;
    page = fig.UserData.currentPage;
    startIdx = (page - 1) * perPage + 1;
    endIdx = min(startIdx + perPage - 1, nChannels);
    visibleIdx = startIdx:endIdx;

    axesHandles = fig.UserData.axesHandles;

    % Clear axes
    for ax = axesHandles
        cla(ax);
        ax.Title.String = '';
        ax.Color = 'white';
        ax.XColor = 'black';
        ax.YColor = 'black';
        ax.LineWidth = 0.5;
    end

    % PLOTTING ==== Compute global safe Y-limit based on all visible channels ====
    visible_traces = avg_data(visibleIdx, :);
    global_min = min(visible_traces(:));
    global_max = max(visible_traces(:));
    yrange = global_max - global_min;
    if yrange == 0
        pad = 1;  % minimum padding if all values are the same
    else
        pad = 0.1 * yrange;
    end
    global_ylim = [global_min - pad, global_max + pad];

    % Plot traces
    for i = 1:length(axesHandles)
        ax = axesHandles(i);

        if i <= length(visibleIdx) && visibleIdx(i) <= size(avg_data, 1)
            chIdx = visibleIdx(i);
            ch_name = fig.UserData.channelnames_bipolar{chIdx};
            fprintf('[PLOT] Plotting channel %d (%s)...\n', chIdx, ch_name);
            trace = avg_data(chIdx, :);

            % Skip bad channels
            if all(isnan(trace)) || all(trace == 0) || std(trace) < 1e-6
                warning('[⚠️] Channel %d (%s) has a flat or invalid trace. Skipping plot.', chIdx, ch_name);
                cla(ax);
                ax.Title.String = sprintf('%s (ch %d of %d)', ch_name, chIdx, nChannels);
                xlim(ax, [t(1), t(end)]);
                ylim(ax, [-1 1]);
                continue;
            end

            % Check if it's the main channel
            main_label = strtrim(fig.UserData.marking_channel_used);
            is_main_channel = strcmp(ch_name, main_label);

            % Set trace color and background
            if is_main_channel
                trace_color = [0 0 0];           % black
                ax.Color = [0.9 1.0 0.9];        % light green
            else
                trace_color = [0.7 0.7 0.7];     % gray
                ax.Color = 'white';
            end

            % Plot the trace
            safe_ylim = global_ylim;
            lineHandle = plot(ax, t, trace, 'Color', trace_color, 'LineWidth', 2);
            lineHandle.HitTest = 'off';

            % Show previously marked onset/peak
            if ~isnan(fig.UserData.onsetTimes(chIdx))
                xline(ax, fig.UserData.onsetTimes(chIdx), '--', 'Color', [0 0.6 0], 'LineWidth', 2);
            end
            if ~isnan(fig.UserData.peakTimes(chIdx))
                xline(ax, fig.UserData.peakTimes(chIdx), '--', 'Color', [0.5 0 0.8], 'LineWidth', 2);
            end

            % Axes settings
            ax.ButtonDownFcn = @(src, event) handle_click_on_avg(fig, src, event);
            ax.HitTest = 'on';
            ax.PickableParts = 'all';
            xlim(ax, [t(1), t(end)]);
            if numel(safe_ylim) == 2 && all(isfinite(safe_ylim)) && diff(safe_ylim) > 0
                ylim(ax, safe_ylim);
            else
                ylim(ax, [-1 1]);
                warning('[⚠️] Invalid Y-limits detected for channel %d. Using fallback [-1 1].', chIdx);
            end
            title(ax, sprintf('%s (ch %d of %d)', ch_name, chIdx, nChannels));
            xlabel(ax, 'Time (s)');
            ylabel(ax, 'Amplitude');
            grid(ax, 'on');

            % Save subplot image (optional)
            if isfield(fig.UserData, 'parent_dir')
                save_dir = fullfile(fig.UserData.parent_dir, 'avg_IED_subplots');
                if ~exist(save_dir, 'dir')
                    mkdir(save_dir);
                end
                fig_filename = sprintf('%s_%s_%s_ch%02d_%s.png', ...
                    fig.UserData.subject_id, fig.UserData.run_id, fig.UserData.ied_id, ...
                    chIdx, ch_name);
            end

        else
            % Handle unused axes
            cla(ax);
            chIdx = (fig.UserData.currentPage - 1) * fig.UserData.tracesPerPage + i;
            ax.Title.String = sprintf('Channel %d (empty)', chIdx);
            ax.XLim = [t(1), t(end)];
            if numel(global_ylim) == 2 && all(isfinite(global_ylim)) && diff(global_ylim) > 0
                ax.YLim = global_ylim;
            else
                ax.YLim = [-1 1];
                warning('[⚠️] Invalid Y-limits in EMPTY axis (%d). Using fallback [-1 1].', i);
            end
        end
    end


    % Update label
    drawnow;
    update_subject_info_label(fig);

    % Enable click detection only if user is in marking mode
    if fig.UserData.isMarkingOnset || fig.UserData.isMarkingPeak
    else
        set(fig, 'WindowButtonDownFcn', []);  % optional: clear when not marking
    end
end

%% function 7 (correct)
% To make all IEDs consistently aligned in time (mainchannel) before averaging
% input: - avgTemplate25 (template), and adjusted/Original IED times from
% the rest of IEDs that have not been manually marked in the firs GUI
% ✅ Output:
% fig.UserData.templateAlignedTimes: [nIEDs x 1] aligned sample indices
% Saved .txt file: *_aligned_mainchan.txt with aligned IED times (in seconds)
% Manual IEDs are left untouched
% Remaining IEDs are shifted via cross-correlation (main channel only)

function align_all_ieds(fig)

    if isempty(fig.UserData.avgTemplate25)
        uialert(fig, 'Compute the initial average template first.', 'Missing Template');
        return;
    end

    wait = waitbar(0, 'Aligning all IEDs... Please wait.');

    try
        % Align all original IEDs using cross-correlation with avgTemplate25
        fs = fig.UserData.fs;
        eeg = fig.UserData.eeg_data;
        avg_template = fig.UserData.avgTemplate25;
        manual_times = fig.UserData.manualIEDTimes;     % 25 manually marked IEDs
        adjusted_times = fig.UserData.adjustedIEDTimes; % all detected IEDs (uploaded via second button)

        % Combine both sets of IEDs, avoiding duplicates
        skip_flags = fig.UserData.manualSkipFlags;

        if length(manual_times) ~= length(adjusted_times)
            error('Manual and adjusted IED timing files must have the same number of entries.');
        end

        % === Get the main marking channel ===
        channelnames = fig.UserData.channelnames_bipolar;
        mono_map = containers.Map(fig.UserData.channelnames_mono, 1:numel(fig.UserData.channelnames_mono));
        main_label = strtrim(fig.UserData.marking_channel_used);

        parts = split(main_label, '-');
        ch1 = parts{1}; ch2 = parts{2};
        if ~isKey(mono_map, ch1) || ~isKey(mono_map, ch2)
            error('Main marking channel parts not found.');
        end
        idx1 = mono_map(ch1);
        idx2 = mono_map(ch2);

        % Parameters
        n_ieds = numel(fig.UserData.adjustedIEDTimes);  % ✅ adjusted/orig IED file
        template_len = length(avg_template);
        
        % aligned_times: single column vector for main channel only, 
        % A vector of absolute IED times in seconds for all IEDs (both manual and cross-correlation adjusted).
        aligned_times = nan(n_ieds, 1);  

        pre_samples = round(0.5 * fs);    % 1250 samples before the center
        post_samples = round(0.5 * fs);   % 1250 samples after the center
        search_margin = round(0.2 * fs);  % 500 samples (±200 ms): search ±200 msec around center

        % === Main alignment loop
        for i = 1:n_ieds
            if skip_flags(i)
                % This IED was not manually marked — align it
                center = round(adjusted_times(i) * fs); % center = adjusted/original IED times
                % Define search window for segment of EEG is considered for cross-correlation
                start_idx = center - search_margin - pre_samples;  % start_idx = center - 500 - 1250 = center - 1750
                end_idx   = center + search_margin + post_samples; % end_idx   = center + 500 + 1250 = center + 1750
                % full_segment = eeg(idx1, center - 1750 : center + 1750);
                % Length = center + 1750 - (center - 1750) + 1 = 3501 samples
                % At 2500 Hz → 3501 samples/ 2500 = 1.4004 seconds
                % This is centered around the IED time ±0.7 sec
                % You are preparing a large enough segment (±0.7 s) to: Later slide the template (which is 1.0 second long, 2500 samples) across a centered ±200 ms window (±500 samples)
                % That’s why you need 0.7 s on both sides of the center to allow for this sliding.

                if start_idx < 1 || end_idx > size(eeg, 2)
                    continue;
                end

                full_segment = eeg(idx1, start_idx:end_idx) - eeg(idx2, start_idx:end_idx);

                % Slide template across center ± search_margin
                max_corr = -inf;
                best_lag = 0;

                for lag = -search_margin:search_margin % ±500 samples (±200 ms):  all 401 lags for alignment, from -200 ms to +200 ms.
                    seg_start = search_margin + 1 + lag;
                    seg_end   = seg_start + template_len - 1;

                    if seg_start < 1 || seg_end > length(full_segment)
                        continue;
                    end

                    window = full_segment(seg_start:seg_end);
                    [cross_corr, lags] = xcov(window, avg_template, 0, 'coeff'); % This slides the template25 over ±100 ms relative to the center.
                    [c, max_idx] = max(cross_corr);
                    best_lag = lags(max_idx);    % in samples

                    if ~isnan(c) && c > max_corr
                        max_corr = c;
                        best_lag = lag;
                    end
                end

                % Update aligned time (shifted from adjusted time)
                aligned_times(i) = (center + best_lag) / fs;  % convert back to seconds
            else
                % Keep the manual marking as-is
                aligned_times(i) = manual_times(i);
            end

            waitbar(i / n_ieds, wait, sprintf('Aligning: %d/%d IEDs...', i, n_ieds));
        end

        % === Store result
        fig.UserData.templateAlignedTimes = round(aligned_times * fs); % in samples
        fig.UserData.templateAlignedTimes_sec = aligned_times;
        fprintf('[✔] Stored %d aligned IEDs × %d channels in templateAlignedTimes.\n', size(aligned_times, 1), size(aligned_times, 2));
        fprintf('[ℹ️] Includes %d manually marked IEDs and %d cross-correlation aligned IEDs.\n', sum(~fig.UserData.manualSkipFlags), sum(fig.UserData.manualSkipFlags));

        aligned_secs = aligned_times / fs; % all IEDs after alignments using cross-corr (in seconds)
        avg_times = nanmean(aligned_secs, 2);
        save_path = fullfile(fig.UserData.output_dir, ...
            [fig.UserData.subject_id '_' fig.UserData.run_id '_' fig.UserData.ied_id '_final_aligned.txt']);

        fid = fopen(save_path, 'w');
        fprintf(fid, '# Channel used for IED alignment: %s\n', fig.UserData.marking_channel_used);
        fprintf(fid, '%.2f\n', avg_times);
        fclose(fid);

        fprintf('[💾] Saved aligned IED times to: %s\n', save_path);

        % === Plot final average automatically
        uialert(fig, 'Main channel IED alignment complete. File saved.', 'Done');

        catch ME
        waitbar(1, wait, 'Error during alignment!');
        pause(1);
        close(wait);
        rethrow(ME);
    end

    waitbar(1, wait, 'Alignment complete.');
    pause(0.5);
    %     close(wait);

end

%% zoom in and out of the EEG traces
function changeTimeWindow(fig)
    dd = findobj(fig, 'Tag', 'TimeWindowMenu');
    val_str = dd.Value;  % e.g., '±0.5s'
    tokens = regexp(val_str, '±([\d.]+)s', 'tokens');

    if ~isempty(tokens)
        win_sec = str2double(tokens{1}{1});
        fig.UserData.win_sec = win_sec;
        fprintf('[UI] Time window updated to ±%.2f sec\n', win_sec);

        % === Automatically refresh relevant plots ===

        % Plot if current template is avgAllChannels (multi-channel plot)
        if isfield(fig.UserData, 'currentTemplateKey') && strcmp(fig.UserData.currentTemplateKey, 'avgAllChannels')
            show_avg_all_channels(fig);  % will recompute & replot with new window
        end

        % Plot the single-channel final average template (if that's active)
        if isfield(fig.UserData, 'currentTemplateKey') && strcmp(fig.UserData.currentTemplateKey, 'finalAvgTemplate')
            plot_avg_template(fig, 'finalAvgTemplate');
        end

        % Plot main channel average of all IEDs (if it was computed and shown last)
        if isfield(fig.UserData, 'currentTemplateKey') && strcmp(fig.UserData.currentTemplateKey, 'mainAvgAllIEDs')
            plot_Avg_allIEDs_OfMainCh(fig);
        end

    else
        warning('Could not parse time window string: %s', val_str);
    end
end


%%
function save_final_avg_template(fig)
    fs = fig.UserData.fs;
    template = fig.UserData.finalAvgTemplate;
    if isempty(template)
        warning('No finalAvgTemplate found to save.');
        return;
    end

    [~, nSamples] = size(template);
    win_sec = fig.UserData.win_sec;
    t = linspace(-win_sec, win_sec, nSamples);

    % Build filenames and paths
    [ied_folder, ~, ~] = fileparts(fig.UserData.ied_txt_path);
    [parent_dir, ~, ~] = fileparts(ied_folder);
    save_dir = fullfile(parent_dir, 'avg_onset_peak_times');
    if ~exist(save_dir, 'dir'), mkdir(save_dir); end

    base = sprintf('%s_%s_%s', ...
        fig.UserData.subject_id, fig.UserData.run_id, fig.UserData.ied_id);

    % ----------- PLOT -------------
    fig_avg = figure('Color', 'w');
    hold on;

    % Plot all channels in light gray
    h1 = plot(t, template', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.5); % All channels in gray

    % Highlight the channel used for IED timing adjustment
    main_label = fig.UserData.marking_channel_used;
    all_labels = fig.UserData.channelnames_bipolar;
    main_idx = find(strcmp(all_labels, main_label));
    if ~isempty(main_idx)
        hold on;
        h2 = plot(t, template(main_idx, :), 'r', 'LineWidth', 2);  % Highlighted in red
        legend('Other channels', sprintf('Main: %s', main_label), 'Location', 'northeast');
        legend([h1(1), h2], {'Other channels', sprintf('Main: %s', main_label)}, 'Location', 'northeast');
    end

    % Final touches
    title('Final Avg IED Template (All Channels)');
    xlabel('Time (s)'); ylabel('Amplitude');
    grid on; xlim([-win_sec win_sec]);
    set(gca, 'FontSize', 12);
    hold off;

    % Save
    png_path = fullfile(save_dir, [base '_final_avg_template.png']);
    saveas(fig_avg, png_path);
%     close(fig_avg);
    fprintf('[🖼] PNG export saved: %s\n', png_path);
end


%% function 9
function enable_onset_marking(fig)
    fig.UserData.isMarkingOnset = true;
    fig.UserData.isMarkingPeak = false;
    fig.Pointer = 'crosshair';
    uialert(fig, 'Click on a waveform to mark the ONSET', 'Mark Onset');   
end

%% function 10
function enable_peak_marking(fig)
    fig.UserData.isMarkingPeak = true;
    fig.UserData.isMarkingOnset = false;
    fig.Pointer = 'crosshair';
    uialert(fig, 'Click on a waveform to mark the PEAK.', 'Mark Peak');
end

%% function 11
function handle_click_on_avg(fig, src, ~)
    if ~fig.UserData.isMarkingOnset && ~fig.UserData.isMarkingPeak
        disp('[⚠️] Not in marking mode. Exiting.');
        return;
    end

    click_ax = src;  % use the actual clicked axis!

    if ~isa(click_ax, 'matlab.ui.control.UIAxes')
        return;
    end

    cp = click_ax.CurrentPoint;
    x_click = cp(1, 1);  % X in seconds

    % Determine which channel was clicked
    ax_list = fig.UserData.axesHandles;
    ch_offset = (fig.UserData.currentPage - 1) * fig.UserData.tracesPerPage;
    ax_idx = find(ax_list == click_ax, 1);

    if isempty(ax_idx)
        return;
    end
    ch_idx = ch_offset + ax_idx;

    % Save the time
    if fig.UserData.isMarkingOnset
        alreadyMarked = ~isnan(fig.UserData.onsetTimes(ch_idx));
        if alreadyMarked
            choice = uiconfirm(fig, ...
                sprintf('Onset already marked on channel %d.\nReplace with new time?', ch_idx), ...
                'Replace Onset?', ...
                'Options', {'Yes', 'No'}, ...
                'DefaultOption', 2);
            if strcmp(choice, 'No')
                return;
            end
        end
        fig.UserData.onsetTimes(ch_idx) = x_click;
        clr = [0 0.6 0];
        label_type = 'onset';

    elseif fig.UserData.isMarkingPeak
        alreadyMarked = ~isnan(fig.UserData.peakTimes(ch_idx));
        if alreadyMarked
            choice = uiconfirm(fig, ...
                sprintf('Peak already marked on channel %d.\nReplace with new time?', ch_idx), ...
                'Replace Peak?', ...
                'Options', {'Yes', 'No'}, ...
                'DefaultOption', 2);
            if strcmp(choice, 'No')
                return;
            end
        end
        fig.UserData.peakTimes(ch_idx) = x_click;
        clr = [0.5 0 0.8];
        label_type = 'peak';

    else
        label_type = 'unknown';
        return;
    end

    % Draw a vertical line at that point
    hold(click_ax, 'on');
    % Remove existing xline of the same type (if any)
    existing_lines = findall(click_ax, 'Type', 'ConstantLine', 'Color', clr);
    delete(existing_lines);

    xline(click_ax, x_click, '--', 'Color', clr, 'LineWidth', 2);
    hold(click_ax, 'off');

    % Keep marking mode ON (continuous marking mode)
    fig.Pointer = 'crosshair';  % keep crosshair
    % Do NOT reset isMarkingOnset / isMarkingPeak


    fprintf('Marked %s at %.3f sec on channel %d\n', label_type, x_click, ch_idx);
end

%% function 12
function save_onset_peak_times(fig)
    onset_offsets = fig.UserData.onsetTimes;  % relative times like -0.035s
    peak_offsets  = fig.UserData.peakTimes;
    aligned_times = fig.UserData.templateAlignedTimes_sec;  % NEW: absolute aligned times in sec
    ch_names = fig.UserData.channelnames_bipolar;
    avg_data = fig.UserData.finalAvgTemplate;

    if all(isnan(onset_offsets)) && all(isnan(peak_offsets))
        uialert(fig, 'No onset or peak times marked.', 'Nothing to Save');
        return;
    end

    subj = fig.UserData.subject_id;
    run  = fig.UserData.run_id;
    ied  = fig.UserData.ied_id;

    [ied_folder, ~, ~] = fileparts(fig.UserData.ied_txt_path);
    [parent_dir, ~, ~] = fileparts(ied_folder);
    save_dir = fullfile(parent_dir, 'avg_onset_peak_times');
    if ~exist(save_dir, 'dir'), mkdir(save_dir); end

    % Compute true times (per-IED + offset), then average across IEDs
    avg_onset_time = nan(size(ch_names));
    avg_peak_time  = nan(size(ch_names));
    peak_amp = nan(size(ch_names));

    win_sec = fig.UserData.win_sec;
    t = linspace(-win_sec, win_sec, size(avg_data, 2));

    for i = 1:numel(ch_names)
        if isnan(onset_offsets(i)) && isnan(peak_offsets(i))
            continue;
        end

        % Add offset to each aligned IED
        if ~isnan(onset_offsets(i))
            true_onsets = aligned_times + onset_offsets(i);
            avg_onset_time(i) = mean(true_onsets, 'omitnan');
        end
        if ~isnan(peak_offsets(i))
            true_peaks = aligned_times + peak_offsets(i);
            avg_peak_time(i) = mean(true_peaks, 'omitnan');

            % Get peak amplitude from avg waveform
            [~, idx] = min(abs(t - peak_offsets(i)));
            peak_amp(i) = avg_data(i, idx);
        end
    end

    % Save all 3 formats
    base_name = sprintf('%s_%s_%s_avg_onset_peak_times', subj, run, ied);
    mat_file = fullfile(save_dir, [base_name '.mat']);
    txt_file = fullfile(save_dir, [base_name '.txt']);
    xls_file = fullfile(save_dir, [base_name '.xlsx']);

    save(mat_file, 'avg_onset_time', 'avg_peak_time', 'peak_amp', 'ch_names');

    % TXT
    fid = fopen(txt_file, 'w');
    fprintf(fid, '# Averaged Real Onset/Peak Times\n');
    fprintf(fid, '# Channel\tAvgOnset(s)\tAvgPeak(s)\tPeakAmp\n');
    for i = 1:numel(ch_names)
        fprintf(fid, '%s\t%s\t%s\t%s\n', ch_names{i}, ...
            format_nan(avg_onset_time(i)), format_nan(avg_peak_time(i)), format_nan(peak_amp(i)));
    end
    fclose(fid);

    % XLSX
    T = table(ch_names(:), peak_amp(:), avg_onset_time(:), avg_peak_time(:), ...
        'VariableNames', {'Channel', 'PeakAmplitude', 'AvgOnsetTime', 'AvgPeakTime'});
    writetable(T, xls_file);

    msg = sprintf('Saved:\n%s\n%s\n%s', mat_file, txt_file, xls_file);
    uialert(fig, msg, 'Saved!');
end

function s = format_nan(val)
    if isnan(val)
        s = 'NaN';
    else
        s = sprintf('%.4f', val);
    end
end




%% function 13 
% Create Button Functions
function go_to_next_page(fig)
    maxPage = ceil(fig.UserData.totalChannels / fig.UserData.tracesPerPage);
    if fig.UserData.currentPage < maxPage
        fig.UserData.currentPage = fig.UserData.currentPage + 1;
        if isfield(fig.UserData, 'currentTemplateKey') && ...
           isfield(fig.UserData, fig.UserData.currentTemplateKey)
            plot_avg_template(fig, fig.UserData.currentTemplateKey);  % ✅ use current key
        else
            for ax = fig.UserData.axesHandles
                cla(ax);
                title(ax, '');
            end
        end
    end
end


%% function 14
function go_to_previous_page(fig)
    if fig.UserData.currentPage > 1
        fig.UserData.currentPage = fig.UserData.currentPage - 1;
        if isfield(fig.UserData, 'currentTemplateKey') && ...
           isfield(fig.UserData, fig.UserData.currentTemplateKey)
            plot_avg_template(fig, fig.UserData.currentTemplateKey);  % ✅ use current key
        else
            for ax = fig.UserData.axesHandles
                cla(ax);
                title(ax, '');
            end
        end
    end
end


%% function 15
function update_subject_info_label(fig)
    % Basic subject/run/IED info
    if isfield(fig.UserData, 'subject_id') && isfield(fig.UserData, 'run_id') && isfield(fig.UserData, 'ied_id')
        label = sprintf('Subject: %s | Run: %s | IED Type: %s', ...
            fig.UserData.subject_id, ...
            fig.UserData.run_id, ...
            fig.UserData.ied_id);
    else
        label = 'Subject Info: N/A';
    end

    % Add marking channel if available
    if isfield(fig.UserData, 'marking_channel_used') && ~isempty(fig.UserData.marking_channel_used)
        label = sprintf('%s | 📍 Marked in: %s', label, fig.UserData.marking_channel_used);
    end

    % Show current channel range
    if isfield(fig.UserData, 'currentPage') && isfield(fig.UserData, 'tracesPerPage') && isfield(fig.UserData, 'totalChannels')
        startIdx = (fig.UserData.currentPage - 1) * fig.UserData.tracesPerPage + 1;
        endIdx = min(startIdx + fig.UserData.tracesPerPage - 1, fig.UserData.totalChannels);
        label = sprintf('%s | Channels %d–%d of %d', ...
            label, startIdx, endIdx, fig.UserData.totalChannels);
    end

    % Assign to label in GUI
    fig.UserData.subjectInfoLabel.Text = label;
end

%% function 16
function [onset, peak] = load_onset_peak_txt(txtfile, ch_names)

    fid = fopen(txtfile, 'r');
    lines = {};
    while ~feof(fid)
        line = strtrim(fgetl(fid));
        if isempty(line) || startsWith(line, '#')
            continue;
        end
        lines{end+1} = line; %#ok<AGROW>
    end
    fclose(fid);

    n = numel(ch_names);
    onset = nan(1, n);
    peak  = nan(1, n);

    for i = 1:length(lines)
        parts = strsplit(lines{i}, '\t');
        if numel(parts) < 3, continue; end
        ch = parts{1};
        ch_idx = find(strcmp(ch_names, ch));
        if isempty(ch_idx), continue; end

        onset_val = str2double(parts{2});
        peak_val  = str2double(parts{3});
        if ~isnan(onset_val)
            onset(ch_idx) = onset_val;
        end
        if ~isnan(peak_val)
            peak(ch_idx) = peak_val;
        end
    end
end

%% function 17
% Superimpose all channel traces for one subject into one
function plot_all_superimposed_IEDs_sorted(fig, sort_by)
    % Inputs:
    %   fig      - your GUI handle with all UserData
    %   sort_by  - 'onset' or 'peak'

    if strcmp(sort_by, 'onset')
        latencies = fig.UserData.onsetTimes;
    elseif strcmp(sort_by, 'peak')
        latencies = fig.UserData.peakTimes;
    else
        uialert(fig, 'Invalid sort criterion. Must be "onset" or "peak".', ...
            'Sort Error');
        return;
    end
    % Add to the start of plot_all_superimposed_IEDs_sorted
    if all(isnan(latencies))
        uialert(fig, 'No onset or peak times found. Please mark them first.', ...
            'Missing Data');
        return;
    end

    avg_data = fig.UserData.finalAvgTemplate;  % or avgTemplate25
    fs = fig.UserData.fs;
    win_sec = fig.UserData.win_sec;
    t = linspace(-win_sec, win_sec, size(avg_data, 2));
    nChannels = size(avg_data, 1);

    % Choose which latency to sort by
    if strcmp(sort_by, 'onset')
        latencies = fig.UserData.onsetTimes;
    elseif strcmp(sort_by, 'peak')
        latencies = fig.UserData.peakTimes;
    else
        error('Invalid sort_by option. Use ''onset'' or ''peak''.');
    end

    % Handle NaNs: put them at the end
    [~, sorted_idx] = sortrows([isnan(latencies(:)), latencies(:)]);

    cmap = parula(nChannels);  % use perceptual color map

    figure('Name', ['IEDs sorted by ' sort_by], 'Color', 'w'); hold on;

    for plot_order = 1:nChannels
        ch = sorted_idx(plot_order);
        trace = avg_data(ch, :);
        if all(isnan(trace))
            continue;
        end
        label = sprintf('%s (%.2f s)', fig.UserData.channelnames_bipolar{ch}, latencies(ch));
        plot(t, trace, 'Color', cmap(plot_order, :), ...
            'DisplayName', label, 'LineWidth', 2);
    end

    title(sprintf('Superimposed IEDs sorted by %s latency | %s | %s | %s', ...
        sort_by, fig.UserData.subject_id, fig.UserData.run_id, fig.UserData.ied_id));
    xlabel('Time (s)');
    ylabel('Amplitude');
    legend('show', 'Interpreter', 'none', 'Location', 'eastoutside');
    grid on;

    % Build save path
    save_dir = fullfile(fig.UserData.output_dir, 'superimposed_plots');
    if ~exist(save_dir, 'dir')
        mkdir(save_dir);
    end

    % Save figure
    filename = sprintf('%s_%s_%s_superimposed_%s_sorted.png', ...
        fig.UserData.subject_id, fig.UserData.run_id, ...
        fig.UserData.ied_id, sort_by);
    exportgraphics(gcf, fullfile(save_dir, filename));
    fprintf('[💾] Saved superimposed figure to: %s\n', filename);

end

%% function 18
function [eeg_data, sampling_rate, channel_labels, ch_name_to_index] = load_eeg_bin_with_labels(eeg_bin_path)
%LOAD_EEG_BIN_WITH_LABELS Load EEG data and channel labels from ICE-format .bin files.
%
%   This function loads a binary EEG file recorded in ICE project format and
%   extracts associated channel labels from the accompanying 
%   'ProcessingLog_&_ChannelLabels.txt' file located in the same directory.
%
%   The binary file is assumed to have a filename containing metadata, such as
%   number of samples, number of channels, and sampling rate (e.g., 
%   'XXXXXSamples_64C_1000Hz.bin').
%
%   INPUTS:
%       eeg_bin_path : string
%           Full path to the EEG .bin file.
%
%   OUTPUTS:
%       eeg_data : matrix [n_channels x n_samples]
%           The EEG signal data loaded from the .bin file.
%
%       sampling_rate : scalar
%           Sampling rate extracted from the filename (in Hz).
%
%       channel_labels : cell array of strings
%           Labels of each EEG channel, extracted from the associated log file.
%
%       ch_name_to_index : containers.Map
%           A mapping from channel name to index in the EEG matrix (used for lookup).
%
%   EXAMPLE:
%       [eeg, fs, labels, ch_map] = load_eeg_bin_with_labels('/path/ICE030_Run3_XXXXXSamples_64C_1000Hz.bin');
%
%   NOTE:
%       - The associated log file must be in the same folder and include the phrase
%         'ProcessingLog_&_ChannelLabels.txt'.
%       - Channel labels are extracted from the section starting with:
%           'c_Labels (Column):'
%
%   Written by: Tahereh Rashnavadi, March 27, 2025


    % ===== Extract metadata from filename =====
    [~, bin_filename, ~] = fileparts(eeg_bin_path);
    tokens = regexp(bin_filename, '(\d+)Samples_(\d+)C_(\d+)Hz', 'tokens', 'once');
    if isempty(tokens)
        error('Could not extract metadata from filename: %s', bin_filename);
    end
    n_samples = str2double(tokens{1});
    n_channels = str2double(tokens{2});
    sampling_rate = str2double(tokens{3});

    % ===== Load EEG binary file =====
    fprintf('[INFO] Loading EEG binary: %s\n', eeg_bin_path);
    fid = fopen(eeg_bin_path, 'r');
    eeg_data = fread(fid, [n_channels, Inf], 'float32');
    fclose(fid);

    if size(eeg_data, 2) ~= n_samples
        warning('[WARNING] Sample count mismatch! Expected %d samples, got %d', ...
            n_samples, size(eeg_data, 2));
    end

    % ===== Locate the associated log file =====
    eeg_dir = fileparts(eeg_bin_path);

    % Get subject_run prefix (e.g., ICE030_Run3)
    match_prefix = regexp(bin_filename, '^(ICE\d+_Run\d+[a-zA-Z]?)', 'match', 'once');
    if isempty(match_prefix)
        error('Could not extract subject/run prefix from .bin filename: %s', bin_filename);
    end

    % Search for log file with this prefix
    log_file_struct = dir(fullfile(eeg_dir, [match_prefix '*ProcessingLog_&_ChannelLabels.txt']));
    if isempty(log_file_struct)
        error('No log file found with prefix "%s" in directory: %s', match_prefix, eeg_dir);
    end
    log_path = fullfile(eeg_dir, log_file_struct(1).name);

    % ===== Extract channel labels =====
    fid = fopen(log_path, 'r');
    lines = textscan(fid, '%s', 'Delimiter', '\n', 'Whitespace', '');
    fclose(fid);
    lines = lines{1};

    channel_labels = {};
    collecting = false;
    for i = 1:length(lines)
        line = strtrim(lines{i});
        if startsWith(line, 'c_Labels (Column):')
            collecting = true;
            continue;
        end
        if collecting
            if isempty(line) || contains(line, ':')
                break;
            else
                channel_labels = [channel_labels, strsplit(line)];
            end
        end
    end

    if length(channel_labels) ~= n_channels
        error('Mismatch between channel count in .bin (%d) and .txt file (%d)', ...
              n_channels, length(channel_labels));
    end

    ch_name_to_index = containers.Map(channel_labels, 1:n_channels);

end