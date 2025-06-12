

%% Launch first GUI to manually find the peak in the main channel
% By Tahere Rashnavadi
% Last update: June 11, 2025

%%
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

    % Add layout: 5 rows, 2 columns
    glayout = uigridlayout(fig, [5, 2]);
    glayout.RowHeight = {'fit', 'fit', 'fit', 'fit', '1x'};
    glayout.ColumnWidth = {'fit', '1x'};

    % === Left panel containing vertically stacked controls ===
    leftPanel = uipanel(glayout);
    leftPanel.Layout.Row = [1 5];
    leftPanel.Layout.Column = 1;

    % Vertical layout inside left panel
    leftLayout = uigridlayout(leftPanel, [10, 1]);
    leftLayout.RowHeight = repmat({'fit'}, 1, 10);
    leftLayout.ColumnWidth = {'1x'};

    % EEG file selection button
    eegButton = uibutton(leftLayout, 'Text', 'Select EEG .bin File', ...
        'ButtonPushedFcn', @(btn,event) select_eeg_file(fig));

    % IED file selection button
    iedButton = uibutton(leftLayout, 'Text', 'Select IED Manual Timings .txt File', ...
        'ButtonPushedFcn', @(btn,event) select_ied_file(fig));
    
    % Page buttons for EEG channels: Next and Previous buttons in leftLayout
    nextButton = uibutton(leftLayout, 'Text', 'Next EEG Channel ➤', ...
        'ButtonPushedFcn', @(btn, event) go_to_next_page(fig));

    prevButton = uibutton(leftLayout, 'Text', '◀ Prev. EEG Channel', ...
        'ButtonPushedFcn', @(btn, event) go_to_previous_page(fig));

    % Compute IED Average Template
    avgButton = uibutton(leftLayout, ...
    'Text', 'Compute Avg IED Template (25 IEDs)', ...
    'ButtonPushedFcn', @(btn,event) compute_avg_template(fig));

    % Button: Align All IEDs Using Template
    alignAllButton = uibutton(leftLayout, ...
    'Text', 'Align All IEDs (Cross-corr)', ...
    'ButtonPushedFcn', @(btn,event) align_all_ieds(fig));

    % Compute Final Average IED
    finalAvgButton = uibutton(leftLayout, ...
    'Text', 'Plot Final Average (All IEDs)', ...
    'ButtonPushedFcn', @(btn,event) compute_final_avg(fig));

    % Mark Onset & Peak (on Avg IED)
    markOnsetButton = uibutton(leftLayout, ...
        'Text', 'Mark Onset on Avg', ...
        'ButtonPushedFcn', @(btn,event) enable_onset_marking(fig));

    markPeakButton = uibutton(leftLayout, ...
        'Text', 'Mark Peak on Avg', ...
        'ButtonPushedFcn', @(btn,event) enable_peak_marking(fig));

    % Save Onset/Peak
    saveOnsetPeakButton = uibutton(leftLayout, ...
    'Text', 'Save Onset & Peak', ...
    'ButtonPushedFcn', @(btn,event) save_onset_peak(fig));

    fig.UserData.avgTemplate25 = [];      % Initial avg template (channels x time)
    fig.UserData.templateAlignedTimes = [];  % Refined times via cross-corr (IEDs x channels)
    fig.UserData.finalAvgTemplate = [];   % Final avg waveform after alignment

    fig.UserData.isMarkingOnset = false;
    fig.UserData.isMarkingPeak = false;
    fig.UserData.onsetTimes = [];  % Will be initialized after EEG is loaded
    fig.UserData.peakTimes  = [];

    % === Fixed plot panel for paging ===
    plotPanel = uipanel(glayout);
    plotPanel.Layout.Row = [1 5];
    plotPanel.Layout.Column = 2;
    plotPanel.Title = 'EEG traces';
    fig.UserData.plotPanel = plotPanel;

    % Create 5 axes and store them
    axesHandles = gobjects(1, 5);
    for i = 1:5
        ax = uiaxes(plotPanel);
        ax.Position = [50, 600 - i*110, 800, 100];  % Manually positioned
        ax.XLabel.String = 'Time (s)';
        ax.YLabel.String = 'Amplitude';
        axesHandles(i) = ax;
    end
    fig.UserData.axesHandles = axesHandles;

    % Store paging control info in UserData
    fig.UserData.currentPage = 1;
    fig.UserData.tracesPerPage = 5;
    fig.UserData.totalChannels = 0;  % Will be set after data is loaded
    fig.UserData.allEEGData = [];    % All EEG traces (nTime x nIEDs x nChannels)
    fig.UserData.channelNames = {};
    fig.UserData.t_ied = [];         % Time vector
    fig.UserData.fs = [];            % Sampling rate

    fig.UserData.subjectInfoLabel = uilabel(fig, ...
        'Text', 'Subject Info: N/A', ...
        'FontWeight', 'bold', ...
        'FontSize', 14, ...
        'HorizontalAlignment', 'left', ...
        'Position', [300 figHeight - 65 figWidth - 320 30]);  % <-- Moved right, trimmed left margin

    fig.UserData.adjustedIEDTimes = [];  % From manual .txt file

    main_map = load('main_channels_map.mat', 'main_channels_map');
    fig.UserData.main_channels_map = main_map.main_channels_map;
    fig.UserData.main_channels = {};  % Default fallback

end

%%
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


    fprintf('[ℹ] Loading channel info from: %s\n', full_log_path);

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

%%
function select_ied_file(fig)

    disp('Please select the manually adjusted IED timing file (.txt)');

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

    % Step 3: Read times, skip non-numeric lines ("skipped by user")
    fid = fopen(ied_txt_path, 'r');
    lines = {};
    marking_channel = 'Unknown';

    line_num = 0;
    while ~feof(fid)
        line = strtrim(fgetl(fid));
        line_num = line_num + 1;

        if line_num == 1 && contains(line, 'Channel used for IED timing adjustment')
            % Extract channel from line like: "# Channel used for IED timing adjustment: dLpIN1-dLpIN2"
            tokens = regexp(line, 'adjustment:\s*(\S+)', 'tokens');
            if ~isempty(tokens)
                marking_channel = tokens{1}{1};
            end
            continue;  % Skip first line from IED times
        end

        if ~startsWith(line, '#') && ~contains(lower(line), 'skipped') && ~isempty(line)
            lines{end+1} = line; %#ok<AGROW>
        end
    end
    fclose(fid);

    % Store and display marking channel
    fig.UserData.marking_channel_used = marking_channel;
    fig.UserData.markingChannelLabel.Text = sprintf('📍 Marked in: %s', marking_channel);

    % Convert to numeric
    manual_ied_times = [];
    for i = 1:numel(lines)
        parts = str2num(lines{i});  %#ok<ST2NM>
        if numel(parts) >= 2
            manual_ied_times(end+1,1) = parts(2);  % Keep only the second column (IED time)
        end
    end

    if isempty(manual_ied_times)
        uialert(fig, '❌ No valid IED timings found in the file.', 'Empty File');
        return;
    end

    % Step 4: Save to UserData
    fig.UserData.adjustedIEDTimes = manual_ied_times(:);  % column vector
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

    save_dir = fullfile(parent_dir, 'avg_onset_peak_times');
    basename = sprintf('%s_%s_%s', fig.UserData.subject_id, fig.UserData.run_id, fig.UserData.ied_id);
    mat_path = fullfile(save_dir, [basename '_onset_peak_times.mat']);

    if exist(mat_path, 'file')
        loaded = load(mat_path);
        if isfield(loaded, 'onset') && isfield(loaded, 'peak')
            fig.UserData.onsetTimes = loaded.onset;
            fig.UserData.peakTimes  = loaded.peak;
            fprintf('[📂] Loaded onset/peak from MAT: %s\n', mat_path);
        else
            warning('MAT file found but onset/peak fields missing: %s', mat_path);
            n = numel(fig.UserData.channelnames_bipolar);
            fig.UserData.onsetTimes = nan(1, n);
            fig.UserData.peakTimes  = nan(1, n);
        end
    else
        fprintf('[ℹ️] No MAT onset/peak file found. Starting fresh.\n');
        n = numel(fig.UserData.channelnames_bipolar);
        fig.UserData.onsetTimes = nan(1, n);
        fig.UserData.peakTimes  = nan(1, n);
    end
    
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
    if length(fig.UserData.adjustedIEDTimes) >= 1
        compute_avg_template(fig);
        % plot_avg_template is already called inside compute_avg_template
    end

    % === Step 10: Attempt to auto-load saved onset/peak if it exists ===
    % Build expected filename
    parent_dir = fileparts(ied_txt_path);  % This is IED_Cleaned
    save_dir = fullfile(parent_dir, '..', 'avg_onset_peak_times');
    expected_txt = fullfile(save_dir, sprintf('%s_%s_%s_avg_onset_peak_times.txt', ...
        txt_subject, txt_run, txt_ied));

    if exist(expected_txt, 'file')
        fprintf('[📂] Found existing onset/peak file: %s\n', expected_txt);
        [onset, peak] = load_onset_peak_txt(expected_txt, fig.UserData.channelnames_bipolar);
        fig.UserData.onsetTimes = onset;
        fig.UserData.peakTimes  = peak;
        fprintf('[✅] Onset/peak times loaded.\n');
    else
        fprintf('[ℹ️] No saved onset/peak file found. Starting fresh.\n');
    end

end
%%
function compute_avg_template(fig)
    % Parameters
    fs = fig.UserData.fs;
    eeg = fig.UserData.eeg_data;
    adjusted_times_sec = fig.UserData.adjustedIEDTimes;
    channelnames = fig.UserData.channelnames_bipolar;
    mono_map = containers.Map(fig.UserData.channelnames_mono, 1:numel(fig.UserData.channelnames_mono));

    pre_samples = round(0.5 * fs);
    post_samples = round(0.5 * fs);
    win_len = pre_samples + post_samples + 1;
    n_channels = numel(channelnames);
    n_ieds = numel(adjusted_times_sec);
    
    % Initialize output
    avg_template = zeros(n_channels, win_len);
    valid_count = zeros(n_channels, 1);  % Track how many IEDs contribute to each channel

    % Convert adjusted times to sample indices
    sample_times = round(adjusted_times_sec * fs);

    for ch = 1:n_channels
        label = channelnames{ch};
        parts = split(label, '-');
        ch1 = parts{1}; ch2 = parts{2};

        % Skip if channels not found
        if ~isKey(mono_map, ch1) || ~isKey(mono_map, ch2)
            continue;
        end
        idx1 = mono_map(ch1);
        idx2 = mono_map(ch2);

        segments = zeros(n_ieds, win_len);
        count = 0;

        for i = 1:n_ieds
            center = sample_times(i);

            start_idx = center - pre_samples;
            end_idx   = center + post_samples;

            if start_idx < 1 || end_idx > size(eeg, 2)
                continue;  % Skip if out of bounds
            end

            segment = eeg(idx1, start_idx:end_idx) - eeg(idx2, start_idx:end_idx);
            count = count + 1;
            segments(count, :) = segment;
        end

        if count > 0
            avg_template(ch, :) = mean(segments(1:count, :), 1);
            valid_count(ch) = count;
        end
    end

    fig.UserData.avgTemplate25 = avg_template;
    fprintf('[✔] Avg template computed using %d manually marked IEDs per channel.\n', n_ieds);

    % Optional: show how many IEDs contributed per channel
%     disp('Valid IED counts per channel (non-zero only):');
%     disp(find(valid_count > 0));

    % Plot
    plot_avg_template(fig, 'avgTemplate25');
end

%%
function plot_avg_template(fig, template_key)
    % template_key: 'avgTemplate25' or 'finalAvgTemplate'
    
    if ~isfield(fig.UserData, template_key) || isempty(fig.UserData.(template_key))
        uialert(fig, 'Template not found. Compute it first.', 'Missing Template');
        return;
    end

    avg_data = fig.UserData.(template_key);  % [nChannels x nTimePoints]
    fs = fig.UserData.fs;
    nChannels = size(avg_data, 1);
    nTime = size(avg_data, 2);
    t = linspace(-0.5, 0.5, nTime);  % assuming ±0.5s window

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

    % Plot traces
    for i = 1:length(axesHandles)
        ax = axesHandles(i);
        if i <= length(visibleIdx)
            chIdx = visibleIdx(i);
            trace = avg_data(chIdx, :);
            ch_name = fig.UserData.channelnames_bipolar{chIdx};

            % Handle missing or malformed channel names
            if contains(ch_name, '-')
                parts = split(ch_name, '-');
            else
                parts = {ch_name, ''};
            end
            ch1 = parts{1}; ch2 = parts{2};
            main_chs = fig.UserData.main_channels;

            % Styling
            if any(strcmp(ch1, main_chs)) || any(strcmp(ch2, main_chs))
                trace_color = [1 0 1];          % magenta
                ax.Color = [0.9 1 0.9];         % light green
            else
                trace_color = [0 0 0.3];        % dark blue
                ax.Color = 'white';
            end

            % Plot
            if all(isnan(trace))
                continue;
            end
            % Plot the trace
            lineHandle = plot(ax, t, trace, 'Color', trace_color, 'LineWidth', 1.5);

            % Show previously marked onset/peak if available
            if ~isnan(fig.UserData.onsetTimes(chIdx))
                xline(ax, fig.UserData.onsetTimes(chIdx), '--', 'Color', [0 0.6 0], 'LineWidth', 2);  % green onset
            end
            if ~isnan(fig.UserData.peakTimes(chIdx))
                xline(ax, fig.UserData.peakTimes(chIdx), '--', 'Color', [0.5 0 0.8], 'LineWidth', 2);  % purple peak
            end

            % Make sure the line does not block mouse clicks on the axes
            lineHandle.HitTest = 'off';  % Allow clicks to reach the axes

            % Configure axes for click capture
            ax.ButtonDownFcn = @(src, event) handle_click_on_avg(fig, src, event);
            ax.HitTest = 'on';
            ax.PickableParts = 'all';

            xlim(ax, [t(1), t(end)]);
            title(ax, sprintf('%s (ch %d of %d)', ch_name, chIdx, nChannels));
            xlabel(ax, 'Time (s)');
            ylabel(ax, 'Amplitude');
            grid(ax, 'on');
        else
            % Clear unused axes
            cla(ax);
            ax.Title.String = '';
            ax.ButtonDownFcn = [];
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

%% To make all IEDs consistently aligned in time (per channel) before averaging
% input: - avgTemplate25 (template), - Original IED times
% ouptu: templateAlignedTimes → [nIEDs × nChannels] of refined sample indices
% Aim: To make all IEDs consistently aligned in time (per channel) before averaging
% Aligns all originally marked IEDs to a consistent shape using cross-correlation with avgTemplate25

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
        adjusted_times_sec = fig.UserData.adjustedIEDTimes;
        channelnames = fig.UserData.channelnames_bipolar;
        mono_map = containers.Map(fig.UserData.channelnames_mono, 1:numel(fig.UserData.channelnames_mono));
    
        % Parameters
        n_ieds = numel(adjusted_times_sec);
        n_channels = numel(channelnames);
        template_len = size(avg_template, 2);
    
        pre_samples = round(0.5 * fs);
        post_samples = round(0.5 * fs);
        search_margin = round(0.2 * fs);  % search ±200 msec around center
    
        % Convert times to samples
        center_samples = round(adjusted_times_sec * fs);  % [nIEDs x 1]
    
        % Initialize output
        aligned_times = nan(n_ieds, n_channels);  % each column: one channel's aligned IEDs

        % === Main alignment loop
        for ch = 1:n_channels
            label = channelnames{ch};
            parts = split(label, '-');
            ch1 = parts{1}; ch2 = parts{2};
    
            if ~isKey(mono_map, ch1) || ~isKey(mono_map, ch2)
                continue;
            end
            idx1 = mono_map(ch1);
            idx2 = mono_map(ch2);
            template = avg_template(ch, :);
    
            for i = 1:n_ieds
                center = center_samples(i);
    
                % Define search window
                start_idx = center - search_margin - pre_samples;
                end_idx   = center + search_margin + post_samples;
    
                if start_idx < 1 || end_idx > size(eeg, 2)
                    continue;
                end
    
                full_segment = eeg(idx1, start_idx:end_idx) - eeg(idx2, start_idx:end_idx);
    
                % Slide template across center ± search_margin
                max_corr = -inf;
                best_lag = 0;
    
                for lag = -search_margin:search_margin
                    seg_start = search_margin + 1 + lag;
                    seg_end   = seg_start + template_len - 1;
    
                    if seg_start < 1 || seg_end > length(full_segment)
                        continue;
                    end
    
                    window = full_segment(seg_start:seg_end);
                    r = corrcoef(window, template);
                    c = r(1, 2);
    
                    if ~isnan(c) && c > max_corr
                        max_corr = c;
                        best_lag = lag;
                    end
                end
    
                % Update aligned time for this channel and IED
                aligned_times(i, ch) = center + best_lag;
            end
            % Optional: update waitbar per channel
            waitbar(ch / n_channels, wait, sprintf('Aligning: %d/%d channels...', ch, n_channels));
        end

        % === Store result
        fig.UserData.templateAlignedTimes = round(aligned_times);

        % === Plot final average automatically
%         compute_final_avg(fig);
        uialert(fig, 'Alignment using template completed. Now click "Plot Final Average (All IEDs)" to view.', 'Done');

        catch ME
        waitbar(1, wait, 'Error during alignment!');
        pause(1);
        close(wait);
        rethrow(ME);
    end

    waitbar(1, wait, 'Alignment complete.');
    pause(0.5);
    close(wait);    
end

%% To produce clean final per-channel average IED waveforms
% input: templateAlignedTimes, eeg_data
% output:  eeg_data	finalAvgTemplate → [nChannels × time] matrix
% Aim: To produce clean final per-channel average IED waveforms

function compute_final_avg(fig)
    % This assumes that fig.UserData.templateAlignedTimes is already computed
    % and contains a matrix: [nIEDs x nChannels] of refined sample indices

    eeg = fig.UserData.eeg_data;
    fs = fig.UserData.fs;
    aligned_samples = fig.UserData.templateAlignedTimes;
    channelnames = fig.UserData.channelnames_bipolar;
    mono_map = containers.Map(fig.UserData.channelnames_mono, 1:numel(fig.UserData.channelnames_mono));

    pre_samples = round(0.5 * fs);
    post_samples = round(0.5 * fs);
    win_len = pre_samples + post_samples + 1;
    n_channels = size(aligned_samples, 2);
    n_ieds = size(aligned_samples, 1);

    final_avg = zeros(n_channels, win_len);

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
        for i = 1:n_ieds
            center = aligned_samples(i, ch);
            start_idx = center - pre_samples;
            end_idx = center + post_samples;
            if start_idx >= 1 && end_idx <= size(eeg, 2)
                segment = eeg(idx1, start_idx:end_idx) - eeg(idx2, start_idx:end_idx);
                segments(end+1, :) = segment; %#ok<AGROW>
            end
        end

        if ~isempty(segments)
            final_avg(ch, :) = mean(segments, 1);
        end
    end

    fig.UserData.finalAvgTemplate = final_avg;

    fprintf('[✔] Final average template computed from aligned IEDs.\n');
    plot_avg_template(fig, 'finalAvgTemplate');  
end

%%
function enable_onset_marking(fig)
    fig.UserData.isMarkingOnset = true;
    fig.UserData.isMarkingPeak = false;
    fig.Pointer = 'crosshair';
    uialert(fig, 'Click on a waveform to mark the ONSET', 'Mark Onset');   
end

%%
function enable_peak_marking(fig)
    fig.UserData.isMarkingPeak = true;
    fig.UserData.isMarkingOnset = false;
    fig.Pointer = 'crosshair';
    uialert(fig, 'Click on a waveform to mark the PEAK.', 'Mark Peak');
end

%%
function handle_click_on_avg(fig, src, ~)
    disp('[🧠] handle_click_on_avg called.');

    if ~fig.UserData.isMarkingOnset && ~fig.UserData.isMarkingPeak
        disp('[⚠️] Not in marking mode. Exiting.');
        return;
    end

    click_ax = src;  % use the actual clicked axis!
    disp(['[🖱] Clicked object: ', class(click_ax)]);

    if ~isa(click_ax, 'matlab.ui.control.UIAxes')
        disp('[❌] Clicked object is not a UIAxes.');
        return;
    end

    cp = click_ax.CurrentPoint;
    x_click = cp(1, 1);  % X in seconds
    disp(['[📍] Clicked X position (seconds): ', num2str(x_click)]);

    % Determine which channel was clicked
    ax_list = fig.UserData.axesHandles;
    ch_offset = (fig.UserData.currentPage - 1) * fig.UserData.tracesPerPage;
    ax_idx = find(ax_list == click_ax, 1);

    if isempty(ax_idx)
        disp('[❌] Could not match clicked axis to one of the plotted axes.');
        return;
    end
    ch_idx = ch_offset + ax_idx;
    disp(['[🔍] Clicked on channel index: ', num2str(ch_idx)]);

    % Save the time
    if fig.UserData.isMarkingOnset
        fig.UserData.onsetTimes(ch_idx) = x_click;
        clr = [0 0.6 0];
        disp('[✅] Marked ONSET.');
        label_type = 'onset';
    elseif fig.UserData.isMarkingPeak
        fig.UserData.peakTimes(ch_idx) = x_click;
        clr = [0.5 0 0.8];
        disp('[✅] Marked PEAK.');
        label_type = 'peak';
    else
        label_type = 'unknown';
    end

    % Draw a vertical line at that point
    hold(click_ax, 'on');
    xline(click_ax, x_click, '--', 'Color', clr, 'LineWidth', 2);
    hold(click_ax, 'off');

    % Keep marking mode ON (continuous marking mode)
    fig.Pointer = 'crosshair';  % keep crosshair
    % Do NOT reset isMarkingOnset / isMarkingPeak


    fprintf('Marked %s at %.3f sec on channel %d\n', label_type, x_click, ch_idx);
end

%%
function save_onset_peak(fig)
    onset = fig.UserData.onsetTimes;
    peak  = fig.UserData.peakTimes;
    ch_names = fig.UserData.channelnames_bipolar;

    if all(isnan(onset)) && all(isnan(peak))
        uialert(fig, 'No onset or peak times marked.', 'Nothing to Save');
        return;
    end

    % Extract subject/run/IED info from UserData
    subj = fig.UserData.subject_id;
    run  = fig.UserData.run_id;
    ied  = fig.UserData.ied_id;

    % Get folder where manual IED file was loaded from
    [ied_folder, ~, ~] = fileparts(fig.UserData.ied_txt_path);  % e.g., /.../ICE014/manual_IED_marking
    [parent_dir, ~, ~] = fileparts(ied_folder);  % go one level up → /.../ICE014

    % Create new subfolder for saving
    save_dir = fullfile(parent_dir, 'avg_onset_peak_times');
    if ~exist(save_dir, 'dir')
        mkdir(save_dir);
    end

    % Construct filenames
    base_name = sprintf('%s_%s_%s_avg_onset_peak_times', subj, run, ied);
    mat_file = fullfile(save_dir, [base_name '.mat']);
    txt_file = fullfile(save_dir, [base_name '.txt']);

    % Save .mat
    save(mat_file, 'onset', 'peak', 'ch_names');

    % Write .txt
    fid = fopen(txt_file, 'w');
    fprintf(fid, '# Onset and Peak per Channel\n');
    fprintf(fid, '# Channel\tOnset (s)\tPeak (s)\n');
    for i = 1:numel(ch_names)
        onset_str = sprintf('%.4f', onset(i));
        peak_str  = sprintf('%.4f', peak(i));
        if isnan(onset(i)), onset_str = 'NaN'; end
        if isnan(peak(i)),  peak_str  = 'NaN'; end
        fprintf(fid, '%s\t%s\t%s\n', ch_names{i}, onset_str, peak_str);
    end
    fclose(fid);

    msg = sprintf('Saved onset & peak to:\n%s\n%s', mat_file, txt_file);
    uialert(fig, msg, 'Saved!');
end


%% Create Button Functions
function go_to_next_page(fig)
    maxPage = ceil(fig.UserData.totalChannels / fig.UserData.tracesPerPage);
    if fig.UserData.currentPage < maxPage
        fig.UserData.currentPage = fig.UserData.currentPage + 1;
        if isfield(fig.UserData, 'finalAvgTemplate') && ~isempty(fig.UserData.finalAvgTemplate)
            plot_avg_template(fig, 'finalAvgTemplate');
        elseif isfield(fig.UserData, 'avgTemplate25') && ~isempty(fig.UserData.avgTemplate25)
            plot_avg_template(fig, 'avgTemplate25');
        else
            % Optional: Clear axes or show a message
            for ax = fig.UserData.axesHandles
                cla(ax);
                title(ax, '');
            end
        end

    end
end

%%
function go_to_previous_page(fig)
    if fig.UserData.currentPage > 1
        fig.UserData.currentPage = fig.UserData.currentPage - 1;
        if isfield(fig.UserData, 'finalAvgTemplate') && ~isempty(fig.UserData.finalAvgTemplate)
            plot_avg_template(fig, 'finalAvgTemplate');
        elseif isfield(fig.UserData, 'avgTemplate25') && ~isempty(fig.UserData.avgTemplate25)
            plot_avg_template(fig, 'avgTemplate25');
        else
            % Optional: Clear axes or show a message
            for ax = fig.UserData.axesHandles
                cla(ax);
                title(ax, '');
            end
        end

    end
end

%%
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

%%
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

%%
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