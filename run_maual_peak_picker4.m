



%%
function manual_ied_loader_gui()
    % === NEW: Get screen size and center the GUI ===
    screenSize = get(0, 'ScreenSize');  % [left bottom width height]
    figWidth = 1000;
    figHeight = 700;
    left = (screenSize(3) - figWidth) / 2;
    bottom = (screenSize(4) - figHeight) / 2;

    % Create the main GUI window (centered and larger)
    fig = uifigure('Name', 'Manual IED Loader GUI', ...
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

    % EEG channel selection button
    channelButton = uibutton(leftLayout, 'Text', 'Select EEG Channels', ...
        'ButtonPushedFcn', @(btn,event) select_channels(fig));

    % IED file selection button
    iedButton = uibutton(leftLayout, 'Text', 'Select IED .txt File', ...
        'ButtonPushedFcn', @(btn,event) select_ied_file(fig));

    % Page buttons for EEG channels: Next and Previous buttons in leftLayout
    nextButton = uibutton(leftLayout, 'Text', 'Next EEG Channel ➤', ...
        'ButtonPushedFcn', @(btn, event) go_to_next_page(fig));

    prevButton = uibutton(leftLayout, 'Text', '◀ Prev. EEG Channel', ...
        'ButtonPushedFcn', @(btn, event) go_to_previous_page(fig));

    % Page buttons for IEDs
    nextIEDButton = uibutton(leftLayout, 'Text', 'Next IED ▶', ...
        'ButtonPushedFcn', @(btn, event) go_to_next_ied(fig));

    prevIEDButton = uibutton(leftLayout, 'Text', '◀ Prev. IED', ...
        'ButtonPushedFcn', @(btn, event) go_to_previous_ied(fig));

    % Confirm/save peak button (same width as others now)
    confirmButton = uibutton(leftLayout, ...
        'Text', 'Save Peak Time', ...
        'ButtonPushedFcn', @(btn, event) save_peak_time(fig));
    fig.UserData.confirmButton = confirmButton;

    saveAllButton = uibutton(leftLayout, ...
        'Text', 'Save All Peaks', ...
        'ButtonPushedFcn', @(btn, event) save_all_confirmed_peaks(fig));

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
    % subject information on the GUI
    fig.UserData.subjectInfoLabel = uilabel(fig, ...
    'Text', 'Subject Info: N/A', ...
    'FontWeight', 'bold', ...
    'FontSize', 14, ...
    'Position', [350 figHeight - 80 400 25]);  % Adjust position as needed

    main_map = load('main_channels_map.mat', 'main_channels_map');
    fig.UserData.main_channels_map = main_map.main_channels_map;
    fig.UserData.main_channels = {};  % Default fallback

end

%%
function select_eeg_file(fig)

    % File picker
    [bin_file, bin_path] = uigetfile('*.bin', 'Select the EEG binary file (.bin)');

    if isequal(bin_file, 0)
        disp('❌ EEG file selection cancelled.')
        return;
    end

    % Store full path in figure's UserData
    fig.UserData.eeg_bin_path = fullfile(bin_path, bin_file);
    fprintf('[✔] EEG file selected: %s\n', fig.UserData.eeg_bin_path);

    % === FIX: Make sure this line exists ===
    [eeg_data, fs, channel_labels, ch_map] = load_eeg_bin_with_labels(fig.UserData.eeg_bin_path);
    
    fig.UserData.eeg_data = eeg_data;
    fig.UserData.fs = fs;

    disp('EEG file loaded successfully.');
end

%%
function select_channels(fig)
    % Step 1: Prompt user to load processing log text file
    [log_file, log_path] = uigetfile('*.txt', '📄 Select processing log (.txt) file');
    if isequal(log_file, 0)
        uialert(fig, '❌ Log file selection cancelled.', 'Cancelled');
        return;
    end
    full_log_path = fullfile(log_path, log_file);

    % Step 2: Read file line-by-line
    fid = fopen(full_log_path, 'r');
    lines = {};
    while ~feof(fid)
        lines{end+1} = fgetl(fid); %#ok<AGROW>
    end
    fclose(fid);

    % Step 3: Extract sampling rate, number of channels, and labels
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

    fprintf('Loaded sampling rate: %d Hz\n', s_Rate);
    fprintf('Found %d monopolar channels, %d bipolar pairs.\n', ...
        numel(c_Labels), numel(bipolar_labels));

    % Placeholder: confirmation for now
    disp(['Extracted ', num2str(numel(bipolar_labels)), ' bipolar traces. Sampling Rate: ', num2str(s_Rate), ' Hz']);
end

%%
function select_ied_file(fig)

    disp('Please select the adjusted IED timing file (.txt)');

    % Step 1: Let user select .txt file
    [ied_file, ied_path] = uigetfile('*.txt', 'Select adjusted IED timing file');
    if isequal(ied_file, 0)
        disp('[⚠] IED file selection cancelled.');
        return;
    end
    ied_txt_path = fullfile(ied_path, ied_file);
    fig.UserData.ied_txt_path = ied_txt_path;

    % Step 2: Parse subject/run/IED from filename
    [~, ied_filename_noext, ~] = fileparts(ied_file);
    tokens = regexp(ied_filename_noext, '(ICE\d+)_([a-zA-Z0-9]+)_([a-zA-Z0-9]+)_adjusted', 'tokens');

    disp(['Full IED filename: ' ied_file]);

    if isempty(tokens)
        warning('❌ Could not parse subject/run/IED info from filename: %s', ied_file);
        return;
    end
    tokens = tokens{1};  % extract matched tokens
    txt_subject = tokens{1};
    txt_run     = tokens{2};
    txt_ied     = tokens{3};

    % Step 3: Compare with .bin EEG file (must already be loaded!)
    if isfield(fig.UserData, 'eeg_bin_path')
        [~, bin_name, ~] = fileparts(fig.UserData.eeg_bin_path);
        
        % Parse only subject and run from .bin file
        bin_tokens = regexp(bin_name, '(ICE\d+)_([a-zA-Z0-9]+)_IED', 'tokens');  % no 3rd token
        if isempty(bin_tokens)
            warning('❌ Could not parse subject/run info from EEG filename: %s', fig.UserData.eeg_bin_path);
            return;
        end
        bin_tokens  = bin_tokens{1};
        bin_subject = bin_tokens{1};
        bin_run     = bin_tokens{2};

        % === Match subject and base run ID only ===
        % Strip trailing letter in run_id (e.g., "Run1a" → "Run1")
        run_base_txt = regexp(txt_run, '^Run\d+', 'match', 'once');
        run_base_bin = bin_run;  % Already base form in .bin

        if ~strcmp(txt_subject, bin_subject) || ~strcmp(run_base_txt, run_base_bin)
            warning(['❌ Mismatch between EEG .bin and IED .txt files:\n' ...
                '    EEG: %s_%s\n    IED: %s_%s'], ...
                bin_subject, bin_run, txt_subject, txt_run);
            return;
        end
    end

    % Step 4: Load IED times
    ied_times_sec = load(ied_txt_path);
    fig.UserData.original_ied_times_sec = ied_times_sec(:);  % ensure column
    fig.UserData.original_ied_samples = round(ied_times_sec * fig.UserData.fs);
    fig.UserData.confirmed_peak_times = nan(length(fig.UserData.original_ied_samples), 1);  % Initialize empty vector
    fig.UserData.currentIED = 1;
    % Reset to page 1 and trigger the first update
    fig.UserData.currentPage = 1;
    fig.UserData.totalIEDs = length(fig.UserData.original_ied_samples);
    % Step 5: Save identifiers
    fig.UserData.subject_id = txt_subject;
    fig.UserData.run_id     = txt_run;
    fig.UserData.ied_id     = txt_ied;

    update_plot(fig);  % Already updates the label inside

    % Compose key and load main channels
    key = sprintf('%s_%s', txt_subject, txt_ied);
    if isKey(fig.UserData.main_channels_map, key)
        fig.UserData.main_channels = fig.UserData.main_channels_map(key);
        fprintf('Main channels for %s loaded.\n', key);
    else
        fig.UserData.main_channels = {};  % fallback: no highlight
        warning('No main channels found for %s', key);
    end

    % Sanity check before plotting
    if isfield(fig.UserData, 'eeg_data') && isfield(fig.UserData, 'channelnames_bipolar')
        fig.UserData.totalChannels = length(fig.UserData.channelnames_bipolar);
        update_plot(fig);
    else
        uialert(fig, 'Missing EEG data or channel labels. Load EEG and channel info first.', 'Missing Data');
    end

    % go through different IEDs
    fig.UserData.currentIED = 1;
    fig.UserData.totalIEDs = length(fig.UserData.original_ied_samples);

    % Success message
    fprintf('IED file selected: %s\n', ied_txt_path);
    fprintf('Parsed: subject_id = %s, run_id = %s, ied_id = %s\n', ...
        txt_subject, txt_run, txt_ied);
end

%%
function update_plot(fig)
    eeg = fig.UserData.eeg_data;
    fs = fig.UserData.fs;
    ied_samples = fig.UserData.original_ied_samples;
    t = (-0.5 : 1/fs : 0.5);  % 1 second window
    win_samples = (length(t) - 1) / 2;
    
    % Use the currently selected IED
    ied_idx = fig.UserData.currentIED;
    center = ied_samples(ied_idx);

    if center - win_samples < 1 || center + win_samples > size(eeg, 2)
        disp('[⚠] First IED is out of bounds.');
        return;
    end

    % Paging parameters: Compute vertical size needed
    n_bipolar = length(fig.UserData.channelnames_bipolar);
    perPage = fig.UserData.tracesPerPage;
    page = fig.UserData.currentPage;

    startIdx = (page - 1) * perPage + 1;
    endIdx = min(startIdx + perPage - 1, n_bipolar);
    visibleIdx = startIdx:endIdx;

    % Get fixed axes
    axesHandles = fig.UserData.axesHandles;

    % Create mapping from channel names to indices
    ch_map = containers.Map(fig.UserData.channelnames_mono, 1:numel(fig.UserData.channelnames_mono));

    % Clear all axes
    for ax = axesHandles
        cla(ax);
        ax.Title.String = '';
    end

    % Plot only channels for current page
    fig.UserData.axesList = [];
    for plotIdx = 1:length(visibleIdx)
        chIdx = visibleIdx(plotIdx);
        label = fig.UserData.channelnames_bipolar{chIdx};
        parts = split(label, '-');
        ch1 = parts{1}; ch2 = parts{2};

        if ~isKey(ch_map, ch1) || ~isKey(ch_map, ch2)
            continue;
        end

        idx1 = ch_map(ch1);
        idx2 = ch_map(ch2);
        trace = eeg(idx1, center - win_samples : center + win_samples) - ...
            eeg(idx2, center - win_samples : center + win_samples);

        ax = axesHandles(plotIdx);
        main_chs = fig.UserData.main_channels;
        if any(strcmp(ch1, main_chs)) || any(strcmp(ch2, main_chs))
            trace_color = [1 0 1];  % Magenta: [R G B]
        else
            trace_color = [0 0 0.3];  % Dark blue for others
        end

        if ~isfield(fig.UserData, 'highlightLabel') || ~isvalid(fig.UserData.highlightLabel)
            fig.UserData.highlightLabel = uilabel(fig.UserData.plotPanel, ...
                'Text', '🟣 = trace includes main EEG channel', ...
                'FontAngle', 'italic', ...
                'FontSize', 12, ...
                'Position', [20, 5, 300, 20]);
        end

        plot(ax, t, trace, 'Color', trace_color, 'LineWidth', 1.5);
        title(ax, sprintf('%s (IED %d of %d)', label, ied_idx, fig.UserData.totalIEDs));
        xlabel(ax, 'Time (s)');
        ylabel(ax, 'Amplitude');
        ylim(ax, 'auto');
        grid(ax, 'on');

        fig.UserData.axesList = [fig.UserData.axesList, ax];
    end

    % Update top subject info label
    fig.UserData.subjectInfoLabel.Text = ...
        sprintf('Subject: %s | Run: %s | IED Type: %s (%d of %d)', ...
        fig.UserData.subject_id, ...
        fig.UserData.run_id, ...
        fig.UserData.ied_id, ...
        fig.UserData.currentIED, ...
        fig.UserData.totalIEDs);

    % Add the Red Vertical Line: Shared red line for all plots
    ied_idx = fig.UserData.currentIED;
    peak_time = fig.UserData.confirmed_peak_times(ied_idx);

    if isnan(peak_time)
        x_init = 0;  % default if user hasn't picked
    else
        x_init = peak_time - fig.UserData.original_ied_times_sec(ied_idx);  % relative to IED center
    end
    line_handle = gobjects(numel(fig.UserData.axesList), 1);
    for i = 1:numel(fig.UserData.axesList)
        ax = fig.UserData.axesList(i);
        hold(ax, 'on');
        line_handle(i) = xline(ax, x_init, 'r-', 'LineWidth', 2, 'DisplayName', 'IED Peak');
    end

    fig.UserData.peakLineHandles = line_handle;
    fig.UserData.peakLineTime = x_init;

    % Enable line dragging
    set(fig, 'WindowButtonDownFcn', @(~,~) start_drag(fig));
    set(fig, 'WindowButtonUpFcn', @(~,~) stop_drag(fig));
    set(fig, 'WindowButtonMotionFcn', @(~,~) drag_line(fig));

end

%% Add Dragging Functions
function start_drag(fig)
    fig.UserData.isDragging = true;
end

function stop_drag(fig)
    fig.UserData.isDragging = false;
    fig.UserData.confirmed_peak_times(fig.UserData.currentIED) = ...
    fig.UserData.original_ied_times_sec(fig.UserData.currentIED) + fig.UserData.peakLineTime;
end

function drag_line(fig)
    if isfield(fig.UserData, 'isDragging') && fig.UserData.isDragging
        cp = get(fig, 'CurrentPoint');
        fig_pos = fig.Position;
        x_mouse = cp(1) / fig_pos(3);  % normalized
        x_val = -0.5 + x_mouse;  % assumes xlim = [-0.5 0.5]
        for lh = fig.UserData.peakLineHandles(:)'
            lh.Value = x_val;
        end
        fig.UserData.peakLineTime = x_val;
    end
end


%% Saving the Time to .txt
function save_peak_time(fig)
    if ~isfield(fig.UserData, 'currentIED') || ~isfield(fig.UserData, 'peakLineTime')
        uialert(fig, 'No IED or peak selected.', 'Error');
        return;
    end

    ied_idx = fig.UserData.currentIED;
    rel_time = fig.UserData.peakLineTime;

    % Compute absolute corrected time
    abs_time = fig.UserData.original_ied_times_sec(ied_idx) + rel_time;

    % Store the marked time (relative to center)
    fig.UserData.confirmed_peak_times(ied_idx) = abs_time;

    fig.UserData.confirmed_peak_times(fig.UserData.currentIED) = ...
    fig.UserData.original_ied_times_sec(fig.UserData.currentIED) + fig.UserData.peakLineTime;

    % ✅ Show absolute time in the confirmation message
    uialert(fig, ...
        sprintf('Peak saved for IED %d at %.3f s (absolute)', ied_idx, abs_time), ...
        'Success', 'Icon', 'success');
end


%%
function save_all_confirmed_peaks(fig)
    if ~isfield(fig.UserData, 'confirmed_peak_times') || ...
       ~isfield(fig.UserData, 'original_ied_times_sec')
        uialert(fig, 'Missing confirmed peaks or original times.', 'Error');
        return;
    end

    confirmed = fig.UserData.confirmed_peak_times;
    original = fig.UserData.original_ied_times_sec;

    if length(confirmed) ~= length(original)
        uialert(fig, 'Mismatch between peak and IED time vectors.', 'Error');
        return;
    end

    % Compute absolute peak times
    adjusted_peak_times = original + confirmed;

    % Warn if any NaNs still exist (unconfirmed IEDs)
    unconfirmed = isnan(adjusted_peak_times);
    if any(unconfirmed)
        choice = uiconfirm(fig, ...
            sprintf('There are %d unmarked IEDs. Save anyway?', sum(unconfirmed)), ...
            'Unconfirmed IEDs', ...
            'Options', {'Yes', 'Cancel'}, ...
            'DefaultOption', 1, 'CancelOption', 2);
        if strcmp(choice, 'Cancel')
            return;
        end
    end

    % Get subject/run/IED info from UserData
    subj = fig.UserData.subject_id;
    run  = fig.UserData.run_id;
    ied  = fig.UserData.ied_id;
    orig_path = fig.UserData.ied_txt_path;

    % Step up to subject-level folder (e.g., /Volumes/.../IED_times/ICE023/)
    subject_dir = fileparts(fileparts(orig_path));  % one level up from 'adjusted_IED_times'
    save_dir = fullfile(subject_dir, 'manual_IED_marking');

    % Create folder if it doesn't exist
    if ~exist(save_dir, 'dir')
        mkdir(save_dir);
    end

    % Construct filename and full path
    save_name = sprintf('%s_%s_%s_manual.txt', subj, run, ied);
    save_full_path = fullfile(save_dir, save_name);

    % Save the adjusted times
    writematrix(adjusted_peak_times, save_full_path);

    % Notify user
    uialert(fig, sprintf('Adjusted IED peak times saved to:\n%s', save_full_path), ...
        'Saved!', 'Icon', 'success');
end


%%
function plot_current_page(fig)
    allData = fig.UserData.allEEGData;
    t_ied = fig.UserData.t_ied;
    ch_names = fig.UserData.channelNames;
    axesHandles = fig.UserData.axesHandles;

    totalCh = fig.UserData.totalChannels;
    perPage = fig.UserData.tracesPerPage;
    page = fig.UserData.currentPage;
    
    startIdx = (page - 1) * perPage + 1;
    endIdx = min(startIdx + perPage - 1, totalCh);
    
    % Clear all axes first
    for ax = axesHandles
        cla(ax);
        ax.Title.String = '';
    end

    % Plot the channels for this page
    for i = startIdx:endIdx
        ax = axesHandles(i - startIdx + 1);
        data = squeeze(allData(:, 1, i));  % first IED
        plot(ax, t_ied, data, 'k');
        ax.Title.String = sprintf('%s (IED 1)', ch_names{i});
    end
end

%% Create Button Functions
function go_to_next_page(fig)
    maxPage = ceil(fig.UserData.totalChannels / fig.UserData.tracesPerPage);
    if fig.UserData.currentPage < maxPage
        fig.UserData.currentPage = fig.UserData.currentPage + 1;
        update_plot(fig);
    end
end

%%
function go_to_previous_page(fig)
    if fig.UserData.currentPage > 1
        fig.UserData.currentPage = fig.UserData.currentPage - 1;
        update_plot(fig);
    end
end

%% go_to_next_ied.m:
function go_to_next_ied(fig)
    if fig.UserData.currentIED < fig.UserData.totalIEDs
        fig.UserData.currentIED = fig.UserData.currentIED + 1;
        update_plot(fig);
    end
end

%% go_to_previous_ied.m:
function go_to_previous_ied(fig)
    if fig.UserData.currentIED > 1
        fig.UserData.currentIED = fig.UserData.currentIED - 1;
        update_plot(fig);
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

