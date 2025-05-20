



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
    leftLayout = uigridlayout(leftPanel, [5, 1]);
    leftLayout.RowHeight = {'fit', 'fit', 'fit', 'fit', 'fit'};
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

    % Dropdown to select bipolar trace
    bipolarDropdown = uidropdown(leftLayout, ...
        'Items', {'[Load channels first]'}, ...
        'ValueChangedFcn', @(dd,event) update_plot(fig));
    fig.UserData.bipolarDropdown = bipolarDropdown;

    % Confirm/save peak button (same width as others now)
    confirmButton = uibutton(leftLayout, ...
        'Text', 'Save Peak Time', ...
        'ButtonPushedFcn', @(btn, event) save_peak_time(fig));
    fig.UserData.confirmButton = confirmButton;

    % === Right side: EEG trace panel ===
    scrollPanel = uipanel(glayout);
    scrollPanel.Layout.Row = [1 5];  % spans all rows
    scrollPanel.Layout.Column = 2;
    scrollPanel.Scrollable = 'on';
    scrollPanel.Title = 'EEG traces';
    fig.UserData.scrollPanel = scrollPanel;

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
    fig.UserData.fs = fs;  % Optional re-store in case it wasn't done earlier

    % Optional: show confirmation alert
    disp('[✔] EEG file loaded successfully.');

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

    fprintf('[✔] Loaded sampling rate: %d Hz\n', s_Rate);
    fprintf('[✔] Found %d monopolar channels, %d bipolar pairs.\n', ...
        numel(c_Labels), numel(bipolar_labels));

    % Placeholder: confirmation for now
    disp(['[✔] Extracted ', num2str(numel(bipolar_labels)), ' bipolar traces. Sampling Rate: ', num2str(s_Rate), ' Hz']);

    % Update dropdown and slider in GUI
    fig.UserData.bipolarDropdown.Items = bipolar_labels;
    fig.UserData.bipolarDropdown.Value = bipolar_labels{1};

end

%%
function select_ied_file(fig)

    disp('[INFO] Please select the adjusted IED timing file (.txt)');

    % Step 1: Let user select .txt file
    [ied_file, ied_path] = uigetfile('*.txt', '📄 Select adjusted IED timing file');
    if isequal(ied_file, 0)
        disp('[⚠] IED file selection cancelled.');
        return;
    end
    ied_txt_path = fullfile(ied_path, ied_file);
    fig.UserData.ied_txt_path = ied_txt_path;

    % Step 2: Parse subject/run/IED from filename
    [~, ied_filename_noext, ~] = fileparts(ied_file);
    tokens = regexp(ied_filename_noext, '(ICE\d+)_([a-zA-Z0-9]+)_([a-zA-Z0-9]+)_adjusted', 'tokens');

    disp(['[DEBUG] Full IED filename: ' ied_file]);
    disp(['[DEBUG] Parsed filename part: ' ied_filename_noext]);

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
        bin_tokens = bin_tokens{1};
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
    fig.UserData.original_ied_samples = round(ied_times_sec * fig.UserData.fs);

    % Step 5: Save identifiers
    fig.UserData.subject_id = txt_subject;
    fig.UserData.run_id     = txt_run;
    fig.UserData.ied_id     = txt_ied;

    % ✅ Success message
    fprintf('IED file selected: %s\n', ied_txt_path);
    fprintf('Parsed: subject_id = %s, run_id = %s, ied_id = %s\n', ...
        txt_subject, txt_run, txt_ied);
end


%%
function update_plot(fig)
    disp('[DEBUG] update_plot called');
    disp(['Num IEDs: ' num2str(length(fig.UserData.original_ied_samples))]);
    disp(['EEG size: ' mat2str(size(fig.UserData.eeg_data))]);

    eeg = fig.UserData.eeg_data;
    fs = fig.UserData.fs;
    ied_samples = fig.UserData.original_ied_samples;
    t = (-0.5 : 1/fs : 0.5);  % 1 second window
    win_samples = (length(t) - 1) / 2;

    % === Use only the first IED for this view ===
    ied_idx = 1;
    center = ied_samples(ied_idx);

    if center - win_samples < 1 || center + win_samples > size(eeg, 2)
        disp('[⚠] First IED is out of bounds.');
        return;
    end

    % Clear previous axes
    delete(allchild(fig.UserData.scrollPanel));

    % Compute vertical size needed
    n_bipolar = length(fig.UserData.channelnames_bipolar);
    panel_height = max(800, 140 * n_bipolar);
    fig.UserData.scrollPanel.Position(4) = panel_height;

    % Create mapping from channel names to indices
    ch_map = containers.Map(fig.UserData.channelnames_mono, 1:numel(fig.UserData.channelnames_mono));

    % Loop through all bipolar pairs
    for i = 1:n_bipolar
        label = fig.UserData.channelnames_bipolar{i};
        parts = split(label, '-');
        ch1 = parts{1};
        ch2 = parts{2};

        if ~isKey(ch_map, ch1) || ~isKey(ch_map, ch2), continue; end
        idx1 = ch_map(ch1);
        idx2 = ch_map(ch2);

        trace = eeg(idx1, center - win_samples : center + win_samples) - ...
                eeg(idx2, center - win_samples : center + win_samples);

        ax = uiaxes(fig.UserData.scrollPanel);
        ax.Position = [10, panel_height - i*140, 950, 120];
        plot(ax, t, trace, 'k');
        title(ax, sprintf('%s (IED %d)', label, ied_idx));
        xlabel(ax, 'Time (s)');
        ylabel(ax, 'Amplitude');
        ylim(ax, [-100 100]);  % Optional for scaling
        grid(ax, 'on');

        if i == 1
            fig.UserData.axesList = [];  % reset on every call
        end
        fig.UserData.axesList = [fig.UserData.axesList, ax];


    end

    % Add the Red Vertical Line: Shared red line for all plots
    x_init = 0;  % initial at time zero
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
    t_sec = fig.UserData.peakLineTime;
    [file,path] = uiputfile('*.txt','Save peak time as');
    if isequal(file,0)
        disp('User canceled saving.');
        return;
    end
    writematrix(t_sec, fullfile(path, file));
    disp(['Saved peak time to: ' fullfile(path, file)]);
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

