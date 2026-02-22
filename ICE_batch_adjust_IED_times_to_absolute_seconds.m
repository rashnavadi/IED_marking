%% Batch adjust IED times for ALL subjects/runs/IED types
% Fix: segment start is in fMRI volumes; offset must be startVol * TR seconds.
% Output: one-column *_adjusted.txt in Tara/IED_times/<SUBJ>/adjusted_IED_times/

clear; clc;

% --- Paths ---
root_in  = '/Volumes/MIND/ICE/original_ICE/4_Data_and_Analysis';
root_out = '/Volumes/MIND/ICE/Tara/IED_times';

% --- fMRI TR (seconds) ---
TR = 1.5;

% --- Load segmentation volumes ---
S = load('segments.mat', 'segments');
segments = S.segments;

% --- Find all ICE subject folders ---
subj_dirs = dir(fullfile(root_in, 'ICE*'));
subj_dirs = subj_dirs([subj_dirs.isdir]);

if isempty(subj_dirs)
    error('No ICE subject folders found under: %s', root_in);
end

n_files_total = 0;
n_saved = 0;
n_skipped = 0;

for s = 1:numel(subj_dirs)
    subject_id = subj_dirs(s).name;

    % Skip if subject has no segmentation info
    if ~isfield(segments, subject_id)
        fprintf('[SKIP] %s: no segments entry.\n', subject_id);
        continue;
    end

    % Input IED folder
    ied_dir = fullfile(root_in, subject_id, '3_EEG', '3_Events', 'IED');
    if ~isfolder(ied_dir)
        fprintf('[SKIP] %s: no IED folder at %s\n', subject_id, ied_dir);
        continue;
    end

    % Find all IED timing txt files (original, not already adjusted)
    files = dir(fullfile(ied_dir, '*.txt'));
    files = files(~contains({files.name}, '_adjusted')); % avoid reprocessing

    if isempty(files)
        fprintf('[SKIP] %s: no IED timing txt files.\n', subject_id);
        continue;
    end

    % Output folder per subject
    out_dir = fullfile(root_out, subject_id, 'adjusted_IED_times');
    if ~exist(out_dir, 'dir'), mkdir(out_dir); end

    for k = 1:numel(files)
        fname = files(k).name;
        fpath = fullfile(ied_dir, fname);
        n_files_total = n_files_total + 1;

        % Parse: ICE###_RunX[a/b/c]_IEDY.txt
        tok = regexp(fname, '^(ICE\d+)_((Run\d+[a-z]?))_(IED\d+)\.txt$', 'tokens', 'once');
        if isempty(tok)
            fprintf('[SKIP] %s: filename format not recognized.\n', fpath);
            n_skipped = n_skipped + 1;
            continue;
        end

        subj_from_file = tok{1};
        run_name       = tok{2};   % e.g., Run1a
        ied_type       = tok{3};   % e.g., IED1

        % Sanity: should match current folder subject
        if ~strcmp(subj_from_file, subject_id)
            fprintf('[SKIP] subject mismatch folder=%s file=%s\n', subject_id, subj_from_file);
            n_skipped = n_skipped + 1;
            continue;
        end

        % Make sure segmentation exists for this run
        if ~isfield(segments.(subject_id), run_name)
            fprintf('[SKIP] %s %s %s: no segmentation for run.\n', subject_id, run_name, ied_type);
            n_skipped = n_skipped + 1;
            continue;
        end

        seg_range = segments.(subject_id).(run_name); % [startVol endVol]
        startVol  = seg_range(1);

        % Convert volume offset to seconds (the key fix)
        offset_sec = startVol * TR;

        % Read file robustly (your files can be 3 columns; first col is time)
        A = readmatrix(fpath);
        if isempty(A)
            fprintf('[SKIP] empty file: %s\n', fpath);
            n_skipped = n_skipped + 1;
            continue;
        end

        t_seg_sec = A(:,1);
        t_seg_sec = t_seg_sec(~isnan(t_seg_sec)); % drop NIL / NaN

        if isempty(t_seg_sec)
            fprintf('[SKIP] no valid times in: %s\n', fpath);
            n_skipped = n_skipped + 1;
            continue;
        end

        % Apply offset
        t_abs_sec = t_seg_sec + offset_sec;

        % Save: ICE016_Run1a_IED1_adjusted.txt (one column)
        out_file = fullfile(out_dir, sprintf('%s_%s_%s_adjusted.txt', subject_id, run_name, ied_type));
        writematrix(t_abs_sec, out_file, 'Delimiter', ' ');

        n_saved = n_saved + 1;

        fprintf('[OK] %s | %s %s %s | startVol=%d => offset=%.3f sec | n=%d\n', ...
            out_file, subject_id, run_name, ied_type, startVol, offset_sec, numel(t_abs_sec));
    end
end

fprintf('\n==== DONE ====\n');
fprintf('Total candidate files: %d\n', n_files_total);
fprintf('Saved adjusted files : %d\n', n_saved);
fprintf('Skipped              : %d\n', n_skipped);