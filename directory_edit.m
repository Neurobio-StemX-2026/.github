% Step 1: copy paste the path to the downloaded dataset
data_dir = 'DCM_Data';

% Step 2: RUN!
% The groups in your dataset
groups = {'ASD', 'Control'};
% ======================================================

% Initialize SPM (just to ensure paths are loaded)
spm('defaults', 'fmri');

fprintf('Starting batch splitting in: %s\n', data_dir);

% 1. Loop through Groups (ASD, Control)
for g = 1:length(groups)
    group_name = groups{g};
    group_path = fullfile(data_dir, group_name);
    
    % Get list of subject folders (ignoring . and ..)
    subjects = dir(group_path);
    subjects = subjects([subjects.isdir]); % Keep only directories
    subjects = subjects(~ismember({subjects.name}, {'.', '..'}));
    
    fprintf('\n--- Processing Group: %s (%d subjects) ---\n', group_name, length(subjects));
    
    % 2. Loop through Subjects
    for s = 1:length(subjects)
        sub_id = subjects(s).name;
        func_dir = fullfile(group_path, sub_id, 'func');
        
        % Check if func directory exists
        if ~exist(func_dir, 'dir')
            fprintf('Skipping %s (No func folder)\n', sub_id);
            continue;
        end
        
        % 3. Find the Functional File
        % We look for *.nii or *.nii.gz
        file_search = dir(fullfile(func_dir, '*func_preproc.nii*'));
        
        if isempty(file_search)
             % Try simpler name if the long one isn't found
             file_search = dir(fullfile(func_dir, 'func.nii*'));
        end
        
        if isempty(file_search)
            fprintf('WARNING: No functional file found for %s\n', sub_id);
            continue;
        end
        
        filename = file_search(1).name;
        full_file_path = fullfile(func_dir, filename);
        
        fprintf('Processing %s...', sub_id);
        
        % 4. Handle .gz files (SPM split requires .nii)
        [~, ~, ext] = fileparts(filename);
        if strcmpi(ext, '.gz')
            % Unzip it
            try
                gunzip(full_file_path);
                % Update path to the unzipped .nii file
                filename = filename(1:end-3); % remove .gz
                full_file_path = fullfile(func_dir, filename);
            catch ME
                fprintf(' FAILED to unzip: %s\n', ME.message);
                continue;
            end
        end
        
        % 5. Execute SPM File Split
        try
            % This creates func_0001.nii, func_0002.nii etc. in the same folder
            spm_file_split(full_file_path);
            fprintf(' Done! (Split into 3D volumes)\n');
            
            % Optional: Delete the large 4D .nii file to save space? 
            % (Keep the original .gz if you want backups)
            % delete(full_file_path); 
            
        catch ME
            fprintf(' SPM ERROR: %s\n', ME.message);
        end
    end
end

fprintf('\nBatch processing complete!\n');
