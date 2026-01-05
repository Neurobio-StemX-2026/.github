% This batch script analyses the resting state fMRI dataset
% available from the SPM website using DCM:
%   http://www.fil.ion.ucl.ac.uk/spm/data/spDCM/
% as described in the SPM manual:
%   http://www.fil.ion.ucl.ac.uk/spm/doc/spm12_manual.pdf#Chap:DCM_rsfmri
%__________________________________________________________________________
% Copyright (C) 2020 Wellcome Centre for Human Neuroimaging
% Adeel Razi & Guillaume Flandin
% ------------------------------------------------------------------------------
data_dir = 'C:\Users\LOQ\Desktop\Uni Work\IBA STEMx\ABIDE\IBA Dataset';

% The groups in your dataset
groups = {'ASD', 'Control'};

% 1. Loop through Groups (ASD, Control)
for g = 1:length(groups)
    group_name = groups{g};
    group_path = fullfile(data_dir, group_name);
    
    % Get list of subject folders
    subjects = dir(group_path);
    subjects = subjects([subjects.isdir]); 
    subjects = subjects(~ismember({subjects.name}, {'.', '..'}));
    
    fprintf('\n--- Processing Group: %s (%d subjects) ---\n', group_name, length(subjects));
    
    % 2. Loop through Subjects
    for s = 1:length(subjects)
        sub_id = subjects(s).name;
        sub_path = fullfile(group_path, sub_id);
        func_dir = fullfile(sub_path, 'func');
        
        if ~exist(func_dir, 'dir')
            fprintf('Skipping %s (No func folder)\n', sub_id);
            continue;
        end
        
        % 3. Find the Functional File
        file_search = dir(fullfile(func_dir, '*func_preproc.nii*'));
        if isempty(file_search)
             file_search = dir(fullfile(func_dir, 'func.nii*'));
        end
        
        if isempty(file_search)
            fprintf('WARNING: No functional file found for %s\n', sub_id);
            continue;
        end
        
        filename = file_search(1).name;
        
        % Initialize SPM
        spm('Defaults','fMRI');
        spm_jobman('initcfg');
        RT = 2; % repetition time for ABIDE
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % INITIAL GLM FOR EXTRACTING of ROI TIME SERIES
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        glmdir = fullfile(sub_path, 'glm'); 
        if ~exist(glmdir, 'dir'), mkdir(glmdir); end
        
        f = fullfile(func_dir, filename); 
        
        clear matlabbatch;
        matlabbatch{1}.spm.stats.fmri_spec.dir          = cellstr(glmdir);
        matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'scans';
        matlabbatch{1}.spm.stats.fmri_spec.timing.RT    = RT;
        matlabbatch{1}.spm.stats.fmri_spec.sess.scans   = cellstr(f);
        matlabbatch{1}.spm.stats.fmri_spec.mthresh      = -Inf;
        matlabbatch{2}.spm.stats.fmri_est.spmmat        = cellstr(fullfile(glmdir,'SPM.mat'));
        
        spm_jobman('run',matlabbatch);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ROI extraction
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        clear matlabbatch;
        matlabbatch{1}.spm.util.voi.spmmat  = cellstr(fullfile(glmdir,'SPM.mat'));
        matlabbatch{1}.spm.util.voi.adjust  = NaN;
        matlabbatch{1}.spm.util.voi.session = 1;
        matlabbatch{1}.spm.util.voi.name    = 'PCC';
        matlabbatch{1}.spm.util.voi.roi{1}.sphere.centre = [0 -52 26];
        matlabbatch{1}.spm.util.voi.roi{1}.sphere.radius = 8;
        matlabbatch{1}.spm.util.voi.roi{2}.mask.image    = cellstr(fullfile(glmdir,'mask.nii'));
        matlabbatch{1}.spm.util.voi.expression = 'i1 & i2';
        
        matlabbatch{2} = matlabbatch{1};
        matlabbatch{2}.spm.util.voi.name = 'mPFC';
        matlabbatch{2}.spm.util.voi.roi{1}.sphere.centre = [3 54 -2];
        
        matlabbatch{3} = matlabbatch{1};
        matlabbatch{3}.spm.util.voi.name = 'LIPC';
        matlabbatch{3}.spm.util.voi.roi{1}.sphere.centre = [-50 -63 32];
        
        matlabbatch{4} = matlabbatch{1};
        matlabbatch{4}.spm.util.voi.name = 'RIPC';
        matlabbatch{4}.spm.util.voi.roi{1}.sphere.centre = [48 -69 35];
        
        spm_jobman('run',matlabbatch);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % SPECIFY & ESTIMATE DCM
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        clear DCM;
        cd(glmdir)
        
        load('VOI_PCC_1.mat');  DCM.xY(1) = xY;
        load('VOI_mPFC_1.mat'); DCM.xY(2) = xY; 
        load('VOI_LIPC_1.mat'); DCM.xY(3) = xY; 
        load('VOI_RIPC_1.mat'); DCM.xY(4) = xY; 
        
        v = length(DCM.xY(1).u); 
        n = length(DCM.xY);      
        DCM.v = v; DCM.n = n;
        
        DCM.Y.dt  = RT;
        DCM.Y.X0  = DCM.xY(1).X0;
        DCM.Y.Q   = spm_Ce(ones(1,n)*v);
        for i = 1:DCM.n
            DCM.Y.y(:,i)  = DCM.xY(i).u;
            DCM.Y.name{i} = DCM.xY(i).name;
        end
        
        DCM.U.u    = zeros(v,1);
        DCM.U.name = {'null'};         
        DCM.a      = ones(n,n);
        DCM.b      = zeros(n,n,0);
        DCM.c      = zeros(n,0);
        DCM.d      = zeros(n,n,0);
        DCM.TE     = 0.04;
        DCM.delays = repmat(RT,DCM.n,1);
        
        DCM.options.nonlinear  = 0;
        DCM.options.two_state  = 0;
        DCM.options.stochastic = 1;
        DCM.options.analysis   = 'CSD';
        DCM.options.induced    = 1;
        
        % MODIFICATION: Set Name to Subject ID for graph titles
        DCM.name = sub_id; 
        
        str = sprintf('DCM_%s_DMN.mat', sub_id);
        save(fullfile(glmdir, str), 'DCM');
        
        % Estimate the DCM
        DCM = spm_dcm_fmri_csd(fullfile(glmdir, str));
        
        % MODIFICATION: Completion message and 2-second pause
        fprintf('\nFinished processing Subject: %s\n', sub_id);
        pause(2);
        
    end
end
fprintf('\nBatch processing complete!\n');