function create_ppi_regressors(inp)

% Create PPI regressors per run for one VOI. inp must contain
%   voi_label_image e.g.  'atlas-BAISins_space-MNI152NLin6Asym_res-02_dseg.nii'
%   voi_label_index e.g.  1
%   voi_name        e.g.  dAI_L
%
% We also need the first first level stats dir
%   spm_dir
%
% And a place to capture all the VOI, PPI files
%   ppi_dir

% Effects of interest contrasts have been added to orig first level run.
% They will be one per session labeled
%    Effects of Interest - Session N
% etc in SPM.xCon(N).name

% Find the SPM.mat and load it for some info
spm_mat = fullfile(inp.spm_dir,'SPM.mat');
load(spm_mat,'SPM');

% Loop over runs. These are numbered by SPM's count, NOT by the original
% HCT run number label given at scan time.
for rct = 1:numel(SPM.Sess)

    clear matlabbatch

    % Find our effects of interest for this session r
    connames = {SPM.xCon.name};
    effcon = find(strcmp(connames,['Effects of Interest - Session ' num2str(rct)]));

    % Now we load those results to start building PPI regressors
    matlabbatch{1}.spm.util.voi.spmmat = {spm_mat};
    matlabbatch{1}.spm.util.voi.session = rct;
    matlabbatch{1}.spm.util.voi.adjust = effcon;
    matlabbatch{1}.spm.util.voi.name = inp.voi_name;
    matlabbatch{1}.spm.util.voi.roi{1}.label.image = {inp.voi_label_image};
    matlabbatch{1}.spm.util.voi.roi{1}.label.list = inp.voi_label_index;
    matlabbatch{1}.spm.util.voi.expression = 'i1';

    matlabbatch{2}.spm.stats.ppi.spmmat = {spm_mat};
    matlabbatch{2}.spm.stats.ppi.type.ppi.voi = ...
        {[inp.spm_dir filesep 'VOI_' inp.voi_name '_' num2str(rct) '.mat']};
    matlabbatch{2}.spm.stats.ppi.type.ppi.u = [2 1 1];
    matlabbatch{2}.spm.stats.ppi.name = ['Heart_' inp.voi_name '_sess' num2str(rct)];
    matlabbatch{2}.spm.stats.ppi.disp = 1;

    matlabbatch{3}.spm.stats.ppi.spmmat = {spm_mat};
    matlabbatch{3}.spm.stats.ppi.type.ppi.voi = ...
        {[inp.spm_dir filesep 'VOI_' inp.voi_name '_' num2str(rct) '.mat']};
    matlabbatch{3}.spm.stats.ppi.type.ppi.u = [3 1 1];
    matlabbatch{3}.spm.stats.ppi.name = ['Counting_' inp.voi_name '_sess' num2str(rct)];
    matlabbatch{3}.spm.stats.ppi.disp = 1;

    save(fullfile(inp.ppi_dir, ...
        ['spmbatch_ppi_' inp.voi_name '_sess' num2str(rct) '.mat']),'matlabbatch')
    spm_jobman('run',matlabbatch);


end

% Copy all files to ppi_dir
copyfile([inp.spm_dir filesep 'VOI*'],inp.ppi_dir);
copyfile([inp.spm_dir filesep 'PPI*'],inp.ppi_dir);

% At this point, what we have is the PPI time series for a specific ROI.

