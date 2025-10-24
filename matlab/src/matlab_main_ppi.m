function matlab_main_ppi(inp)

spm('defaults','fmri')

fwhm_mm = str2double(inp.fwhm_mm);

% Mean image
clear job
job.comp{1}.def = {inp.deffwd_nii};
job.comp{2}.id.space = {which(inp.refimg_nii)};
job.out{1}.pull.fnames = {inp.meanfmri_nii};
job.out{1}.pull.savedir.saveusr = {inp.out_dir};
job.out{1}.pull.interp = 1;
job.out{1}.pull.mask = 0;
job.out{1}.pull.fwhm = [0 0 0];
spm_deformations(job);
[~,n,e] = fileparts(inp.meanfmri_nii);
inp.wmeanfmri_nii = fullfile(inp.out_dir,['w' n e]);


% Full runs
for fmri_nii = {'fmri1_nii','fmri2_nii','fmri3_nii','fmri4_nii'}

    if strcmp(inp.(fmri_nii{1}),'NONE')
        inp.(['w' fmri_nii{1}]) = 'NONE';
        inp.(['sw' fmri_nii{1}]) = 'NONE';
    else

        % Warp
        clear job
        job.comp{1}.def = {inp.deffwd_nii};
        job.comp{2}.id.space = {which(inp.refimg_nii)};
        job.out{1}.pull.fnames = {inp.(fmri_nii{1})};
        job.out{1}.pull.savedir.saveusr = {inp.out_dir};
        job.out{1}.pull.interp = 1;
        job.out{1}.pull.mask = 0;
        job.out{1}.pull.fwhm = [0 0 0];
        spm_deformations(job);

        % Get filename of warped image
        [~,n,e] = fileparts(inp.(fmri_nii{1}));
        inp.(['w' fmri_nii{1}]) = fullfile(inp.out_dir,['w' n e]);

        % Smooth
        clear matlabbatch
        matlabbatch{1}.spm.spatial.smooth.data = {inp.(['w' fmri_nii{1}])};
        matlabbatch{1}.spm.spatial.smooth.fwhm = [fwhm_mm fwhm_mm fwhm_mm];
        matlabbatch{1}.spm.spatial.smooth.dtype = 0;
        matlabbatch{1}.spm.spatial.smooth.im = 0;
        matlabbatch{1}.spm.spatial.smooth.prefix = 's';
        spm_jobman('run',matlabbatch);

        % Get filename of smoothed warped image
        [~,n,e] = fileparts(inp.(['w' fmri_nii{1}]));
        inp.(['sw' fmri_nii{1}]) = fullfile(inp.out_dir,['s' n e]);

    end

end

% First level stats and contrasts, no PPI yet. We use unsmoothed images so
% that ROIs are as specified, not larger due to smoothing
spm_dir = first_level_stats_hct_unsmoothed(inp);

% Set up more inputs for PPI stuff (ROI set is hardcoded)
inp.spm_dir = spm_dir;
inp.ppi_dir = fullfile(inp.out_dir,'ppi');
if ~exist(inp.ppi_dir,'dir'), mkdir(inp.ppi_dir), end
inp.voi_label_image = 'atlas-BAISins_space-MNI152NLin6Asym_res-02_dseg.nii';

% Get ROI list
roilist = readtable('atlas-BAISins_dseg.tsv','FileType','text');

% VOI and PPIs for each ROI
for roi = roilist.index
    inp.voi_label_index = roi;
    inp.voi_name = roilist.label{roilist.index==roi};
    create_ppi_regressors(inp);
end

% First level stats with PPI for each ROI
for roi = roilist.index
    inp.voi_name = roilist.label{roilist.index==roi};
    first_level_stats_hct_manppi(inp)
end

