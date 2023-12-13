function matlab_main(inp)

spm('defaults','fmri')

% Apply cat12 warp
disp('Warp')
clear job
job.comp{1}.def = {inp.deffwd_nii};
job.comp{2}.id.space = {which(inp.refimg_nii)};
job.out{1}.pull.fnames = {
        inp.meanfmri_nii
        inp.fmri1_nii
        inp.fmri2_nii
        inp.fmri3_nii
        inp.fmri4_nii
        };
job.out{1}.pull.savedir.saveusr = {inp.out_dir};
job.out{1}.pull.interp = 1;
job.out{1}.pull.mask = 0;
job.out{1}.pull.fwhm = [0 0 0];
spm_deformations(job);

% Get filenames of warped images
[~,n,e] = fileparts(inp.fmri1_nii);
inp.wfmri1_nii = fullfile(inp.out_dir,['w' n e]);
[~,n,e] = fileparts(inp.fmri2_nii);
inp.wfmri2_nii = fullfile(inp.out_dir,['w' n e]);
[~,n,e] = fileparts(inp.fmri3_nii);
inp.wfmri3_nii = fullfile(inp.out_dir,['w' n e]);
[~,n,e] = fileparts(inp.fmri4_nii);
inp.wfmri4_nii = fullfile(inp.out_dir,['w' n e]);
[~,n,e] = fileparts(inp.meanfmri_nii);
inp.wmeanfmri_nii = fullfile(inp.out_dir,['w' n e]);

% Smooth warped fmri timeseries
disp('Smoothing')
fwhm_mm = str2double(inp.fwhm_mm);
clear matlabbatch
matlabbatch{1}.spm.spatial.smooth.data = {
	inp.wfmri1_nii
	inp.wfmri2_nii
	inp.wfmri3_nii
	inp.wfmri4_nii
	};
matlabbatch{1}.spm.spatial.smooth.fwhm = [fwhm_mm fwhm_mm fwhm_mm];
matlabbatch{1}.spm.spatial.smooth.dtype = 0;
matlabbatch{1}.spm.spatial.smooth.im = 0;
matlabbatch{1}.spm.spatial.smooth.prefix = 's';
spm_jobman('run',matlabbatch);

% Get filenames of smoothed warped images
[~,n,e] = fileparts(inp.wfmri1_nii);
inp.swfmri1_nii = fullfile(inp.out_dir,['s' n e]);
[~,n,e] = fileparts(inp.wfmri2_nii);
inp.swfmri2_nii = fullfile(inp.out_dir,['s' n e]);
[~,n,e] = fileparts(inp.wfmri3_nii);
inp.swfmri3_nii = fullfile(inp.out_dir,['s' n e]);
[~,n,e] = fileparts(inp.wfmri4_nii);
inp.swfmri4_nii = fullfile(inp.out_dir,['s' n e]);

% First level stats and contrasts
first_level_stats_hct(inp);
