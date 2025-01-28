function matlab_main(inp)

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

% First level stats and contrasts
first_level_stats_hct(inp);

