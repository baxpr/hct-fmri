function first_level_stats_hct(inp)

% Block design, four predictors: anticipate, heart, sun, fixation
% Some 5-sec rest sections are left out in the model (motion is expected)
% Contrasts of interest will be
%      heart vs fixation
%      sun vs fixation
%      heart vs sun

tag = 'hct';

% Filter param
hpf_sec = str2double(inp.hpf_sec);

% Save motion params as .mat
for r = 1:4
	mot = readtable(inp.(['motpar' num2str(r) '_txt']),'FileType','text');
	mot = zscore(table2array(mot(:,1:6)));
	writematrix(mot, fullfile(inp.out_dir,['motpar' num2str(r) '.txt']))
end

% Get TRs and check
N = nifti(inp.swfmri1_nii);
tr = N.timing.tspace;
for r = 2:4
	N = nifti(inp.(['swfmri' num2str(r) '_nii']));
	if abs(N.timing.tspace-tr) > 0.001
		error('TR not matching for run %d',r)
	end
end
fprintf('ALERT: USING TR OF %0.3f sec FROM FMRI NIFTI\n',tr)

% Load condition timing info
timings = get_timings(inp.eprime_csv);


%% Design
clear matlabbatch
matlabbatch{1}.spm.stats.fmri_spec.dir = ...
	{fullfile(inp.out_dir,['spm_' tag])};
matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
matlabbatch{1}.spm.stats.fmri_spec.timing.RT = tr;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 1;
matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
matlabbatch{1}.spm.stats.fmri_spec.mthresh = -Inf;
matlabbatch{1}.spm.stats.fmri_spec.mask = {[spm('dir') '/tpm/mask_ICV.nii']};
matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';

for r = 1:4

	% Session-specific scans, regressors, params
	matlabbatch{1}.spm.stats.fmri_spec.sess(r).scans = ...
		{inp.(['swfmri' num2str(r) '_nii'])};
	matlabbatch{1}.spm.stats.fmri_spec.sess(r).multi = {''};
	matlabbatch{1}.spm.stats.fmri_spec.sess(r).regress = ...
		struct('name', {}, 'val', {});
	matlabbatch{1}.spm.stats.fmri_spec.sess(r).multi_reg = ...
		{fullfile(inp.out_dir,['motpar' num2str(r) '.txt'])};
	matlabbatch{1}.spm.stats.fmri_spec.sess(r).hpf = hpf_sec;
	
    % Conditions
    c = 0;
    for cond = {'anticipate','heart','counting','fixation'}
    	c = c + 1;
    	matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(c).name = cond{1};
    	matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(c).onset = ...
    		timings{r}.fmri_onset_sec(strcmp(timings{r}.condition,cond));
    	matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(c).duration = ...
    		timings{r}.duration_sec(strcmp(timings{r}.condition,cond));
    	matlabbatch{1}.spm.stats.fmri_spec.sess(r).cond(c).tmod = 0;
    end

end


%% Estimate
matlabbatch{2}.spm.stats.fmri_est.spmmat = ...
	fullfile(matlabbatch{1}.spm.stats.fmri_spec.dir,'SPM.mat');
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;


%% Contrasts
matlabbatch{3}.spm.stats.con.spmmat = ...
	matlabbatch{2}.spm.stats.fmri_est.spmmat;
matlabbatch{3}.spm.stats.con.delete = 1;
c = 0;

c = c + 1;
matlabbatch{3}.spm.stats.con.consess{c}.tcon.name = 'Anticipate gt Fixation';
matlabbatch{3}.spm.stats.con.consess{c}.tcon.weights = [1 0 0 -1];
matlabbatch{3}.spm.stats.con.consess{c}.tcon.sessrep = 'replsc';

c = c + 1;
matlabbatch{3}.spm.stats.con.consess{c}.tcon.name = 'Heart gt Fixation';
matlabbatch{3}.spm.stats.con.consess{c}.tcon.weights = [0 1 0 -1];
matlabbatch{3}.spm.stats.con.consess{c}.tcon.sessrep = 'replsc';

c = c + 1;
matlabbatch{3}.spm.stats.con.consess{c}.tcon.name = 'Counting gt Fixation';
matlabbatch{3}.spm.stats.con.consess{c}.tcon.weights = [0 0 1 -1];
matlabbatch{3}.spm.stats.con.consess{c}.tcon.sessrep = 'replsc';

c = c + 1;
matlabbatch{3}.spm.stats.con.consess{c}.tcon.name = 'Heart gt Counting';
matlabbatch{3}.spm.stats.con.consess{c}.tcon.weights = [0 1 -1 0];
matlabbatch{3}.spm.stats.con.consess{c}.tcon.sessrep = 'replsc';

% Inverse of all existing contrasts since SPM won't show us both sides
numc = numel(matlabbatch{3}.spm.stats.con.consess);
for k = 1:numc
        c = c + 1;
        matlabbatch{3}.spm.stats.con.consess{c}.tcon.name = ...
                ['Neg ' matlabbatch{3}.spm.stats.con.consess{c-numc}.tcon.name];
        matlabbatch{3}.spm.stats.con.consess{c}.tcon.weights = ...
                - matlabbatch{3}.spm.stats.con.consess{c-numc}.tcon.weights;
        matlabbatch{3}.spm.stats.con.consess{c}.tcon.sessrep = 'replsc';
end


%% Review design
matlabbatch{4}.spm.stats.review.spmmat = ...
	matlabbatch{2}.spm.stats.fmri_est.spmmat;
matlabbatch{4}.spm.stats.review.display.matrix = 1;
matlabbatch{4}.spm.stats.review.print = false;

matlabbatch{5}.cfg_basicio.run_ops.call_matlab.inputs{1}.string = ...
        fullfile(inp.out_dir,['first_level_design_' tag '.png']);
matlabbatch{5}.cfg_basicio.run_ops.call_matlab.outputs = cell(1,0);
matlabbatch{5}.cfg_basicio.run_ops.call_matlab.fun = 'spm_window_print';


%% Save batch and run
save(fullfile(inp.out_dir,['spmbatch_first_level_stats_' tag '.mat']),'matlabbatch')
spm_jobman('run',matlabbatch);

% And save contrast names
numc = numel(matlabbatch{3}.spm.stats.con.consess);
connames = table((1:numc)','VariableNames',{'ConNum'});
for k = 1:numc
	try
		connames.ConName{k,1} = ...
			matlabbatch{3}.spm.stats.con.consess{k}.tcon.name;
	catch
		connames.ConName{k,1} = ...
			matlabbatch{3}.spm.stats.con.consess{k}.fcon.name;
	end
end
writetable(connames,fullfile(inp.out_dir,['spm_contrast_names_' tag '.csv']));


%% Results display
% Needed to create the spmT even if we don't get the figure window
xSPM = struct( ...
    'swd', matlabbatch{1}.spm.stats.fmri_spec.dir, ...
    'title', '', ...
    'Ic', 2, ...
    'n', 0, ...
    'Im', [], ...
    'pm', [], ...
    'Ex', [], ...
    'u', 0.005, ...
    'k', 10, ...
    'thresDesc', 'none' ...
    );
[hReg,xSPM] = spm_results_ui('Setup',xSPM);

% Show on the subject MNI anat
spm_sections(xSPM,hReg,inp.biasnorm_nii)

% Jump to global max activation
spm_mip_ui('Jump',spm_mip_ui('FindMIPax'),'glmax');

