function first_level_stats_hct_manppi(inp)

% Re-run the HCT first level stats, but add PPI regressors for a seed ROI.
%
% Additional inputs needed:
%
%  spm_dir       Location of the original non-PPI first level stats result
%  ppi_dir       Location of VOI/PPI files generated with create_ppi_regressors
%  voi_name      Name of the VOI to use as a seed

tag = ['hct_manppi_' inp.voi_name];

% Filter param
hpf_sec = str2double(inp.hpf_sec);

% Select available runs
runs =[];
if ~strcmp(inp.swfmri1_nii,'NONE'), runs = [runs 1]; end
if ~strcmp(inp.swfmri2_nii,'NONE'), runs = [runs 2]; end
if ~strcmp(inp.swfmri3_nii,'NONE'), runs = [runs 3]; end
if ~strcmp(inp.swfmri4_nii,'NONE'), runs = [runs 4]; end


% Save motion params as .mat
for r = runs
	mot = readtable(inp.(['motpar' num2str(r) '_txt']),'FileType','text');
	mot = zscore(table2array(mot(:,1:6)));
	writematrix(mot, fullfile(inp.out_dir,['motpar' num2str(r) '.txt']))
end

% Get TRs and check
N = nifti(inp.(['swfmri' num2str(runs(1)) '_nii']));
tr = N.timing.tspace;
for r = runs(2:end)
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

% get our previous SPM analysis (no PPIs)
load(fullfile(inp.spm_dir,'SPM.mat'),'SPM');

rct = 0;
for r = runs

    rct = rct + 1;

    % Find our PPI regressors for this session (by SPM count, not original
    % scanner run labeling)
    heartppi_file = [inp.ppi_dir filesep 'PPI_Heart_' inp.voi_name '_sess' num2str(rct) '.mat'];
    heartppi = load(heartppi_file,'PPI');
    countppi_file = [inp.ppi_dir filesep 'PPI_Counting_' inp.voi_name '_sess' num2str(rct) '.mat'];
    countppi = load(countppi_file,'PPI');

    % Get the task-specific connectivity regressors using SPM PPI as
    % starting point
    heart = heartppi.PPI.xn .* heartppi.PPI.psy.u;
    heart = conv(full(heart),spm_hrf(SPM.xBF.dt));
    heart = heart(1:numel(heartppi.PPI.xn));
    heart = heart(SPM.xBF.T0:SPM.xBF.T:end);

    counting = countppi.PPI.xn .* countppi.PPI.psy.u;
    counting = conv(full(counting),spm_hrf(SPM.xBF.dt));
    counting = counting(1:numel(countppi.PPI.xn));
    counting = counting(SPM.xBF.T0:SPM.xBF.T:end);

	% Session-specific scans, regressors, params
	matlabbatch{1}.spm.stats.fmri_spec.sess(rct).scans = ...
		cellstr(spm_select('expand',inp.(['swfmri' num2str(r) '_nii'])));
	matlabbatch{1}.spm.stats.fmri_spec.sess(rct).multi = {''};
    matlabbatch{1}.spm.stats.fmri_spec.sess(rct).regress(1) = ...
        struct('name', 'HeartConn', 'val', heart);
    matlabbatch{1}.spm.stats.fmri_spec.sess(rct).regress(2) = ...
        struct('name', 'CountingConn', 'val', counting);
    matlabbatch{1}.spm.stats.fmri_spec.sess(rct).multi_reg = ...
		{fullfile(inp.out_dir,['motpar' num2str(r) '.txt'])};
	matlabbatch{1}.spm.stats.fmri_spec.sess(rct).hpf = hpf_sec;
	
    % Conditions
    c = 0;
    for cond = {'anticipate','heart','counting','fixation'}
    	c = c + 1;
    	matlabbatch{1}.spm.stats.fmri_spec.sess(rct).cond(c).name = cond{1};
    	matlabbatch{1}.spm.stats.fmri_spec.sess(rct).cond(c).onset = ...
    		timings{r}.fmri_onset_sec(strcmp(timings{r}.condition,cond));
    	matlabbatch{1}.spm.stats.fmri_spec.sess(rct).cond(c).duration = ...
    		timings{r}.duration_sec(strcmp(timings{r}.condition,cond));
    	matlabbatch{1}.spm.stats.fmri_spec.sess(rct).cond(c).tmod = 0;
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
matlabbatch{3}.spm.stats.con.consess{c}.tcon.weights = [1 0 0 -1 0 0];
matlabbatch{3}.spm.stats.con.consess{c}.tcon.sessrep = 'replsc';

c = c + 1;
matlabbatch{3}.spm.stats.con.consess{c}.tcon.name = 'Heart gt Fixation';
matlabbatch{3}.spm.stats.con.consess{c}.tcon.weights = [0 1 0 -1 0 0];
matlabbatch{3}.spm.stats.con.consess{c}.tcon.sessrep = 'replsc';

c = c + 1;
matlabbatch{3}.spm.stats.con.consess{c}.tcon.name = 'Counting gt Fixation';
matlabbatch{3}.spm.stats.con.consess{c}.tcon.weights = [0 0 1 -1 0 0];
matlabbatch{3}.spm.stats.con.consess{c}.tcon.sessrep = 'replsc';

c = c + 1;
matlabbatch{3}.spm.stats.con.consess{c}.tcon.name = 'Heart gt Counting';
matlabbatch{3}.spm.stats.con.consess{c}.tcon.weights = [0 1 -1 0 0 0];
matlabbatch{3}.spm.stats.con.consess{c}.tcon.sessrep = 'replsc';

c = c + 1;
matlabbatch{3}.spm.stats.con.consess{c}.tcon.name = ['HeartConn_' inp.voi_name];
matlabbatch{3}.spm.stats.con.consess{c}.tcon.weights = [0 0 0 0 1 0];
matlabbatch{3}.spm.stats.con.consess{c}.tcon.sessrep = 'replsc';

c = c + 1;
matlabbatch{3}.spm.stats.con.consess{c}.tcon.name = ['CountingConn_' inp.voi_name];
matlabbatch{3}.spm.stats.con.consess{c}.tcon.weights = [0 0 0 0 0 1];
matlabbatch{3}.spm.stats.con.consess{c}.tcon.sessrep = 'replsc';

c = c + 1;
matlabbatch{3}.spm.stats.con.consess{c}.tcon.name = ['HeartConn_plus_CountingConn_' inp.voi_name];
matlabbatch{3}.spm.stats.con.consess{c}.tcon.weights = [0 0 0 0 0.5 0.5];
matlabbatch{3}.spm.stats.con.consess{c}.tcon.sessrep = 'replsc';

c = c + 1;
matlabbatch{3}.spm.stats.con.consess{c}.tcon.name = ['HeartConn_gt_CountingConn_' inp.voi_name];
matlabbatch{3}.spm.stats.con.consess{c}.tcon.weights = [0 0 0 0 1 -1];
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

