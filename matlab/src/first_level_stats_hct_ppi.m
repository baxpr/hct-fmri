function first_level_stats_hct_ppi(inp)

% Block design, four predictors: anticipate, heart, sun, fixation
% Some 5-sec rest sections are left out in the model (motion is expected)
% Contrasts of interest will be
%      heart vs fixation
%      sun vs fixation
%      heart vs sun
%
% PPI regressors added for the supplied VOI, for heart and sun conditions.
% Probably only the heart vs sun PPI contrast is interpretable, due to
% expected movement in the rest sections affecting the baseline
% connectivity estimate.
%
% VOIs are created per session. Data can be adjusted, e.g. for movement
% params - note that "effects of interest" is the contrast specified here,
% not "effects of no interest".


%% DRAFT %%%%%%%%%%

matlabbatch{1}.spm.stats.con.spmmat = {'spm_hct/SPM.mat'};
matlabbatch{1}.spm.stats.con.consess{1}.fcon.name = 'Effects of Interest';
matlabbatch{1}.spm.stats.con.consess{1}.fcon.weights = [1 0 0 0
                                                        0 1 0 0
                                                        0 0 1 0
                                                        0 0 0 1];
matlabbatch{1}.spm.stats.con.consess{1}.fcon.sessrep = 'sess';
matlabbatch{1}.spm.stats.con.delete = 0;

matlabbatch{2}.spm.util.voi.spmmat = {'spm_hct/SPM.mat'};
matlabbatch{2}.spm.util.voi.adjust = 9;
matlabbatch{2}.spm.util.voi.session = 1;
matlabbatch{2}.spm.util.voi.name = 'testvoi';
matlabbatch{2}.spm.util.voi.roi{1}.label.image = {'atlas-BAISins_space-MNI152NLin6Asym_res-02_dseg.nii'};
matlabbatch{2}.spm.util.voi.roi{1}.label.list = 1;
matlabbatch{2}.spm.util.voi.expression = 'i1';

matlabbatch{1}.spm.stats.ppi.spmmat = {'spm_hct/SPM.mat'};
matlabbatch{1}.spm.stats.ppi.type.ppi.voi = {'spm_hct/VOI_testvoi_1.mat'};
matlabbatch{1}.spm.stats.ppi.type.ppi.u = [2 1 1];  % For heart
matlabbatch{1}.spm.stats.ppi.name = 'testvoi_2_1_1';
matlabbatch{1}.spm.stats.ppi.disp = 1;

% At this point, what we have is the two PPI time series for this ROI
% (after running PPI once for each condition).

% Andy advises concatenating runs:
% https://andysbrainbook.readthedocs.io/en/latest/SPM/SPM_Short_Course/SPM_PPI.html
%
% See also https://www.fil.ion.ucl.ac.uk/spm/docs/manual/ppi/ppi/
%
% Note that concatenating movement regressors means they will not be
% completely removed, as weights may vary between runs. Better to center
% these and put in their own sections per run?
% 
% kron is used to get predictors for the means though.
%
% Note also that Andy includes both conditions in a single PPI with a 1 -1
% contrast, e.g.
%    ppi.u = [2 1 1; 3 1 -1]  % Heart gt Counting
%
% Instead, I think we'll use this procedure to generate adjusted PPI
% covariates per session, then go back and insert them in the original
% analysis config as two additional regressors. Then we'll need the
% additional contrasts for them as well.


%% END DRAFT %%%%%%%%%%%%%%%%



tag = ['hctppi-' inp.voiname];

% Filter param
hpf_sec = str2double(inp.hpf_sec);

% Select available runs
runs =[];
if ~strcmp(inp.fmri1_nii,'NONE'), runs = [runs 1]; end
if ~strcmp(inp.fmri2_nii,'NONE'), runs = [runs 2]; end
if ~strcmp(inp.fmri3_nii,'NONE'), runs = [runs 3]; end
if ~strcmp(inp.fmri4_nii,'NONE'), runs = [runs 4]; end


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

rct = 0;
for r = runs

    rct = rct + 1;

	% Session-specific scans, regressors, params
	matlabbatch{1}.spm.stats.fmri_spec.sess(rct).scans = ...
		{inp.(['swfmri' num2str(r) '_nii'])};
	matlabbatch{1}.spm.stats.fmri_spec.sess(rct).multi = {''};
	matlabbatch{1}.spm.stats.fmri_spec.sess(rct).regress = ...
		struct('name', {}, 'val', {});
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

