function first_level_stats_hctppi(inp)

% Block design, four predictors: anticipate, heart, sun, fixation
% Some 5-sec rest sections are left out in the model (motion is expected)
% Contrasts of interest will be
%      heart vs fixation
%      sun vs fixation
%      heart vs sun
%
% PPI regressors are added
%    inp.voi_name      % VOI name
%    inp.ppi_dir       % Where to find PPI files
%    inp.ppi_con       % 'Heart_gt_Counting'

% FIXME
% Homemade condition-specific PPI regressors
%  PPI.xn is the supersampled time series
%  PPI.psy.u is the boxcar stimulus function (need conv/shift?)
%  Then how to resample back to TR?
% Or, simpler,
%  PPI.Y is the VOI time series (or Y in the VOI file)
%  SPM.xX.X has condition vectors that are unfiltered and can be applied directly
%  SPM.Sess(s).row,col gives us the elements to extract from SPM.xX.X


%% FIXME FIXME
%
% Scaling the connectivity regressors? Shouldn't units already match across
% subjects due to the scaling of the original data, assuming the PPI isn't
% weirdly rescaled? Connectivity seems like it might call for z score
% scaling, but that would have to be applied to the data and the design
% matrix both to really make sense.
%
% Let's run both ways and review the regressor amplitudes, correlation
% between the two PPI contrasts in the desmtx, and con images for the two
% contrasts to be sure things are making sense.


%%
%heartPPI = load('../../OUTPUTS/ppi_hct/PPI_Heart_dAI_L_sess1.mat','PPI');
%load('../../OUTPUTS/spm_hct/SPM.mat','SPM');

% Get the task predictor for comparison
%blockHeart = SPM.SPM.xX.X(SPM.SPM.Sess(1).row,SPM.SPM.Sess(1).col(2));

% Use SPM's PPI boxcar and supersampled deconvolved time series to get the
% deconvolved PPI for heart periods only
%heart = heartPPI.PPI.xn .* heartPPI.PPI.psy.u;

% Convolve with HRF and trim excess at end
%heart = conv(full(heart),spm_hrf(SPM.SPM.xBF.dt));
%heart = heart(1:numel(heartPPI.PPI.xn));

% Downsample back to scan temporal resolution
%heart = heart(SPM.SPM.xBF.T0:SPM.SPM.xBF.T:end);

%figure(1); clf
%plot([blockHeart zscore(heart)])
%corr([blockHeart zscore(heart)])



%%







tag = ['hctppi_' inp.voi_name];

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
N = nifti(inp.(['fmri' num2str(runs(1)) '_nii']));
tr = N.timing.tspace;
for r = runs(2:end)
	N = nifti(inp.(['fmri' num2str(r) '_nii']));
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
    %ppi_file = [inp.ppi_dir filesep 'PPI_' inp.ppi_con '_' inp.voi_name '_sess' num2str(rct) '.mat'];
    %load(ppi_file,'PPI');
    ppi1_file = [inp.ppi_dir filesep 'PPI_Heart_' inp.voi_name '_sess' num2str(rct) '.mat'];
    ppi1 = load(ppi1_file,'PPI');
    ppi2_file = [inp.ppi_dir filesep 'PPI_Counting_' inp.voi_name '_sess' num2str(rct) '.mat'];
    ppi2 = load(ppi2_file,'PPI');

    % Get the task-specific connectivity regressors using SPM PPI as
    % starting point
    heart = ppi1.PPI.xn .* ppi1.PPI.psy.u;
    heart = conv(full(heart),spm_hrf(SPM.xBF.dt));
    heart = heart(1:numel(ppi1.PPI.xn));
    heart = heart(SPM.xBF.T0:SPM.xBF.T:end);

    counting = ppi2.PPI.xn .* ppi2.PPI.psy.u;
    counting = conv(full(counting),spm_hrf(SPM.xBF.dt));
    counting = counting(1:numel(ppi2.PPI.xn));
    counting = counting(SPM.xBF.T0:SPM.xBF.T:end);

	% Session-specific scans, regressors, params
	matlabbatch{1}.spm.stats.fmri_spec.sess(rct).scans = ...
		cellstr(spm_select('expand',inp.(['fmri' num2str(r) '_nii'])));
	matlabbatch{1}.spm.stats.fmri_spec.sess(rct).multi = {''};
	%matlabbatch{1}.spm.stats.fmri_spec.sess(rct).regress(1) = ...
	%	struct('name', {ppi1.PPI.xY.name}, 'val', {ppi1.PPI.Y});
	%matlabbatch{1}.spm.stats.fmri_spec.sess(rct).regress(2) = ...
	%	struct('name', {ppi1.PPI.name}, 'val', {ppi1.PPI.ppi});
	%matlabbatch{1}.spm.stats.fmri_spec.sess(rct).regress(3) = ...
	%	struct('name', {ppi2.PPI.name}, 'val', {ppi2.PPI.ppi});
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

% Effects of interest contrasts per session to use with PPI generation
c = c + 1;
matlabbatch{3}.spm.stats.con.consess{c}.fcon.name = 'Effects of Interest';
matlabbatch{3}.spm.stats.con.consess{c}.fcon.weights = [1 0 0 0
                                                        0 1 0 0
                                                        0 0 1 0
                                                        0 0 0 1];
matlabbatch{3}.spm.stats.con.consess{c}.fcon.sessrep = 'sess';


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

