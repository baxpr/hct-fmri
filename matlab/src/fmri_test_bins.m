% Test data
trials = readtable('../../INPUTS/fmri_test_bins_trial_report.csv');


% Bin resolution (?)
TR = 1.3;
MR = 40;
bin = TR/MR;

% Switch action vector is 2 for switch, 1 for not-switch, 0 for first trial (?)
swi = nan(height(trials),1);
swi(strcmp(trials.Switch,'Switch')) = 2;
swi(strcmp(trials.Switch,'Stay')) = 1;
swi(strcmp(trials.Switch,'InitialTrial')) = 0;

% DEFINE fMRI regressors for fMRI analysis from hgf fit
epsi23   = [trials.traj_epsi_2 trials.traj_epsi_3];
psi23    = [trials.traj_psi_2 trials.traj_psi_3];
mu3      = trials.traj_muhat_31;
regs          = [epsi23 psi23 mu3 swi];

% get onsets
names            = {'trial', 'missings'};
durations        = {0 0};
onsets_cue       = trials.T1_TrialStart_fMRIsec;
onsets_feedback  = trials.T3_FeedbackOnset_fMRIsec;
offsets_feedback = onsets_feedback + 1;
onsets_missings  = onsets_feedback(strcmp(trials.TrialType,'NoResponse'));
onsets_cue       = onsets_cue(~strcmp(trials.TrialType,'NoResponse'));
onsets_feedback  = onsets_feedback(~strcmp(trials.TrialType,'NoResponse'));
offsets_feedback = offsets_feedback(~strcmp(trials.TrialType,'NoResponse'));

trial = sum([onsets_feedback-onsets_cue offsets_feedback-onsets_feedback],2);
events_c = ceil(trial/bin);
events_f = floor(trial/bin);
diff_c = abs(trial-events_c*bin);
diff_f = abs(trial-events_f*bin);
if diff_c < diff_f
	events = events_c;
else
	events = events_f;
end
nPE = 9;
nPR = events-nPE;
ind = [0; cumsum(events(1:end))];
for i=1:length(onsets_cue)
	ons(ind(i)+1:ind(i+1),1)=[onsets_cue(i):bin:onsets_cue(i)+(events(i)-1)*bin]';
	ons(ind(i)+1:ind(i+1),2)=[ones(nPR(i),1); ones(nPE,1)+1];
	PEs(ind(i)+1:ind(i+1),1:3)=[zeros(nPR(i),3); repmat(regs(i,1:3), nPE,1)];
	PEs(ind(i)+1:ind(i+1),1)=[zeros(nPR(i),1); repmat(regs(i,1), nPE,1)];
	if     i< length(onsets_cue)
		PRs(ind(i)+1:ind(i+1))=[repmat(regs(i), nPR(i),1); repmat(regs(i+1), nPE,1)];
		PRs(ind(i)+1:ind(i+1),1:4)=[repmat(regs(i,3:end), nPR(i),1); repmat(regs(i+1,3:end), nPE,1)];
	elseif i==length(onsets_cue)
		PRs(ind(i)+1:ind(i+1))=[repmat(regs(i), nPR(i),1); repmat(regs(i), nPE,1)];
		PRs(ind(i)+1:ind(i+1),1:4)=[repmat(regs(i,3:end), nPR(i),1); repmat(regs(i,3:end), nPE,1)];
	end
end

onsets{1}       = ons(:,1);
onsets{2}       = onsets_missings;
if isempty(onsets{2}); onsets{2} = NaN; end

%get parametric modulators from model
pmod(1).name{1} = 'E2';
pmod(1).param{1}= [PEs(:,1)];
pmod(1).poly{1} = 1;
pmod(1).name{3} = 'E3';
pmod(1).param{3}= [PEs(:,2)];
pmod(1).poly{3} = 1;
pmod(1).name{4} = 'PSI2';
pmod(1).param{4}= [PRs(:,1)];
pmod(1).poly{4} = 1;
pmod(1).name{5} = 'PSI3';
pmod(1).param{5}= [PRs(:,2)];
pmod(1).poly{5} = 1;
pmod(1).name{6} = 'MU3';
pmod(1).param{6}= [PRs(:,3)];
pmod(1).poly{6} = 1;
pmod(1).name{7} = 'switch';
pmod(1).param{7}= [PRs(:,4)];
pmod(1).poly{7} = 1;
%pmod(1).name{1} = 'da';
%pmod(1).param{1}= [PEs(:,1)];
%pmod(1).poly{1} = 1;
