function timings = get_timings(eprime_csv)

warning('off','MATLAB:table:ModifiedAndSavedVarnames')
eprime = readtable(eprime_csv);

% Run 1

% First grab the start time
run1 = table( ...
    {'Instruction'}, ...
    {'scanstart'}, ...
    eprime.Instruction_RTTime(~isnan(eprime.Instruction_RTTime)) / 1000, ...
    0, ...
    'VariableNames', ...
    {'eprime_label','condition','onset_sec','duration_sec'} ...
    );

% Conditions
run1 = [run1; parsefun(eprime,'Fixation1','fixation')];
run1 = [run1; parsefun(eprime,'Fixation17','anticipate')];
run1 = [run1; parsefun(eprime,'ImageDisplay1','heart')];
run1 = [run1; parsefun(eprime,'RESPONSE','response')];
run1 = [run1; parsefun(eprime,'Fixation2','fixation')];
run1 = [run1; parsefun(eprime,'Fixation10','anticipate')];
run1 = [run1; parsefun(eprime,'RESPONSE1','response')];

% Relative onsets/offsets for fmri
run1.fmri_onset_sec = run1.onset_sec ...
    - run1.onset_sec(strcmp(run1.condition,'scanstart'));
run1.fmri_offset_sec = run1.fmri_onset_sec + run1.duration_sec;

% Run label
run1.run(:) = 1;

run1 = sortrows(run1,'fmri_onset_sec');


% Final return value
timings{1} = run1;


