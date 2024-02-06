function timings = get_timings(eprime_csv)

warning('off','MATLAB:table:ModifiedAndSavedVarnames')
warning('off','MATLAB:table:RowsAddedExistingVars')
eprime = readtable(eprime_csv);

%% Run 1

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

% Temporal order by onset
run1 = sortrows(run1,'fmri_onset_sec');


%% Run 2
run2 = table( ...
    {'WaitForScanner1'}, ...
    {'scanstart'}, ...
    eprime.WaitForScanner1_RTTime(~isnan(eprime.WaitForScanner1_RTTime)) / 1000, ...
    0, ...
    'VariableNames', ...
    {'eprime_label','condition','onset_sec','duration_sec'} ...
    );
run2 = [run2; parsefun(eprime,'Fixation3','fixation')];
run2 = [run2; parsefun(eprime,'Fixation12','anticipate')];
run2 = [run2; parsefun(eprime,'ImageDisplay6','heart')];
run2 = [run2; parsefun(eprime,'RESPONSE2','response')];
run2 = [run2; parsefun(eprime,'Fixation4','fixation')];
run2 = [run2; parsefun(eprime,'Fixation11','anticipate')];
run2 = [run2; parsefun(eprime,'RESPONSE3','response')];
run2.fmri_onset_sec = run2.onset_sec ...
    - run2.onset_sec(strcmp(run2.condition,'scanstart'));
run2.fmri_offset_sec = run2.fmri_onset_sec + run2.duration_sec;
run2.run(:) = 2;
run2 = sortrows(run2,'fmri_onset_sec');


%% Run 3
run3 = table( ...
    {'WaitForScanner2'}, ...
    {'scanstart'}, ...
    eprime.WaitForScanner2_RTTime(~isnan(eprime.WaitForScanner2_RTTime)) / 1000, ...
    0, ...
    'VariableNames', ...
    {'eprime_label','condition','onset_sec','duration_sec'} ...
    );
run3 = [run3; parsefun(eprime,'Fixation5','fixation')];
run3 = [run3; parsefun(eprime,'Fixation13','anticipate')];
run3 = [run3; parsefun(eprime,'ImageDisplay9','heart')];
run3 = [run3; parsefun(eprime,'RESPONSE5','response')];
run3 = [run3; parsefun(eprime,'Fixation6','fixation')];
run3 = [run3; parsefun(eprime,'Fixation14','anticipate')];
run3 = [run3; parsefun(eprime,'RESPONSE4','response')];
run3.fmri_onset_sec = run3.onset_sec ...
    - run3.onset_sec(strcmp(run3.condition,'scanstart'));
run3.fmri_offset_sec = run3.fmri_onset_sec + run3.duration_sec;
run3.run(:) = 3;
run3 = sortrows(run3,'fmri_onset_sec');


%% Run 4
run4 = table( ...
    {'WaitForScanner3'}, ...
    {'scanstart'}, ...
    eprime.WaitForScanner3_RTTime(~isnan(eprime.WaitForScanner3_RTTime)) / 1000, ...
    0, ...
    'VariableNames', ...
    {'eprime_label','condition','onset_sec','duration_sec'} ...
    );
run4 = [run4; parsefun(eprime,'Fixation7','fixation')];
run4 = [run4; parsefun(eprime,'Fixation15','anticipate')];
run4 = [run4; parsefun(eprime,'ImageDisplay12','heart')];
run4 = [run4; parsefun(eprime,'RESPONSE6','response')];
run4 = [run4; parsefun(eprime,'Fixation8','fixation')];
run4 = [run4; parsefun(eprime,'Fixation16','anticipate')];
run4 = [run4; parsefun(eprime,'RESPONSE7','response')];
run4.fmri_onset_sec = run4.onset_sec ...
    - run4.onset_sec(strcmp(run4.condition,'scanstart'));
run4.fmri_offset_sec = run4.fmri_onset_sec + run4.duration_sec;
run4.run(:) = 4;
run4 = sortrows(run4,'fmri_onset_sec');


%% Final return values
timings{1} = run1;
timings{2} = run2;
timings{3} = run3;
timings{4} = run4;



