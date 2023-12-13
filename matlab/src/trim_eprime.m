warning('off','MATLAB:table:ModifiedAndSavedVarnames')
q = readtable('../../OUTPUTS/eprime.csv');

trimq = q(:,{
    'Level'
    'Procedure'
    'Fixation1_OnsetTime'
    'Fixation1_Duration'
    'Fixation17_OnsetTime'
    'Fixation17_Duration'
    'ImageDisplay1_OnsetTime'
    'ImageDisplay1_Duration'
    'FlashImage_OnsetTime'
    'FlashImage_Duration'
    'Fixation2_OnsetTime'
    'Fixation2_Duration'
    'Fixation10_OnsetTime'
    'Fixation10_Duration'
    'Instruction_OffsetTime'
    'WaitForScanner1_OffsetTime'
    'WaitForScanner2_OffsetTime'
    'WaitForScanner3_OffsetTime'
    });


% Run 1

run1 = table(cell(0,1),[],[],'VariableNames',{'condition','onset','duration'});

% heart fixation
inds = find(~isnan(q.Fixation1_OnsetTime));
ons = q.Fixation1_OnsetTime(inds);
dur = q.Fixation1_Duration(inds);
for k = 1:numel(inds)
    run1.condition{end+1,1} = 'fixation';
    run1.onset(end,1) = ons(k) / 1000;
    run1.duration(end,1) = dur(k) / 1000;
end

% heart anticipate
inds = find(~isnan(q.Fixation17_OnsetTime));
ons = q.Fixation17_OnsetTime(inds);
dur = q.Fixation17_Duration(inds);
for k = 1:numel(inds)
    run1.condition{end+1,1} = 'anticipate';
    run1.onset(end,1) = ons(k) / 1000;
    run1.duration(end,1) = dur(k) / 1000;
end

% heart
inds = find(~isnan(q.ImageDisplay1_OnsetTime));
ons = q.ImageDisplay1_OnsetTime(inds);
dur = q.ImageDisplay1_Duration(inds);
for k = 1:numel(inds)
    run1.condition{end+1,1} = 'heart';
    run1.onset(end,1) = ons(k) / 1000;
    run1.duration(end,1) = dur(k) / 1000;
end

% heart response
inds = find(~isnan(q.RESPONSE_OnsetTime));
ons = q.RESPONSE_OnsetTime(inds);
dur = q.RESPONSE_Duration(inds);
for k = 1:numel(inds)
    run1.condition{end+1,1} = 'response';
    run1.onset(end,1) = ons(k) / 1000;
    run1.duration(end,1) = dur(k) / 1000;
end

% counting fixation
inds = find(~isnan(q.Fixation2_OnsetTime));
ons = q.Fixation2_OnsetTime(inds);
dur = q.Fixation2_Duration(inds);
for k = 1:numel(inds)
    run1.condition{end+1,1} = 'fixation';
    run1.onset(end,1) = ons(k) / 1000;
    run1.duration(end,1) = dur(k) / 1000;
end

% counting anticipate
inds = find(~isnan(q.Fixation10_OnsetTime));
ons = q.Fixation10_OnsetTime(inds);
dur = q.Fixation10_Duration(inds);
for k = 1:numel(inds)
    run1.condition{end+1,1} = 'anticipate';
    run1.onset(end,1) = ons(k) / 1000;
    run1.duration(end,1) = dur(k) / 1000;
end

% counting (use previous inds from anticipate condition)
ons = q.Fixation10_OnsetTime(inds) + q.Fixation10_Duration(inds);
dur = q.RESPONSE1_OnsetTime(inds) - ons;
for k = 1:numel(inds)
    run1.condition{end+1,1} = 'counting';
    run1.onset(end,1) = ons(k) / 1000;
    run1.duration(end,1) = dur(k) / 1000;
end

% counting response
inds = find(~isnan(q.RESPONSE1_OnsetTime));
ons = q.RESPONSE1_OnsetTime(inds);
dur = q.RESPONSE1_Duration(inds);
for k = 1:numel(inds)
    run1.condition{end+1,1} = 'response';
    run1.onset(end,1) = ons(k) / 1000;
    run1.duration(end,1) = dur(k) / 1000;
end

sortrows(run1,'onset')


% Heart block is ImageDisplay1, or offset of Fixation17 to onset of
%    RESPONSE. 30 sec.
% Counting block is offset of Fixation10 to onset of RESPONSE1. 28 sec.

% Plain fixation is Fixation1, Fixation2  (13 sec)

% Anticipatory fixation is Fixation17, Fixation10 (2 sec)

% Response block is RESPONSE, RESPONSE1 (5 sec)


% Procedure: Heart
%   Fixation1
%   Fixation17
%   ImageDisplay1
%   RESPONSE

% Procedure: flash

% Procedure: Counting
%   Fixation2
%   Fixation10
%   RESPONSE1
