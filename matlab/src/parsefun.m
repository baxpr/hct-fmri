function info = parsefun(eprime,eprime_label,newlabel)

info = table(cell(0,1),cell(0,1),[],[], ...
    'VariableNames',{'eprime_label','condition','onset_sec','duration_sec'});

inds = find(~isnan(eprime.([eprime_label '_OnsetTime'])));
ons = eprime.([eprime_label '_OnsetTime'])(inds);
dur = eprime.([eprime_label '_Duration'])(inds);
for k = 1:numel(inds)
    info.eprime_label{end+1,1} = eprime_label;
    info.condition{end,1} = newlabel;
    info.onset_sec(end,1) = ons(k) / 1000;
    info.duration_sec(end,1) = dur(k) / 1000;
end

% Also capture counting blocks if that's the current label
countflag = false;

switch eprime_label
    case 'Fixation10'  % Run 1
        ons = eprime.Fixation10_OnsetTime(inds) + eprime.Fixation10_Duration(inds);
        dur = eprime.RESPONSE1_OnsetTime(inds) - ons;
        countflag = true;
    case 'Fixation11'  % Run 2
        ons = eprime.Fixation11_OnsetTime(inds) + eprime.Fixation11_Duration(inds);
        dur = eprime.RESPONSE3_OnsetTime(inds) - ons;
        countflag = true;
    case 'Fixation14'  % Run 3
        ons = eprime.Fixation14_OnsetTime(inds) + eprime.Fixation14_Duration(inds);
        dur = eprime.RESPONSE4_OnsetTime(inds) - ons;
        countflag = true;
    case 'Fixation16'  % Run 4
        ons = eprime.Fixation16_OnsetTime(inds) + eprime.Fixation16_Duration(inds);
        dur = eprime.RESPONSE7_OnsetTime(inds) - ons;
        countflag = true;
end

if countflag
    for k = 1:numel(inds)
        info.eprime_label{end+1,1} = 'NA';
        info.condition{end,1} = 'counting';
        info.onset_sec(end,1) = ons(k) / 1000;
        info.duration_sec(end,1) = dur(k) / 1000;
    end
end

