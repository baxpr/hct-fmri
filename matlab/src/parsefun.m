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
if strcmp(eprime_label,'Fixation10')
    ons = eprime.Fixation10_OnsetTime(inds) + eprime.Fixation10_Duration(inds);
    dur = eprime.RESPONSE1_OnsetTime(inds) - ons;
    for k = 1:numel(inds)
        info.eprime_label{end+1,1} = 'NA';
        info.condition{end,1} = 'counting';
        info.onset_sec(end,1) = ons(k) / 1000;
        info.duration_sec(end,1) = dur(k) / 1000;
    end
end
