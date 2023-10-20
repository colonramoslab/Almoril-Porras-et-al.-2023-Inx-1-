function hm = toHashMap(sm,varargin) 
%function hm = toHashMap(sm,varargin) 

fields = fieldnames(sm);
fields = setdiff(fields, {'eset'});
for j = 1:length(sm)
    for k = 1:length(fields)
        no(j).(fields{k}) = sm(j).(fields{k});
        no(j).experimentNames = {sm(j).eset.expt.fname};
    end
end
no.allowedTransitions = reshape(no.allowedTransitions, [], 1);

hm = toHashMap(no);
