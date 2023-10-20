function addTonToff(expt, fieldname, ramptype, varargin)
%addTonToff(expt, fieldname, ramptype, varargin);

if (length(expt) > 1)
    for j = 1:length(expt)
        expt(j).addTonToff(fieldname, ramptype, varargin{:});
    end
    return
end

ind = find(strcmpi(fieldname, {expt.globalQuantity.fieldname}));
if (isempty(ind))
    disp (['could not find global quantity named ' fieldname]);
    return;
end

if (length(ind) > 1)
    ind = find(strcmp(fieldname, {expt.globalQuantity.fieldname}));
end

if (isempty(ind))
    disp (['found more than one case-insensitive match, but no exact matches for fieldname ' fieldname]);
    return;
end

if (length(ind) > 1)
    disp (['found two global quantities with same name: ', fieldname]);
    return;
end

gqs = expt.globalQuantity(ind).timeOnOffGQs(ramptype, varargin);

for j = 1:length(gqs)
    expt.addGlobalQuantity(gqs(j));
end
