function ad = temporalNavigationAnalysis (esets, temporalFieldName, temporalNavigationOptions, temporalCalcs)
%function ad = temporalNaviationAnalysis (esets, temporalFieldName, temporalNavigationOptions, temporalCalcs)

if (nargin >= 1 && length(esets) > 1)
    existsAndDefault('temporalNavigationOptions', []);
    if (~isfield(temporalNavigationOptions, 'period'))
        maxT = 0;
        for j = 1:length(esets)
            maxT = max(maxT, max(esets(j).gatherField([temporalFieldName '_ton'])));
        end
        temporalNavigationOptions.period = maxT;
    end
    for j = 1:length(esets)
        if (exist('temporalCalcs', 'var') && length(temporalCalcs) >= j)
            ad(j) = temporalNavigationAnalysis (esets(j), temporalFieldName, temporalNavigationOptions, temporalCalcs(j)); %#ok<AGROW>
        else
            ad(j) = temporalNavigationAnalysis(esets(j), temporalFieldName, temporalNavigationOptions); %#ok<AGROW>
        end
    end
    return;
end

tno.timeBinSize = 30; % in seconds
tno.rampType = 'triangle';
tno.minHS = 1;
tno.timerange = [];
tno.fieldname = [];
tno.period = [];
if (nargin == 0)
    ad = tno;
    return;
end

tno.fieldname = temporalFieldName;
ad.fieldname = temporalFieldName;
ad.rampType = tno.rampType;
eset = esets(1);
tno.period = max(eset.gatherField([temporalFieldName '_ton']));


if (nargin > 2 && isstruct(temporalNavigationOptions))
    fn = fieldnames(temporalNavigationOptions);
    for j = 1:length(fn)
        tno.(fn{j}) = temporalNavigationOptions.(fn{j});
    end
end


if (existsAndDefault('temporalCalcs', []))
    tc = temporalCalcs;
else
    tc = temporalCalculations('eset', tno);
end
ad.tno = tno;
ad.tc = tc;
ad.period = tno.period;
etxe = 0:tno.timeBinSize:tno.period;
etx = etxe(1:(end-1)) + tno.timeBinSize/2;
ad.etxf = 0:1:tno.period;
[~,ind] = find(etx > tno.period / 2);
shiftinds = [ind:length(etx), 1:(ind-1)];  
ad.etxs = mod(etx(shiftinds) + tno.period/2, tno.period) - tno.period/2;


if (~isempty(tno.timerange))
    reovalid = tc.reo_time >= tno.timerange(1) & tc.reo_time <= tno.timerange(2);
    hsvalid = tc.hs_htv & tc.hs_time >= tno.timerange(1) & tc.hs_time <= tno.timerange(2);
else
    reovalid = true(size(tc.reo_time));
    hsvalid = tc.hs_htv;
end


[~,nhs, nhs_eb] = meanyvsx(tc.reo_ton(tc.reo_nhs >= tno.minHS & reovalid), tc.reo_nhs(tc.reo_nhs >= tno.minHS & reovalid), etxe);
ad.nhs_vs_ton = nhs(shiftinds);
ad.nhs_vs_ton_eb = nhs_eb(shiftinds);
[~,nhs, nhs_eb] = meanyvsx(tc.reo_toff(tc.reo_nhs >= tno.minHS & reovalid), tc.reo_nhs(tc.reo_nhs >= tno.minHS & reovalid), etxe);
ad.nhs_vs_toff = nhs(shiftinds);
ad.nhs_vs_toff_eb = nhs_eb(shiftinds);

[~,sp, sp_eb] = eset.meanField2vsField1([temporalFieldName '_ton'], 'speed', etxe, 'runs', 'timerange', tno.timerange);
ad.sp_vs_ton = sp(shiftinds);
ad.sp_vs_ton_eb = sp_eb(shiftinds);
[~,sp, sp_eb] = eset.meanField2vsField1([temporalFieldName '_toff'], 'speed', etxe, 'runs', 'timerange', tno.timerange);
ad.sp_vs_toff = sp(shiftinds);
ad.sp_vs_toff_eb = sp_eb(shiftinds);

[~,acc, acc_eb] = meanyvsx(tc.hs_ton(hsvalid), tc.hs_acc(hsvalid), etxe);
ad.hs_acc_vs_ton = acc(shiftinds);
ad.hs_acc_vs_ton_eb = acc_eb(shiftinds);

[~,acc, acc_eb] = meanyvsx(tc.hs_toff(hsvalid), tc.hs_acc(hsvalid), etxe);
ad.hs_acc_vs_toff = acc(shiftinds);
ad.hs_acc_vs_toff_eb = acc_eb(shiftinds);


if (strcmpi(tno.rampType, 'triangle') || strcmpi(tno.rampType, 'exponential'))
    ad.nhs_hist_falling = histc(tc.reo_nhs(logical(tc.reo_falling) & tc.reo_nhs >= tno.minHS & reovalid), 0:1:10)/nnz(tc.reo_falling & tc.reo_nhs >= tno.minHS & reovalid);
    ad.nhs_hist_rising = histc(tc.reo_nhs(logical(tc.reo_rising) & tc.reo_nhs >= tno.minHS& reovalid), 0:1:10)/nnz(tc.reo_rising & tc.reo_nhs >= tno.minHS& reovalid);
    ad.hs_acc_falling = mean(tc.hs_acc(tc.hs_falling));% & tc.hs_htv));
    ad.hs_acc_rising = mean(tc.hs_acc(tc.hs_rising));% & tc.hs_htv));
else
    ad.nhs_hist_falling = [];
    ad.nhs_hist_rising = [];
    ad.hs_acc_falling = [];
    ad.hs_acc_rising = [];
end

 


[rr, rr_eb] = eset.makeReorientationHistogram([temporalFieldName '_ton'], etx, 'minHS', tno.minHS, 'incllastrun', true, 'timerange', tno.timerange);
ad.reo_vs_ton = rr(shiftinds);
ad.reo_vs_ton_eb = rr_eb(shiftinds);
[rr, rr_eb] = eset.makeReorientationHistogram([temporalFieldName '_toff'], etx, 'minHS', tno.minHS, 'incllastrun', true, 'timerange', tno.timerange);
ad.reo_vs_toff = rr(shiftinds);
ad.reo_vs_toff_eb = rr_eb(shiftinds);

[~,val] = eset.meanField2vsField1([temporalFieldName '_ton'], temporalFieldName, etxe, 'timerange', tno.timerange);
ad.val_vs_ton = val(shiftinds);
[~,val] = eset.meanField2vsField1([temporalFieldName '_toff'], temporalFieldName, etxe, 'timerange', tno.timerange);
ad.val_vs_toff = val(shiftinds);

[~,val] = eset.meanField2vsField1([temporalFieldName '_ton'], temporalFieldName, 0:1:tno.period, 'timerange', tno.timerange);
ad.val_vs_ton_highres = val;
[~,val] = eset.meanField2vsField1([temporalFieldName '_toff'], temporalFieldName, 0:1:tno.period, 'timerange', tno.timerange);
ad.val_vs_toff_highres = val;


[~,rm,rm_eb] = meanyvsx(tc.reo_ton(tc.reo_nhs >= tno.minHS & reovalid), abs(tc.reo_mag(tc.reo_nhs >= tno.minHS & reovalid)), etxe);
ad.reo_mag_vs_ton = rm(shiftinds);
ad.reo_mag_vs_ton_eb = rm_eb(shiftinds);

[~,rm,rm_eb] = meanyvsx(tc.reo_toff(tc.reo_nhs >= tno.minHS & reovalid), abs(tc.reo_mag(tc.reo_nhs >= tno.minHS & reovalid)), etxe);
ad.reo_mag_vs_toff = rm(shiftinds);
ad.reo_mag_vs_toff_eb = rm_eb(shiftinds);


[~,rm,rm_eb] = meanyvsx(tc.reo_ton(tc.reo_nhs >= tno.minHS & reovalid), tc.reo_mag(tc.reo_nhs >= tno.minHS & reovalid), etxe);
ad.reo_dir_vs_ton = rm(shiftinds);
ad.reo_dir_vs_ton_eb = rm_eb(shiftinds);

[~,rm,rm_eb] = meanyvsx(tc.reo_toff(tc.reo_nhs >= tno.minHS & reovalid), tc.reo_mag(tc.reo_nhs >= tno.minHS & reovalid), etxe);
ad.reo_dir_vs_toff = rm(shiftinds);
ad.reo_dir_vs_toff_eb = rm_eb(shiftinds);


