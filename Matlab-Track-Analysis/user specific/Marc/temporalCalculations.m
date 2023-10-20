function tc = temporalCalculations (eset, tno)
%
%function tc = temporalCalculations (eset, tno)
%
%tno must have fieldname and ramptype defined
if (length(eset) > 1)
    for j = 1:length(eset)
        disp (['start ' num2str(j)]);
        ts = tic;
        tc(j) = temporalCalculations(eset(j), tno); %#ok<AGROW>
        disp (['done with ' num2str(j) ' - ' num2str(toc(ts))]);
    end
    return
end
ts1 = tic;
tc.hs_htv = [];
tc.hs_htv = logical(eset.gatherSubField('headSwing', 'valid'));
toc(ts1);
tc.reo_nhs = eset.gatherSubField('reorientation', 'numHS');
toc(ts1);
tc.hs_ton = eset.gatherFromSubField('headSwing', [tno.fieldname '_ton'], 'position', 'start');
tc.hs_toff = eset.gatherFromSubField('headSwing', [tno.fieldname '_toff'], 'position', 'start');
tc.hs_time = eset.gatherFromSubField('headSwing', 'eti', 'position', 'start');
toc(ts1);
tc.reo_ton = eset.gatherFromSubField('reorientation', [tno.fieldname '_ton'], 'position', 'start');
tc.reo_toff = eset.gatherFromSubField('reorientation', [tno.fieldname '_toff'], 'position', 'start');
tc.reo_time = eset.gatherFromSubField('reorientation', 'eti', 'position','start');
toc(ts1);
if (strcmpi(tno.rampType, 'triangle') || strcmpi(tno.rampType, 'exponential'))
    tc.hs_rising = logical(eset.gatherFromSubField('headSwing', [tno.fieldname '_rising'], 'position', 'start'));
    tc.hs_falling = logical(eset.gatherFromSubField('headSwing', [tno.fieldname '_falling'], 'position', 'start'));
    tc.reo_rising = logical(eset.gatherFromSubField('reorientation', [tno.fieldname '_rising'], 'position', 'start'));
    tc.reo_falling = logical(eset.gatherFromSubField('reorientation', [tno.fieldname '_falling'], 'position', 'start'));
end
toc(ts1);
if (strcmpi(tno.rampType, 'square'))
    tc.hs_high = eset.gatherFromSubField('headSwing', [tno.fieldname '_high'], 'position', 'start');
    tc.hs_low = eset.gatherFromSubField('headSwing', [tno.fieldname '_low'], 'position', 'start');
    tc.reo_high = eset.gatherFromSubField('reorientation', [tno.fieldname '_high'], 'position', 'start');
    tc.reo_low = eset.gatherFromSubField('reorientation', [tno.fieldname '_low'], 'position', 'start');
end
toc(ts1);
tc.reo_mag = diff(unwrap([eset.gatherSubField('reorientation', 'nextDir');eset.gatherSubField('reorientation', 'prevDir')]));

tc.hs_acc = logical(eset.gatherSubField('headSwing', 'accepted'));
tc.hs_mag = eset.gatherSubField('headSwing','maxTheta');
toc(ts1);
disp('done')