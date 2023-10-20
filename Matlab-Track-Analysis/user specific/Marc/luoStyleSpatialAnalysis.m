function analyzedData = luoStyleSpatialAnalysis (eset, spatial_navigation_options)
%function spatialNavigationMaggotFigures (esets, spatial_navigation_options)
%
% spatial_navigation_options -- 1 set of options for all esets


sno.angleBinSize = 30; % in degrees
sno.deltaAngleBinSize = 15;
sno.preferredDirection = 0;
sno.reoBinSize = 90;
sno.hsBinSize = 60;
sno.hsBinSpacing = 90;
sno.minHS = 0;
sno.minHSTheta = 20;
sno.validname = [];
sno.relativeDirField = [];
sno.dirOffsetField = [];
sno.timeBinSize = 10;
sno.validoperation = @(x) logical(setNonFiniteToZero(x));

if (nargin == 0)
    if (nargout >= 1)
        analyzedData = sno;
    end
    disp ('spatialNavigationMaggotAnalysis (esets, spatial_navigation_options)');
    return;
end
existsAndDefault('spatial_navigation_options', []);

if (length(eset) > 1)
    for j = 1:length(eset)
        disp (['start ' num2str(j)]); ts1 = tic;
        analyzedData(j) = spatialNavigationMaggotAnalysis(eset(j), spatial_navigation_options);  %#ok<AGROW>
        disp (['end ' num2str(j) ' ' num2str(toc(ts1)) ' s']);
    end
    return;
end

if (~isempty(spatial_navigation_options) && isstruct(spatial_navigation_options))
    fn = fieldnames(spatial_navigation_options);
    for j = 1:length(fn)
        sno.(fn{j}) = spatial_navigation_options.(fn{j});
    end
end

if (xor(isempty(sno.relativeDirField), isempty(sno.dirOffsetField)))
    disp('relativeDirField = theta - dirOffsetField');
    error('you must define both relativeDirField and dirOffsetField or neither');
    return;
end

tx = (-180:sno.angleBinSize:(180-sno.angleBinSize)) + sno.preferredDirection;
txc = (-180:sno.angleBinSize:(180)) + sno.preferredDirection;
txe = txc - sno.angleBinSize/2;

ad.tx = tx;
ad.txc = txc;
ad.txe = txe;

rs = eset.gatherSubField('run', 'startTheta');
re = eset.gatherSubField('run', 'endTheta');
rm = eset.gatherSubField('run', 'meanTheta');
rt = eset.gatherSubField('run', 'runTime');
if (~isempty(sno.validname))
    runvalid = sno.validoperation(eset.gatherFromSubField('run', sno.validname, 'position', 'start'));
    rs = rs(runvalid);
    re = re(runvalid);
    rm = rm(runvalid);
    rt = rt(runvalid);
end

if (~isempty(sno.dirOffsetField))
    roffset = eset.gatherFromSubField('run', sno.dirOffsetField, 'position', 'mean');
    roffset = roffset(runvalid);
else
    roffset = zeros(size(rs));
end

if(~isempty(sno.relativeDirField))
    fieldname = sno.relativeDirField;
else
    fieldname = 'theta';
end

rs = mod(rs - roffset + pi, 2*pi) - pi;
re = mod(re - roffset + pi, 2*pi) - pi;
rm = mod(rm - roffset + pi, 2*pi) - pi;


ad.runStartTheta = rs;
ad.runEndTheta = re;
ad.runMeanTheta = rm;
ad.runTime = rt;
ad.deltaRun = diff(unwrap([rs;re]));

rmrm = mod(rm - deg2rad(sno.preferredDirection) + pi, 2*pi) - pi;

towards = abs(rmrm) < pi/4;
away =  abs(rmrm) > 3*pi/4;
left =  abs(rmrm - pi/2) < pi/4;
right =  abs(rmrm + pi/2) < pi/4;

neutral = ~(towards | away);

ad.timeAxis = (0:sno.timeBinSize:(ceil(max(rt)/sno.timeBinSize)*sno.timeBinSize)) + sno.timeBinSize / 2;

h = hist(rt(towards), ad.timeAxis);
n = sum(h);
h = h/n;
heb = h.*(1-h)/sqrt(n);
ad.runTimeTowardsHist = h;
ad.runTimeTowardsHist_eb = heb;

h = hist(rt(away), ad.timeAxis);
n = sum(h);
h = h/n;
heb = h.*(1-h)/sqrt(n);
ad.runTimeAwayHist = h;
ad.runTimeAwayHist_eb = heb;

h = hist(rt(neutral), ad.timeAxis);
n = sum(h);
h = h/n;
heb = h.*(1-h)/sqrt(n);
ad.runTimeNeutralHist = h;
ad.runTimeNeutralHist_eb = heb;

[~,my,se] = eset.meanField2vsField1(fieldname, 'speed', deg2rad(txe), 'polar', true, 'runs','validname', sno.validname, 'validoperation', sno.validoperation);
ad.speedVsDir = my;
ad.speedVsDir_eb = se;


rsrs = mod(rs - deg2rad(sno.preferredDirection) + pi, 2*pi) - pi;

towards = abs(rsrs) < pi/4;
away =  abs(rsrs) > 3*pi/4;
left =  abs(rsrs - pi/2) < pi/4;
right =  abs(rsrs + pi/2) < pi/4;

dt = diff(unwrap([rs;re]));

ad.deltaThetaAxis = -180:sno.deltaAngleBinSize:180;
h = hist(dt(away), deg2rad(ad.deltaThetaAxis));
n = sum(h);
h = h/n;
heb = h.*(1-h)/sqrt(n);
ad.deltaThetaRunHistAway = h;
ad.deltaThetaRunHistAway_eb = heb;

ad.deltaThetaAxis = -180:sno.deltaAngleBinSize:180;
h = hist(dt(towards), deg2rad(ad.deltaThetaAxis));
n = sum(h);
h = h/n;
heb = h.*(1-h)/sqrt(n);
ad.deltaThetaRunHistTowards = h;
ad.deltaThetaRunHistTowards_eb = heb;

ad.deltaThetaAxis = -180:sno.deltaAngleBinSize:180;
h = hist(dt(left), deg2rad(ad.deltaThetaAxis));
n = sum(h);
h = h/n;
heb = h.*(1-h)/sqrt(n);
ad.deltaThetaRunHistLeft = h;
ad.deltaThetaRunHistLeft_eb = heb;

ad.deltaThetaAxis = -180:sno.deltaAngleBinSize:180;
h = hist(dt(right), deg2rad(ad.deltaThetaAxis));
n = sum(h);
h = h/n;
heb = h.*(1-h)/sqrt(n);
ad.deltaThetaRunHistRight = h;
ad.deltaThetaRunHistRight_eb = heb;

h = hist(adjustForPolarHistogram(rm, deg2rad(tx)), deg2rad(tx));
n = sum(h);
h = h/n;
heb = h.*(1-h)/sqrt(n);

ad.runDirectionThetaAxis = tx;
ad.runDirectionHistogram = h;
ad.runDirectionHistogram_eb = heb;


reotx = (-180:sno.reoBinSize:(180-sno.reoBinSize)) + sno.preferredDirection;
reotxc = (-180:sno.reoBinSize:(180)) + sno.preferredDirection;
reotxe = reotxc - sno.reoBinSize/2;
nhstemp = eset.gatherSubField('reorientation', 'numHS');
if (~isempty(sno.validname))
   reovalid = sno.validoperation(eset.gatherFromSubField('reorientation',  sno.validname,'position','start'));
else
   reovalid = true(size(nhstemp));
end
nhs = nhstemp(reovalid);
if (~isempty(sno.dirOffsetField))
    roffset = eset.gatherFromSubField('reorientation', sno.dirOffsetField, 'position', 'mean');
    roffset = roffset(reovalid);
else
    roffset = zeros(size(nhs));
end
ndir = eset.gatherSubField('reorientation', 'nextDir');
ndir= diff(unwrap([roffset;ndir(reovalid)]));
pdir = eset.gatherSubField('reorientation', 'prevDir');
pdir = diff(unwrap([roffset;pdir(reovalid)]));
reodir = diff(unwrap([pdir;ndir]));
rmag = abs(reodir);
%rmag = acos(cos(ndir - pdir));



nd = adjustForPolarHistogram(ndir(nhs >= sno.minHS), deg2rad(tx));  
pd = adjustForPolarHistogram(pdir(nhs >= sno.minHS), deg2rad(reotx));  
rm = rmag(nhs >= sno.minHS);
rd = reodir(nhs >= sno.minHS);

ad.reorientationPreviousDirection = pd;
ad.reorientationDeltaTheta = rd;

pdpd = mod(pd - deg2rad(sno.preferredDirection) + pi, 2*pi) - pi;

towards = abs(pdpd) < pi/4;
away =  abs(pdpd) > 3*pi/4;
left =  abs(pdpd - pi/2) < pi/4;
right =  abs(pdpd + pi/2) < pi/4;

h = hist(rd(right), deg2rad(ad.deltaThetaAxis));
n = sum(h);
h = h/n;
heb = h.*(1-h)/sqrt(n);
ad.deltaThetaReoHistRight = h;
ad.deltaThetaReoHistRight_eb = heb;

h = hist(rd(left), deg2rad(ad.deltaThetaAxis));
n = sum(h);
h = h/n;
heb = h.*(1-h)/sqrt(n);
ad.deltaThetaReoHistLeft = h;
ad.deltaThetaReoHistLeft_eb = heb;

h = hist(rd(towards), deg2rad(ad.deltaThetaAxis));
n = sum(h);
h = h/n;
heb = h.*(1-h)/sqrt(n);
ad.deltaThetaReoHistTowards = h;
ad.deltaThetaReoHistTowards_eb = heb;

h = hist(rd(away), deg2rad(ad.deltaThetaAxis));
n = sum(h);
h = h/n;
heb = h.*(1-h)/sqrt(n);
ad.deltaThetaReoHistAway = h;
ad.deltaThetaReoHistAway_eb = heb;

analyzedData = ad;
return;

[~,my,se] = meanyvsx(pd, rm,  deg2rad(reotxe));
reoMag = my;
reoStd = se;

[~,my,se] = meanyvsx(pd, rd,  deg2rad(reotxe));
reoDir = my;
reoDir_eb = se;

pd = adjustForPolarHistogram(pdir(nhs >= sno.minHS), deg2rad(reotx));  
[~,my,se] = meanyvsx(pd, nhs(nhs >= sno.minHS),  deg2rad(reotxe));
meanNumHS = my;
meanNumHS_eb = se;





hstx = (-180:sno.hsBinSpacing:(180-sno.hsBinSpacing)) + sno.preferredDirection;
hstxc = (-180:sno.hsBinSpacing:(180)) + sno.preferredDirection;
hstxe = hstxc - sno.hsBinSpacing/2;
hbw = sno.hsBinSize/2;

hstail = eset.gatherSubField('headSwing', 'tailDir');
if (~isempty(sno.dirOffsetField))
    hsoffset = eset.gatherFromSubField('headSwing', sno.dirOffsetField, 'position', 'mean');
else
    hsoffset = zeros(size(hstail));
end
hstail = diff(unwrap([hsoffset;hstail]));

hshead = eset.gatherSubField('headSwing', 'headDir');
hshead = diff(unwrap([hsoffset;hshead]));

hsacc = eset.gatherSubField('headSwing', 'accepted');
hssign = eset.gatherSubField('headSwing', 'sign');
htv = logical(eset.gatherSubField('headSwing', 'valid'));
hsmt = eset.gatherSubField('headSwing', 'maxTheta');
htv = htv & (rad2deg(abs(hsmt)) >= sno.minHSTheta);
if (~isempty(sno.validname))
   htv = htv & sno.validoperation(eset.gatherFromSubField('headSwing',  sno.validname, 'position','start'));
end
hsmt = hsmt(htv);
hstail = hstail(htv);
hshead = hshead(htv);
hsacc = hsacc(htv);
toleft = logical(hssign(htv) > 0);
toright = logical(hssign(htv) < 0);
tdeg = rad2deg(hstail);
tdeg = mod(tdeg, 360);
tdeg (tdeg > max(hstxe)) = tdeg(tdeg > max(hstxe)) - 360;
tdeg (tdeg < min(hstxe)) = tdeg(tdeg < min(hstxe)) + 360;

%[~,bin] = histc(tdeg, hstxe);
bin = -ones(size(tdeg));
for k = 1:length(hstx)
    bin(tdeg >= hstx(k)-hbw & tdeg <= hstx(k)+hbw) = k;
end
towards = cos(deg2rad(sno.preferredDirection) - hshead) > cos(deg2rad(sno.preferredDirection) - hstail);
[~,my,se] = meanyvsx(adjustForPolarHistogram(hstail,deg2rad(hstx)), double(toleft), deg2rad(hstxe));
    
allHSDir = my;
allHSDir_eb = se;
for k = 1:length(hstx)
    headSwingAcceptanceRateTowards(k) = mean(hsacc(bin == k & towards));
    headSwingAcceptanceRateTowards_eb(k) = headSwingAcceptanceRateTowards(k) * (1-headSwingAcceptanceRateTowards(k)) / sqrt(nnz(bin ==k & towards));

    headSwingAcceptanceRateAway(k) = mean(hsacc(bin == k & ~towards));
    headSwingAcceptanceRateAway_eb(k) = headSwingAcceptanceRateAway(k) * (1-headSwingAcceptanceRateAway(k)) / sqrt(nnz(bin ==k & ~towards));

    headSwingAcceptanceRateRight(k) = mean(hsacc(bin == k & toright));
    headSwingAcceptanceRateRight_eb(k) = headSwingAcceptanceRateRight(k) * (1-headSwingAcceptanceRateRight(k)) / sqrt(nnz(bin ==k & toright));

    headSwingAcceptanceRateLeft(k) = mean(hsacc(bin == k & toleft));
    headSwingAcceptanceRateLeft_eb(k) = headSwingAcceptanceRateLeft(k) * (1-headSwingAcceptanceRateLeft(k)) / sqrt(nnz(bin ==k & toleft));

    headSwingRejectionBiasRight(k) = nnz(~hsacc(bin == k & toright))/nnz(~hsacc(bin == k)) / mean(toright(bin == k));
    headSwingRejectionBiasRight_eb(k) = sqrt(nnz(~hsacc(bin == k & toright)))/nnz(~hsacc(bin == k)) / mean(toright(bin == k)); %this is a kludge for now
    
    headSwingRejectionBiasLeft(k) = nnz(~hsacc(bin == k & toleft))/nnz(~hsacc(bin == k)) / mean(toleft(bin == k));
    headSwingRejectionBiasLeft_eb(k) = sqrt(nnz(~hsacc(bin == k & toleft)))/nnz(~hsacc(bin == k)) / mean(toleft(bin == k)); %this is a kludge for now
    
    meanRejectedHeadSwingDir(k) = mean(hsmt(bin == k & ~hsacc));
    meanRejectedHeadSwingDir_eb(k) = std(hsmt(bin == k & ~hsacc))/sqrt(nnz(bin == k & ~hsacc));

    meanAcceptedHeadSwingDir(k) = mean(hsmt(bin == k & hsacc));
    meanAcceptedHeadSwingDir_eb(k) = std(hsmt(bin == k & hsacc))/sqrt(nnz(bin == k & hsacc));

   
end    




hstail = eset.gatherSubField('firsths', 'tailDir');
if (~isempty(sno.dirOffsetField))
    hsoffset = eset.gatherFromSubField('firsths', sno.dirOffsetField, 'position', 'mean');
else
    hsoffset = zeros(size(hstail));
end
hstail = diff(unwrap([hsoffset;hstail]));


hshead = eset.gatherSubField('firsths', 'headDir');
hshead = diff(unwrap([hsoffset;hshead]));

hssign = eset.gatherSubField('firsths', 'sign');
hsprevdir = eset.gatherSubField('firsths', 'prevDir');
hsprevdir = diff(unwrap([hsoffset;hsprevdir]));


htv = logical(eset.gatherSubField('firsths', 'valid'));
hsmt = double(eset.gatherSubField('firsths', 'maxTheta'));
htv = htv & (rad2deg(abs(hsmt)) >= sno.minHSTheta);
if (~isempty(sno.validname))
   htv = htv & sno.validoperation(eset.gatherFromSubField('firsths',  sno.validname, 'position','start'));
end
hsmt = hsmt(htv & hssign ~= 0);
hstail = hstail(htv& hssign ~= 0);
hshead = hshead(htv& hssign ~= 0);
hsprevdir = hsprevdir(htv & hssign ~= 0);
hssign = hssign(htv & hssign ~= 0);

towards = cos(deg2rad(sno.preferredDirection) - hshead) > cos(deg2rad(sno.preferredDirection) - hstail);
toleft = hssign > 0;
td = adjustForPolarHistogram(hstail, deg2rad(hstx));
pd = adjustForPolarHistogram(hsprevdir, deg2rad(hstx));
[~,my,se] = meanyvsx(pd, towards, deg2rad(hstxe));
firstHSBias = my;
firstHSBias_eb =se;
[~,my,se] = meanyvsx(td, toleft, deg2rad(hstxe));
firstHSDir = my;
firstHSDir_eb = se;

[~,my,se] = meanyvsx(td, hsmt, deg2rad(hstxe));
firstHSMeanDir = my;
firstHSMeanDir_eb = se;

    


fields = {'tx', 'txc', 'txe', 'hstx', 'hstxc','hstxe', 'thetahist', 'reohist', 'runStartDirectionHist', 'headSwingAcceptanceRateTowards',...
    'headSwingAcceptanceRateAway', 'reoMag', 'reoStd', 'meanNumHS', 'speedVsDir','firstHSBias', 'firstHSDir',...
    'headSwingAcceptanceRateLeft', 'headSwingAcceptanceRateRight','meanrunchange','instantaneousdthetavstheta',...
    'thetahist_eb', 'reohist_eb', 'runStartDirectionHist_eb', 'headSwingAcceptanceRateTowards_eb',...
    'headSwingAcceptanceRateAway_eb', 'speedVsDir_eb','firstHSBias_eb', 'firstHSDir_eb',...
    'headSwingAcceptanceRateLeft_eb', 'headSwingAcceptanceRateRight_eb','meanrunchange_eb','instantaneousdthetavstheta_eb',...
    'allHSDir', 'allHSDir_eb', 'reoDir', 'reoDir_eb', 'reotx', 'reotxe', 'reotxc','headSwingRejectionBiasRight','headSwingRejectionBiasRight_eb','headSwingRejectionBiasLeft','headSwingRejectionBiasLeft_eb',...
    'meanRejectedHeadSwingDir', 'meanRejectedHeadSwingDir_eb','meanAcceptedHeadSwingDir','meanAcceptedHeadSwingDir_eb', 'firstHSMeanDir', 'firstHSMeanDir_eb'};
for j = 1:length(fields)
    try
        eval(['analyzedData.' fields{j} ' = ' fields{j} ';']);
    catch me
        disp(me.getReport);
    end
end



