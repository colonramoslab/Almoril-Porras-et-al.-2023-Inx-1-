function analyzedData = spatialNavigationMaggotAnalysis (eset, spatial_navigation_options)
%function spatialNavigationMaggotFigures (esets, spatial_navigation_options)
%
% spatial_navigation_options -- 1 set of options for all esets
% plot_options -- each eset gets its own options
% plot_options.
%    lineWidth -- width of lines
%    color -- color for line and marker
%    legendEntry -- what to put in legend
%    marker -- marker 
%    plotOptions -- additional options to pass to plotting function

sno.angleBinSize = 30; % in degrees
sno.dtAngleBinSize = 20; % in degrees
sno.hsdtAngleBinSize = 20; % in degrees
sno.preferredDirection = 0;
sno.reoBinSize = 30;
sno.hsBinSize = 90;
sno.hsBinSpacing = 90;
sno.minHS = 0;
sno.minHSTheta = 20;
sno.validname = [];
sno.relativeDirField = [];
sno.dirOffsetField = [];
sno.validoperation = @(x) logical(setNonFiniteToZero(x));
sno.confidenceLevel = 0.95;
sno.autocorr_timerange = [];
sno.runTimeBinSize = 10;

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


dtx = (-180:sno.dtAngleBinSize:(180-sno.dtAngleBinSize));
dtxc = (-180:sno.dtAngleBinSize:180);

hsdtx = sort([(sno.minHSTheta - sno.hsdtAngleBinSize/2):(-sno.hsdtAngleBinSize):0 (sno.minHSTheta + sno.hsdtAngleBinSize/2):(sno.hsdtAngleBinSize):180]);
hsdtx = sort([-hsdtx hsdtx]);
%hsdtx = hsdtxc(1:(end-1));

rs = eset.gatherSubField('run', 'startTheta');
re = eset.gatherSubField('run', 'endTheta');
rt = eset.gatherSubField('run', 'runTime');
rmt = eset.gatherSubField('run', 'meanTheta');
if (~isempty(sno.validname))
    runvalid = sno.validoperation(eset.gatherFromSubField('run', sno.validname, 'position', 'start'));
    rs = rs(runvalid);
    re = re(runvalid);
    rt = rt(runvalid);
    rmt = rmt(runvalid);
end
dt = diff(unwrap([rs;re]));
rs = adjustForPolarHistogram(rs, deg2rad(tx));
changeLessThan90 = logical(abs(dt) < pi/2);
[~,my, se] = meanyvsx(rs(changeLessThan90), dt(changeLessThan90), deg2rad(txe));
meanrunchange = my;
meanrunchange_eb = se;

analyzedData.runStartTheta = rs;
analyzedData.runEndTheta = re;
analyzedData.runMeanTheta = rmt;
analyzedData.runTime = rt;
minruntime = median(eset.gatherSubField('so', 'minRunTime'));
analyzedData.runTimeAxis = minruntime + (sno.runTimeBinSize/2)+ (0:sno.runTimeBinSize:(sno.runTimeBinSize*ceil(600/sno.runTimeBinSize)));
rth = hist(rt(cos(rmt - deg2rad(sno.preferredDirection)) > 1/sqrt(2)), analyzedData.runTimeAxis);
analyzedData.runTimeHistTowards = (rth/sum(rth));
analyzedData.runTimeHistTowards_eb = sqrt((rth/sum(rth)).*(1-(rth/sum(rth)))./sum(rth));
rth = hist(rt(-cos(rmt - deg2rad(sno.preferredDirection)) > 1/sqrt(2)), analyzedData.runTimeAxis);
analyzedData.runTimeHistAway = (rth/sum(rth));
analyzedData.runTimeHistAway_eb = sqrt((rth/sum(rth)).*(1-(rth/sum(rth)))./sum(rth));


if(~isempty(sno.relativeDirField))
    fieldname = sno.relativeDirField;
else
    fieldname = 'theta';
end
[~,my,se]= eset.meanField2vsField1(fieldname, 'lrdtheta', deg2rad(txe), 'runs', 'polar', true,'validname', sno.validname, 'validoperation', sno.validoperation);
instantaneousdthetavstheta = my;
instantaneousdthetavstheta_eb = se;

v = eset.gatherField('vel', 'validname', sno.validname, 'validoperation', sno.validoperation);
M = [cosd(sno.preferredDirection) sind(sno.preferredDirection); -sind(sno.preferredDirection) cosd(sno.preferredDirection)];
v = M*v;
s = sqrt(sum(v.^2));
navind = mean(v,2)/mean(s);
navind_eb = std(v,0,2)/mean(s)/sqrt(length(s)*eset.expt(1).track(1).dr.interpTime/(2*eset.autocorr_tau));
for k = 1:length(eset.expt)
    v = eset.expt(k).gatherField('vel', 'validname', sno.validname, 'validoperation', sno.validoperation);
    M = [cosd(sno.preferredDirection) sind(sno.preferredDirection); -sind(sno.preferredDirection) cosd(sno.preferredDirection)];
    v = M*v;
    s = sqrt(sum(v.^2));
    navind_expt(:,k) = mean(v,2)/mean(s);
    navind_expt_eb(:,k) = std(v,0,2)/mean(s)/sqrt(length(s)*eset.expt(1).track(1).dr.interpTime/(2*eset.autocorr_tau));
end

[h,eb] = eset.makeHistogram(fieldname, deg2rad(tx), 'runs', 'polar', true,'validname', sno.validname, 'validoperation', sno.validoperation);
thetahist_eb = eb ./ mean(h);
thetahist = h./mean(h);



[h,eb] = eset.makeReorientationHistogram(fieldname, deg2rad(tx), 'minHS',sno.minHS, 'polar', true,'validname', sno.validname, 'validoperation', sno.validoperation);
reohist = h;
reohist_eb = eb;



%{
temptx = (-180:90:135) + sno.preferredDirection;
[h,eb] = eset.makeReorientationHistogram(fieldname, deg2rad(temptx), 'minHS',sno.minHS, 'polar', true,'validname', sno.validname, 'validoperation', sno.validoperation);
reoRateIndex = (h(1)-h(3))/(h(1) + h(3));
reoRateIndex_eb = sqrt((eb(1).^2 + eb(3).^2)/2)/(h(1)+h(3));
%}


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

allth = eset.gatherField(fieldname, 'runs','validname',sno.validname, 'validoperation', sno.validoperation);

reo_rate_model = fitReorientationRateTheta(eset, allth, pdir);
reoRateIndex = reo_rate_model.params(2)/reo_rate_model.params(1) * cos(reo_rate_model.params(3) - deg2rad(sno.preferredDirection));
[minx maxx] = confidenceRange(reo_rate_model.params, reo_rate_model.hess, @(x) x(2,:)./x(1,:).*cos(x(3) - deg2rad(sno.preferredDirection)), sno.confidenceLevel);
reoRateIndex_ci = [minx maxx];
reoRateIndex_eb =  stdOfFunVal(reo_rate_model.params, reo_rate_model.hess, @(x) x(2,:)./x(1,:).*cos(x(3) - deg2rad(sno.preferredDirection)));

reoRateIndexPerp = reo_rate_model.params(2)/reo_rate_model.params(1) * sin(reo_rate_model.params(3) - deg2rad(sno.preferredDirection));
[minx maxx] = confidenceRange(reo_rate_model.params, reo_rate_model.hess, @(x) x(2,:)./x(1,:).*sin(x(3) - deg2rad(sno.preferredDirection)), sno.confidenceLevel);
reoRateIndexPerp_eb = 0.5*(maxx-minx);

reo_prevdir = pdir;%(nhs >= sno.minHS);
reo_dtheta = reodir;%(nhs >= sno.minHS);
reo_numHS = nhs;
reo_dir_model = fitReorientationAngleDistributionWithHSInfo(pdir, ndir, nhs);

analyzedData.nixedParams = {{}, {'C'}, {'B', 'E'}, {'D', 'E'}, {'B','C','E'}};
analyzedData.altmodeldescriptions = {'orig model', 'no left/right bias', 'no size bias', 'no skew', 'no bias'};
for j = 1:length(analyzedData.nixedParams)
    clear startValues fixedValues
    for k = 1:length(reo_dir_model.params)
        startValues.(reo_dir_model.paramkey{k}) = reo_dir_model.params(k);
    end
    fixedValues = [];
    for k = 1:length(analyzedData.nixedParams{j})
        startValues.(analyzedData.nixedParams{j}{k}) = 0;
        fixedValues.(analyzedData.nixedParams{j}{k}) = 0;
    end
    m1 = fitReorientationAngleDistribution(pdir(nhs > 0), diff(unwrap([pdir(nhs>0);ndir(nhs>0)])), fixedValues); %nix start values
    m2 = fitReorientationAngleDistribution(pdir(nhs > 0), diff(unwrap([pdir(nhs>0);ndir(nhs>0)])), fixedValues, startValues); %with start values
    if (m1.logLikelihood > m2.logLikelihood)
        analyzedData.altmodel(j) = m1;
    else
        analyzedData.altmodel(j) = m2;
    end
end
analyzedData.deltaLL = [analyzedData.altmodel.logLikelihood] - reo_dir_model.logLikelihood;
% 'A'    'B'    'C'    'U'    'sigma1'    'sigma2'    'skew'    'theta_0'
% P(dt | theta) = a*N(0,sigma1,dt) + (1/2 - a/2 + g) * S(U, sigma2, skew, dt) +
%    (1/2 - a/2 - g) * S(U, sigma2, skew, -dt)
%
%a(theta) = A + B*cos(theta - theta_0); (|B| < A);
%g(theta) = C(sin(theta - theta_0)) 
%other fit parameters, U, sigma1, sigma2, skew, theta_0 
dirIndexFun = @(x) 2*x(3).*cos(deg2rad(sno.preferredDirection) - x(7));
reoDirIndex = dirIndexFun(reo_dir_model.params);
[minx maxx] = confidenceRange(reo_dir_model.params, reo_dir_model.hessian, dirIndexFun, sno.confidenceLevel);
reoDirIndex_ci = [minx maxx];
reoDirIndex_eb = stdOfFunVal(reo_dir_model.params, reo_dir_model.hessian, dirIndexFun);

% P(dt | theta) = (0.5 - g) * S(U, sigma, skew, dt) + (0.5 + g) * S(U,
% sigma, skew, -dt)
%
% U = A - B*cos(theta - theta_0);
% skew = D - E*cos(theta - theta_0);
% g = C*(sin(theta-theta_0))
%other fit parameters, sigma, theta_0 
%
%model.params = [A, B, C, D, E, sigma, theta_0];
magIndexFun = @(x) -(skewNormalMeanSq(x(1) - x(2)*cos((deg2rad(sno.preferredDirection) - x(7))), x(6), x(4) - x(5)*cos((deg2rad(sno.preferredDirection) - x(7)))) - ...
    skewNormalMeanSq(x(1) + x(2)*cos((deg2rad(sno.preferredDirection) - x(7))), x(6), x(4) + x(5)*cos((deg2rad(sno.preferredDirection) - x(7))))) ./ ...
    skewNormalMeanSq(x(1), x(6), x(4));
reoMagIndex = magIndexFun(reo_dir_model.params);
[minx maxx] = confidenceRange(reo_dir_model.params, reo_dir_model.hessian, magIndexFun, sno.confidenceLevel);
reoMagIndex_ci = [minx maxx];
reoMagIndex_eb = stdOfFunVal(reo_dir_model.params, reo_dir_model.hessian, magIndexFun);

%{

reoDirIndex = 2*reo_dir_model.params(1)*cos(reo_dir_model.params(5)-deg2rad(sno.preferredDirection)).*reo_dir_model.params(2)./sqrt(reo_dir_model.params(2).^2 + reo_dir_model.params(4).^2);
[minx maxx] = confidenceRange(reo_dir_model.params, reo_dir_model.hessian, @(x) 2*x(1,:).*cos(x(5,:)-deg2rad(sno.preferredDirection)).*x(2,:)./sqrt(x(2,:).^2 + x(4,:).^2), sno.confidenceLevel);
reoDirIndex_eb = 0.5*(maxx-minx);
reoMagIndex = reo_dir_model.params(3)./sqrt(reo_dir_model.params(2).^2 + reo_dir_model.params(4).^2)*cos(reo_dir_model.params(5)-deg2rad(sno.preferredDirection));
[minx maxx] = confidenceRange(reo_dir_model.params, reo_dir_model.hessian, @(x) x(3,:)./sqrt(x(2,:).^2 + x(4,:).^2).*cos(x(5,:)-deg2rad(sno.preferredDirection)), sno.confidenceLevel);
reoMagIndex_eb = 0.5*(maxx-minx);
%}

nd = adjustForPolarHistogram(ndir(nhs >= sno.minHS), deg2rad(tx));  
runStartDirectionHist = hist(nd, deg2rad(tx));
h = runStartDirectionHist;
n = sum(h);

runStartDirectionHist_eb = sqrt(h/n).*sqrt(1-h/n).*sqrt(n) ./ mean(runStartDirectionHist);
runStartDirectionHist = runStartDirectionHist./mean(runStartDirectionHist);
 
topref = cos(reo_prevdir - deg2rad(sno.preferredDirection)) > 1/sqrt(2);
frompref = cos(reo_prevdir - deg2rad(sno.preferredDirection)) < -1/sqrt(2);
leftofpref = sin(reo_prevdir - deg2rad(sno.preferredDirection)) > 1/sqrt(2);
rightofpref = sin(reo_prevdir - deg2rad(sno.preferredDirection)) < -1/sqrt(2);

indset = {topref, frompref, leftofpref, rightofpref};


dtxf = deg2rad(-180:1:180);
ddtxf = deg2rad(diff(dtxf(1:2)));
dt = adjustForPolarHistogram(reo_dtheta, deg2rad(dtx)); 
for j = 1:length(indset)
    reo_dtheta_dist_hs{j} = hist(dt(indset{j} & nhs > 0), deg2rad(dtx));
    reo_dtheta_dist_hs_eb{j} = sqrt(reo_dtheta_dist_hs{j}.*(1-reo_dtheta_dist_hs{j}/nnz(indset{j} & nhs > 0)));
    reo_dtheta_dist_nohs{j} = hist(dt(indset{j} & nhs == 0), deg2rad(dtx));
    reo_dtheta_dist_nohs_eb{j} = sqrt(reo_dtheta_dist_nohs{j}.*(1-reo_dtheta_dist_nohs{j}/nnz(indset{j} & nhs == 0)));
end
for j = 1:length(indset)
%    reo_dtheta_fit_m{j} = zeros(size(dtxf));
    reo_dtheta_fit_r{j} = zeros(size(dtxf));
    reo_dtheta_fit_l{j} = zeros(size(dtxf));
    for k = find(indset{j} & nhs >= sno.minHS)
        [r,l] = reo_dir_model.pdfOfdThetaBreakdown(reo_dir_model.params, reo_prevdir(k), dtxf);
        %reo_dtheta_fit_m{j} = reo_dtheta_fit_m{j} + m*deg2rad(sno.dtAngleBinSize);
        reo_dtheta_fit_r{j} = reo_dtheta_fit_r{j} + r*deg2rad(sno.dtAngleBinSize);
        reo_dtheta_fit_l{j} = reo_dtheta_fit_l{j} + l*deg2rad(sno.dtAngleBinSize);
    end
    reo_dtheta_fit{j} = reo_dtheta_fit_r{j} + reo_dtheta_fit_l{j};
end
%{
reo_dtheta_dist_topref = hist(dt(topref), deg2rad(dtx));
reo_dtheta_dist_topref_eb = sqrt(reo_dtheta_dist_topref.*(1-reo_dtheta_dist_topref/nnz(topref)));
reo_dtheta_dist_frompref = hist(dt(frompref), deg2rad(dtx));
reo_dtheta_dist_frompref_eb = sqrt(reo_dtheta_dist_frompref.*(1-reo_dtheta_dist_frompref/nnz(topref)));
reo_dtheta_dist_leftofpref = hist(dt(leftofpref), deg2rad(dtx));
reo_dtheta_dist_leftofpref_eb = sqrt(reo_dtheta_dist_leftofpref.*(1-reo_dtheta_dist_leftofpref/nnz(topref)));
reo_dtheta_dist_rightofpref = hist(dt(rightofpref), deg2rad(dtx));
reo_dtheta_dist_rightofpref_eb = sqrt(reo_dtheta_dist_rightofpref.*(1-reo_dtheta_dist_rightofpref/nnz(topref)));

reo_dtheta_fit_topref = zeros(size(dtxf));
for j = find(topref)
    reo_dtheta_fit_topref = reo_dtheta_fit_topref + reo_dir_model.pdfOfdTheta(reo_dir_model.params, reo_prevdir(j), dtxf)*deg2rad(sno.dtAngleBinSize);
end
reo_dtheta_fit_frompref = zeros(size(dtxf));
for j = find(frompref)
    reo_dtheta_fit_frompref = reo_dtheta_fit_frompref + reo_dir_model.pdfOfdTheta(reo_dir_model.params, reo_prevdir(j), dtxf)*deg2rad(sno.dtAngleBinSize);
end
reo_dtheta_fit_leftofpref = zeros(size(dtxf));
for j = find(leftofpref)
    reo_dtheta_fit_leftofpref = reo_dtheta_fit_leftofpref + reo_dir_model.pdfOfdTheta(reo_dir_model.params, reo_prevdir(j), dtxf)*deg2rad(sno.dtAngleBinSize);
end
reo_dtheta_fit_rightofpref = zeros(size(dtxf));
for j = find(rightofpref)
    reo_dtheta_fit_rightofpref = reo_dtheta_fit_rightofpref + reo_dir_model.pdfOfdTheta(reo_dir_model.params, reo_prevdir(j), dtxf)*deg2rad(sno.dtAngleBinSize);
end
%}

dtxf = rad2deg(dtxf);
reobasedirections = [mod(sno.preferredDirection + 180, 360) - 180, mod(sno.preferredDirection, 360) - 180, mod(sno.preferredDirection + 270, 360) - 180, mod(sno.preferredDirection + 90, 360) - 180];
%{
reo_dtheta_dist = {reo_dtheta_dist_topref, reo_dtheta_dist_frompref, reo_dtheta_dist_leftofpref, reo_dtheta_dist_rightofpref};
reo_dtheta_dist_eb = {reo_dtheta_dist_topref_eb, reo_dtheta_dist_frompref_eb, reo_dtheta_dist_leftofpref_eb, reo_dtheta_dist_rightofpref_eb};
reo_dtheta_fit = {reo_dtheta_fit_topref, reo_dtheta_fit_frompref, reo_dtheta_fit_leftofpref, reo_dtheta_fit_rightofpref};
%}
rm = rmag(nhs >= sno.minHS);
meanReoMag = mean(rm);
frompref = frompref(nhs >= sno.minHS);
topref = topref(nhs >= sno.minHS);
reoMagIndex_old = (mean(rm(frompref)) - mean(rm (topref)))/meanReoMag;
reoMagIndex_old_eb = sqrt ((var(rm(topref))*nnz(topref) + var(rm (frompref))*nnz(frompref)))/(nnz(topref)+nnz(frompref))/meanReoMag;


pd = adjustForPolarHistogram(reo_prevdir, deg2rad(reotx)); 
[~,my,se] = meanyvsx(pd(nhs >= sno.minHS), reodir(nhs >= sno.minHS).^2,  deg2rad(reotxe));
reoMag = my;
reoStd = se;
reoMag_eb = se; 




[~,my,se] = meanyvsx(pd(nhs >= sno.minHS), reodir(nhs >= sno.minHS),  deg2rad(reotxe));
reoDir = my;
reoDir_eb = se;
rd = reodir(nhs >= sno.minHS);
rightofpref = rightofpref(nhs >= sno.minHS);
leftofpref = leftofpref(nhs >= sno.minHS);
reoDirIndex_old = (mean(rd(rightofpref)) - mean(rd (leftofpref)))/meanReoMag;
reoDirIndex_old_eb = sqrt ((var(rd(rightofpref))*nnz(rightofpref) + var(rd (leftofpref))*nnz(leftofpref)))/(nnz(rightofpref)+nnz(leftofpref))/meanReoMag;



pd = adjustForPolarHistogram(pdir(nhs >= sno.minHS), deg2rad(reotx));  
[~,my,se] = meanyvsx(pd, nhs(nhs >= sno.minHS),  deg2rad(reotxe));
meanNumHS = my;
meanNumHS_eb = se;


[~,my,se] = eset.meanField2vsField1(fieldname, 'speed', deg2rad(txe), 'polar', true, 'runs','validname', sno.validname, 'validoperation', sno.validoperation);
speedVsDir = my;
speedVsDir_eb = se;



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
 
    headSwingAcceptanceRateRight(k) = mean(hsacc(bin == k & toright));
    headSwingAcceptanceRateRight_eb(k) = sqrt(headSwingAcceptanceRateRight(k) * (1-headSwingAcceptanceRateRight(k))) / sqrt(nnz(bin ==k & toright));

    headSwingAcceptanceRateLeft(k) = mean(hsacc(bin == k & toleft));
    headSwingAcceptanceRateLeft_eb(k) = sqrt(headSwingAcceptanceRateLeft(k) * (1-headSwingAcceptanceRateLeft(k))) / sqrt(nnz(bin ==k & toleft));

    headSwingRejectionBiasRight(k) = nnz(~hsacc(bin == k & toright))/nnz(~hsacc(bin == k)) / mean(toright(bin == k));
    headSwingRejectionBiasRight_eb(k) = sqrt(nnz(~hsacc(bin == k & toright)))/nnz(~hsacc(bin == k)) / mean(toright(bin == k)); %this is a kludge for now
    
    headSwingRejectionBiasLeft(k) = nnz(~hsacc(bin == k & toleft))/nnz(~hsacc(bin == k)) / mean(toleft(bin == k));
    headSwingRejectionBiasLeft_eb(k) = sqrt(nnz(~hsacc(bin == k & toleft)))/nnz(~hsacc(bin == k)) / mean(toleft(bin == k)); %this is a kludge for now
    
    meanRejectedHeadSwingDir(k) = mean(hsmt(bin == k & ~hsacc));
    meanRejectedHeadSwingDir_eb(k) = std(hsmt(bin == k & ~hsacc))/sqrt(nnz(bin == k & ~hsacc));

    meanAcceptedHeadSwingDir(k) = mean(hsmt(bin == k & hsacc));
    meanAcceptedHeadSwingDir_eb(k) = std(hsmt(bin == k & hsacc))/sqrt(nnz(bin == k & hsacc));
end

perp = abs(sind(sno.preferredDirection - tdeg)) > cosd(hbw);

headSwingAcceptanceRateTowards = mean(hsacc(perp & towards));
std(hsacc(perp&towards))/sqrt(nnz(perp&towards));
sqrt(headSwingAcceptanceRateTowards * (1-headSwingAcceptanceRateTowards)) / sqrt(nnz(perp & towards));
headSwingAcceptanceRateTowards_eb = sqrt(headSwingAcceptanceRateTowards * (1-headSwingAcceptanceRateTowards)) / sqrt(nnz(perp & towards));

headSwingAcceptanceRateAway = mean(hsacc(perp & ~towards));
headSwingAcceptanceRateAway_eb = sqrt(headSwingAcceptanceRateAway * (1-headSwingAcceptanceRateAway)) / sqrt(nnz(perp & ~towards));

headSwingAcceptanceRateIndex = (headSwingAcceptanceRateTowards - headSwingAcceptanceRateAway)/(headSwingAcceptanceRateTowards + headSwingAcceptanceRateAway);
headSwingAcceptanceRateIndex_eb = sqrt(var(hsacc(perp & towards)) * nnz(perp & towards) + var(hsacc(perp & ~towards)) * nnz(perp & ~towards))/nnz(perp)/(headSwingAcceptanceRateTowards + headSwingAcceptanceRateAway);



topref = cos(hstail - deg2rad(sno.preferredDirection)) > 1/sqrt(2);
frompref = cos(hstail - deg2rad(sno.preferredDirection)) < -1/sqrt(2);
leftofpref = sin(hstail - deg2rad(sno.preferredDirection)) > 1/sqrt(2);
rightofpref = sin(hstail - deg2rad(sno.preferredDirection)) < -1/sqrt(2);
hsmt = adjustForPolarHistogram(hsmt, deg2rad(hsdtx));
hs_all_dist{1} = hist(hsmt(topref), deg2rad(hsdtx));
hs_all_dist{2} = hist(hsmt(frompref), deg2rad(hsdtx));
hs_all_dist{3} = hist(hsmt(leftofpref), deg2rad(hsdtx));
hs_all_dist{4} = hist(hsmt(rightofpref), deg2rad(hsdtx));

hs_acc_dist{1} = hist(hsmt(topref & hsacc), deg2rad(hsdtx));
hs_acc_dist{2} = hist(hsmt(frompref & hsacc), deg2rad(hsdtx));
hs_acc_dist{3} = hist(hsmt(leftofpref & hsacc), deg2rad(hsdtx));
hs_acc_dist{4} = hist(hsmt(rightofpref & hsacc), deg2rad(hsdtx));


hs_rej_dist{1} = hist(hsmt(topref & ~hsacc), deg2rad(hsdtx));
hs_rej_dist{2} = hist(hsmt(frompref & ~hsacc), deg2rad(hsdtx));
hs_rej_dist{3} = hist(hsmt(leftofpref & ~hsacc), deg2rad(hsdtx));
hs_rej_dist{4} = hist(hsmt(rightofpref & ~hsacc), deg2rad(hsdtx));



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

    
hsacc = eset.gatherSubField('firsths', 'accepted');
hsacc = hsacc(htv);

topref = cos(hstail - deg2rad(sno.preferredDirection)) > 1/sqrt(2);
frompref = cos(hstail - deg2rad(sno.preferredDirection)) < -1/sqrt(2);
leftofpref = sin(hstail - deg2rad(sno.preferredDirection)) > 1/sqrt(2);
rightofpref = sin(hstail - deg2rad(sno.preferredDirection)) < -1/sqrt(2);
hsmt = adjustForPolarHistogram(hsmt, deg2rad(hsdtx));
firsths_all_dist{1} = hist(hsmt(topref), deg2rad(hsdtx));
firsths_all_dist{2} = hist(hsmt(frompref), deg2rad(hsdtx));
firsths_all_dist{3} = hist(hsmt(leftofpref), deg2rad(hsdtx));
firsths_all_dist{4} = hist(hsmt(rightofpref), deg2rad(hsdtx));

firsths_acc_dist{1} = hist(hsmt(topref & hsacc), deg2rad(hsdtx));
firsths_acc_dist{2} = hist(hsmt(frompref & hsacc), deg2rad(hsdtx));
firsths_acc_dist{3} = hist(hsmt(leftofpref & hsacc), deg2rad(hsdtx));
firsths_acc_dist{4} = hist(hsmt(rightofpref & hsacc), deg2rad(hsdtx));


firsths_rej_dist{1} = hist(hsmt(topref & ~hsacc), deg2rad(hsdtx));
firsths_rej_dist{2} = hist(hsmt(frompref & ~hsacc), deg2rad(hsdtx));
firsths_rej_dist{3} = hist(hsmt(leftofpref & ~hsacc), deg2rad(hsdtx));
firsths_rej_dist{4} = hist(hsmt(rightofpref & ~hsacc), deg2rad(hsdtx));


[autocorr_theta_withinruns, autocorr_np_withinruns,autocorr_tx_withinruns] = eset.autocorrelate('theta', 'isangle', true, 'withinRuns', true,'timerange',sno.autocorr_timerange);

[autocorr_theta_inruns, autocorr_np_inruns,autocorr_tx_inruns] = eset.autocorrelate('theta', 'isangle', true, 'inRuns', true,'timerange',sno.autocorr_timerange);

fields = {'tx', 'txc', 'txe', 'hstx', 'hstxc','hstxe', 'dtx', 'dtxc', 'hsdtx', 'hsdtx', 'thetahist', 'reohist', 'runStartDirectionHist', 'headSwingAcceptanceRateTowards',...
    'headSwingAcceptanceRateAway', 'reoMag', 'reoMag_eb', 'reoStd', 'meanNumHS', 'speedVsDir','firstHSBias', 'firstHSDir',... %reoStd is deprecated
    'headSwingAcceptanceRateLeft', 'headSwingAcceptanceRateRight','meanrunchange','instantaneousdthetavstheta',...
    'thetahist_eb', 'reohist_eb', 'runStartDirectionHist_eb', 'headSwingAcceptanceRateTowards_eb',...
    'headSwingAcceptanceRateAway_eb', 'speedVsDir_eb','firstHSBias_eb', 'firstHSDir_eb',...
    'headSwingAcceptanceRateLeft_eb', 'headSwingAcceptanceRateRight_eb','meanrunchange_eb','instantaneousdthetavstheta_eb',...
    'allHSDir', 'allHSDir_eb', 'reoDir', 'reoDir_eb', 'reotx', 'reotxe', 'reotxc','headSwingRejectionBiasRight','headSwingRejectionBiasRight_eb','headSwingRejectionBiasLeft','headSwingRejectionBiasLeft_eb',...
    'meanRejectedHeadSwingDir', 'meanRejectedHeadSwingDir_eb','meanAcceptedHeadSwingDir','meanAcceptedHeadSwingDir_eb', 'firstHSMeanDir', 'firstHSMeanDir_eb',...
    'reoMagIndex', 'reoMagIndex_eb','reoMagIndex_ci','meanReoMag', 'reoDirIndex', 'reoDirIndex_eb','reoDirIndex_ci','reoRateIndex','reoRateIndex_ci','reoRateIndex_eb','reoRateIndexPerp','reoRateIndexPerp_eb','headSwingAcceptanceRateIndex','headSwingAcceptanceRateIndex_eb',...
    'navind', 'navind_eb','reo_prevdir', 'reo_dtheta','reo_rate_model','reo_dir_model','reoMagIndex_old', 'reoMagIndex_old_eb','reoDirIndex_old', 'reoDirIndex_old_eb',...
    'reobasedirections','reo_dtheta_dist_hs', 'reo_dtheta_dist_hs_eb', 'reo_dtheta_dist_nohs', 'reo_dtheta_dist_nohs_eb','reo_dtheta_fit', 'reo_dtheta_fit_m','reo_dtheta_fit_r','reo_dtheta_fit_l','dtxf', ...
    'hs_all_dist', 'hs_acc_dist', 'hs_rej_dist','firsths_all_dist', 'firsths_acc_dist', 'firsths_rej_dist', 'navind_expt', 'navind_expt_eb', ...
    'autocorr_tx_withinruns', 'autocorr_theta_withinruns', 'autocorr_np_withinruns', 'autocorr_tx_inruns', 'autocorr_theta_inruns', 'autocorr_np_inruns','reo_numHS'};%,...
 %   'reo_dir_model_fixed_theta_0', 'reoMagIndex_fixed_theta_0', 'reoMagIndex_fixed_theta_0_eb', 'reoDirIndex_fixed_theta_0', 'reoDirIndex_fixed_theta_0_eb'};
    %'reo_dtheta_dist_topref','reo_dtheta_dist_frompref','reo_dtheta_dist_leftofpref','reo_dtheta_dist_rightofpref', 'reo_dtheta_fit_topref','reo_dtheta_fit_frompref','reo_dtheta_fit_leftofpref','reo_dtheta_fit_rightofpref',...
    %'reo_dtheta_dist_topref_eb','reo_dtheta_dist_frompref_eb','reo_dtheta_
    %dist_leftofpref_eb','reo_dtheta_dist_rightofpref_eb','navind_ci' 
 
for j = 1:length(fields)
    try
        eval(['analyzedData.' fields{j} ' = ' fields{j} ';']);
    catch me
        disp(me.getReport);
    end
end



