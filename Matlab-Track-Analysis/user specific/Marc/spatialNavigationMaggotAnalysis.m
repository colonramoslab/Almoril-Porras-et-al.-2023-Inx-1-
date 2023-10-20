function ad = spatialNavigationMaggotAnalysis (eset, spatial_navigation_options)
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
sno.validoperation = func2str(@(x) logical(setNonFiniteToZero(x)));
sno.confidenceLevel = 0.95;
sno.autocorr_timerange = [];
sno.runTimeBinSize = 10;

if (nargin == 0)
    if (nargout >= 1)
        ad = sno;
    end
    disp ('spatialNavigationMaggotAnalysis (esets, spatial_navigation_options)');
    return;
end
existsAndDefault('spatial_navigation_options', []);

if (length(eset) > 1)
    for j = 1:length(eset)
        disp (['start ' num2str(j)]); ts1 = tic;
        ad(j) = spatialNavigationMaggotAnalysis(eset(j), spatial_navigation_options);  %#ok<AGROW>
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
    return; %#ok<UNRCH>
end

ad.sno = sno;
if (ischar(sno.validoperation))
    sno.validoperation = str2func(sno.validoperation);
end

ad.txf = -180:1:180;
ad.tx = (-180:sno.angleBinSize:(180-sno.angleBinSize)) + sno.preferredDirection;
ad.txc = (-180:sno.angleBinSize:(180)) + sno.preferredDirection;

txrad = deg2rad(ad.txf);
bsa = deg2rad(sno.angleBinSize);
%{

dtx = (-180:sno.dtAngleBinSize:(180-sno.dtAngleBinSize));
dtxc = (-180:sno.dtAngleBinSize:180);

hsdtx = sort([(sno.minHSTheta - sno.hsdtAngleBinSize/2):(-sno.hsdtAngleBinSize):0 (sno.minHSTheta + sno.hsdtAngleBinSize/2):(sno.hsdtAngleBinSize):180]);
hsdtx = sort([-hsdtx hsdtx]);
%hsdtx = hsdtxc(1:(end-1));
%}

%calculate statistics of runs including change within runs and run length
%-------------------------------------------------------------------------%
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
changeLessThan90 = logical(abs(dt) < pi/2);

[~,my, se] = meanyvsx_slidingwindow(rs(changeLessThan90), dt(changeLessThan90), txrad, bsa, 'step', true);
ad.meanrunchange_step = my;
ad.meanrunchange_step_eb = se;
[~,my, se] = meanyvsx_slidingwindow(rs(changeLessThan90), dt(changeLessThan90), txrad, bsa, 'gaussian', true);
ad.meanrunchange_gauss = my;
ad.meanrunchange_gauss_eb = se;


ad.runStartTheta = rs;
ad.runEndTheta = re;
ad.runMeanTheta = rmt;
ad.runTime = rt;
minruntime = median(eset.gatherSubField('so', 'minRunTime'));

ad.runTimeAxis = minruntime + (sno.runTimeBinSize/2)+ (0:sno.runTimeBinSize:(sno.runTimeBinSize*ceil(600/sno.runTimeBinSize)));
rth = hist(rt(cos(rmt - deg2rad(sno.preferredDirection)) > 1/sqrt(2)), ad.runTimeAxis);
ad.runTimeHistTowards = (rth/sum(rth));
ad.runTimeHistTowards_eb = sqrt((rth/sum(rth)).*(1-(rth/sum(rth)))./sum(rth));
rth = hist(rt(-cos(rmt - deg2rad(sno.preferredDirection)) > 1/sqrt(2)), ad.runTimeAxis);
ad.runTimeHistAway = (rth/sum(rth));
ad.runTimeHistAway_eb = sqrt((rth/sum(rth)).*(1-(rth/sum(rth)))./sum(rth));


%calculate instantaneous change within runs vs. heading angle
%-----------------------------------------------------------------------%
if(~isempty(sno.relativeDirField))
    fieldname = sno.relativeDirField;
else
    fieldname = 'theta';
end
[~,my,se]= eset.meanField2vsField1_slidingwindow(fieldname, 'lrdtheta', txrad,bsa,'step', 'runs', 'polar', true,'validname', sno.validname, 'validoperation', sno.validoperation);
ad.instantaneousdthetavstheta_step = my;
ad.instantaneousdthetavstheta_step_eb = se;
[~,my,se]= eset.meanField2vsField1_slidingwindow(fieldname, 'lrdtheta', txrad,bsa,'gaussian', 'runs', 'polar', true,'validname', sno.validname, 'validoperation', sno.validoperation);
ad.instantaneousdthetavstheta_gauss = my;
ad.instantaneousdthetavstheta_gauss_eb = se;


%calculate speed vs. previous angle
%[~,my,se] = eset.meanField2vsField1(fieldname, 'speed', deg2rad(txe), 'polar', true, 'runs','validname', sno.validname, 'validoperation', sno.validoperation);
[~,my,se]= eset.meanField2vsField1_slidingwindow(fieldname, 'speed', txrad,bsa,'step', 'runs', 'polar', true,'validname', sno.validname, 'validoperation', sno.validoperation);
ad.speedVsDir_step = my;
ad.speedVsDir_step_eb = se;
[~,my,se]= eset.meanField2vsField1_slidingwindow(fieldname, 'speed', txrad,bsa,'gaussian', 'runs', 'polar', true,'validname', sno.validname, 'validoperation', sno.validoperation);
ad.speedVsDir_gauss = my;
ad.speedVsDir_gauss_eb = se;

%speedVsDir = my;
%speedVsDir_eb = se;



%calculate navigation index for eset as a whole and for individual
%experiments
%------------------------------------------------------------------------%
v = eset.gatherField('vel', 'validname', sno.validname, 'validoperation', sno.validoperation);
M = [cosd(sno.preferredDirection) sind(sno.preferredDirection); -sind(sno.preferredDirection) cosd(sno.preferredDirection)];
v = M*v;
s = sqrt(sum(v.^2));
ad.navind = mean(v,2)/mean(s);
ad.navind_eb = std(v,0,2)/mean(s)/sqrt(length(s)*eset.expt(1).track(1).dr.interpTime/(2*eset.autocorr_tau));
for k = 1:length(eset.expt)
    v = eset.expt(k).gatherField('vel', 'validname', sno.validname, 'validoperation', sno.validoperation);
    M = [cosd(sno.preferredDirection) sind(sno.preferredDirection); -sind(sno.preferredDirection) cosd(sno.preferredDirection)];
    v = M*v;
    s = sqrt(sum(v.^2));
    ad.navind_expt(:,k) = mean(v,2)/mean(s);
    ad.navind_expt_eb(:,k) = std(v,0,2)/mean(s)/sqrt(length(s)*eset.expt(1).track(1).dr.interpTime/(2*eset.autocorr_tau));
end

%calculate histogram of instantaneous direction
%------------------------------------------------------------------------%
[h,eb] = eset.makeHistogram(fieldname, deg2rad(ad.tx), 'runs', 'polar', true,'validname', sno.validname, 'validoperation', sno.validoperation);
ad.thetahist_eb = eb ./ mean(h);
ad.thetahist = h./mean(h);


%calculate reorientation rate vs. direction
%------------------------------------------------------------------------%
[h,eb] = eset.makeReorientationHistogram(fieldname, deg2rad(ad.tx), 'minHS',sno.minHS, 'polar', true,'validname', sno.validname, 'validoperation', sno.validoperation);
ad.reohist = h;
ad.reohist_eb = eb;


%calculate statistics of reorientations
%-------------------------------------------------------------------------

%initial calculations
ad.reotx = (-180:sno.reoBinSize:(180-sno.reoBinSize)) + sno.preferredDirection;
ad.reotxc = (-180:sno.reoBinSize:(180)) + sno.preferredDirection;
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
%rmag = abs(reodir);
%rmag = acos(cos(ndir - pdir));

allth = eset.gatherField(fieldname, 'runs','validname',sno.validname, 'validoperation', sno.validoperation);


%calculate reorientation models
ad.reo_rate_model = fitReorientationRateTheta(eset, allth, pdir);
ad.reoRateIndex = ad.reo_rate_model.params(2)/ad.reo_rate_model.params(1) * cos(ad.reo_rate_model.params(3) - deg2rad(sno.preferredDirection));
[minx maxx] = confidenceRange(ad.reo_rate_model.params, ad.reo_rate_model.hess, @(x) x(2,:)./x(1,:).*cos(x(3) - deg2rad(sno.preferredDirection)), sno.confidenceLevel);
ad.reoRateIndex_ci = [minx maxx];
ad.reoRateIndex_eb =  stdOfFunVal(ad.reo_rate_model.params, ad.reo_rate_model.hess, @(x) x(2,:)./x(1,:).*cos(x(3) - deg2rad(sno.preferredDirection)));

ad.reoRateIndexPerp = ad.reo_rate_model.params(2)/ad.reo_rate_model.params(1) * sin(ad.reo_rate_model.params(3) - deg2rad(sno.preferredDirection));
[minx maxx] = confidenceRange(ad.reo_rate_model.params, ad.reo_rate_model.hess, @(x) x(2,:)./x(1,:).*sin(x(3) - deg2rad(sno.preferredDirection)), sno.confidenceLevel);
ad.reoRateIndexPerp_eb = 0.5*(maxx-minx);

reo_prevdir = pdir;
reo_dtheta = reodir;%since these were previously used interchangeably, stick them in here, to avoid bugs
%reo_numHS = nhs; 
ad.reo_prevdir = pdir;
ad.reo_dtheta = reodir;
ad.reo_numHS = nhs;
ad.reo_dir_model = fitReorientationAngleDistributionWithHSInfo(pdir, ndir, nhs);

ad.nixedParams = {{}, {'C'}, {'B', 'E'}, {'D', 'E'}, {'B','C','E'}};
ad.altmodeldescriptions = {'orig model', 'no left/right bias', 'no size bias', 'no skew', 'no bias'};
for j = 1:length(ad.nixedParams)
    clear startValues fixedValues
    for k = 1:length(ad.reo_dir_model.params)
        startValues.(ad.reo_dir_model.paramkey{k}) = ad.reo_dir_model.params(k);
    end
    fixedValues = [];
    for k = 1:length(ad.nixedParams{j})
        startValues.(ad.nixedParams{j}{k}) = 0;
        fixedValues.(ad.nixedParams{j}{k}) = 0;
    end
    m1 = fitReorientationAngleDistribution(pdir(nhs > 0), diff(unwrap([pdir(nhs>0);ndir(nhs>0)])), fixedValues,[],false); %nix start values
    m2 = fitReorientationAngleDistribution(pdir(nhs > 0), diff(unwrap([pdir(nhs>0);ndir(nhs>0)])), fixedValues, startValues,false); %with start values
    if (m1.logLikelihood > m2.logLikelihood)
        ad.altmodel(j) = m1;
    else
        ad.altmodel(j) = m2;
    end
end
ad.deltaLL = [ad.altmodel.logLikelihood] - ad.reo_dir_model.logLikelihood;
% 'A'    'B'    'C'    'U'    'sigma1'    'sigma2'    'skew'    'theta_0'
% P(dt | theta) = a*N(0,sigma1,dt) + (1/2 - a/2 + g) * S(U, sigma2, skew, dt) +
%    (1/2 - a/2 - g) * S(U, sigma2, skew, -dt)
%
%a(theta) = A + B*cos(theta - theta_0); (|B| < A);
%g(theta) = C(sin(theta - theta_0)) 
%other fit parameters, U, sigma1, sigma2, skew, theta_0 
dirIndexFun = @(x) 2*x(3).*cos(deg2rad(sno.preferredDirection) - x(7));
ad.reoDirIndex = dirIndexFun(ad.reo_dir_model.params);
[minx maxx] = confidenceRange(ad.reo_dir_model.params, ad.reo_dir_model.hessian, dirIndexFun, sno.confidenceLevel);
ad.reoDirIndex_ci = [minx maxx];
ad.reoDirIndex_eb = stdOfFunVal(ad.reo_dir_model.params, ad.reo_dir_model.hessian, dirIndexFun);

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
ad.reoMagIndex = magIndexFun(ad.reo_dir_model.params);
[minx maxx] = confidenceRange(ad.reo_dir_model.params, ad.reo_dir_model.hessian, magIndexFun, sno.confidenceLevel);
ad.reoMagIndex_ci = [minx maxx];
ad.reoMagIndex_eb = stdOfFunVal(ad.reo_dir_model.params, ad.reo_dir_model.hessian, magIndexFun);


%calculate statistics of run starts (outcome of reorientations)

nd = adjustForPolarHistogram(ndir(nhs >= sno.minHS), deg2rad(ad.tx));  
runStartDirectionHist = hist(nd, deg2rad(ad.tx));
h = runStartDirectionHist;
n = sum(h);

ad.runStartDirectionHist_eb = sqrt(h/n).*sqrt(1-h/n).*sqrt(n) ./ mean(runStartDirectionHist);
ad.runStartDirectionHist = runStartDirectionHist./mean(runStartDirectionHist);
 

%calculate distributions of reorientations for given previous headings
topref = cos(reo_prevdir - deg2rad(sno.preferredDirection)) > 1/sqrt(2);
frompref = cos(reo_prevdir - deg2rad(sno.preferredDirection)) < -1/sqrt(2);
leftofpref = sin(reo_prevdir - deg2rad(sno.preferredDirection)) > 1/sqrt(2);
rightofpref = sin(reo_prevdir - deg2rad(sno.preferredDirection)) < -1/sqrt(2);

indset = {topref, frompref, leftofpref, rightofpref};


ad.dtx = (-180:sno.dtAngleBinSize:(180-sno.dtAngleBinSize));
ad.dtxc = (-180:sno.dtAngleBinSize:180);

dtxf = deg2rad(-180:1:180);
%ddtxf = deg2rad(diff(dtxf(1:2)));
dt = adjustForPolarHistogram(reo_dtheta, deg2rad(ad.dtx)); 
for j = 1:length(indset)
    ad.reo_dtheta_dist_hs{j} = hist(dt(indset{j} & nhs > 0), deg2rad(ad.dtx));
    ad.reo_dtheta_dist_hs_eb{j} = sqrt(ad.reo_dtheta_dist_hs{j}.*(1-ad.reo_dtheta_dist_hs{j}/nnz(indset{j} & nhs > 0)));
    ad.reo_dtheta_dist_nohs{j} = hist(dt(indset{j} & nhs == 0), deg2rad(ad.dtx));
    ad.reo_dtheta_dist_nohs_eb{j} = sqrt(ad.reo_dtheta_dist_nohs{j}.*(1-ad.reo_dtheta_dist_nohs{j}/nnz(indset{j} & nhs == 0)));
end

%calculate predicted distributions from reorientation model
for j = 1:length(indset)
%    reo_dtheta_fit_m{j} = zeros(size(dtxf));
    ad.reo_dtheta_fit_r{j} = zeros(size(dtxf));
    ad.reo_dtheta_fit_l{j} = zeros(size(dtxf));
    for k = find(indset{j} & nhs >= sno.minHS)
        [r,l] = ad.reo_dir_model.pdfOfdThetaBreakdown(ad.reo_dir_model.params, reo_prevdir(k), dtxf);
        %reo_dtheta_fit_m{j} = reo_dtheta_fit_m{j} + m*deg2rad(sno.dtAngleBinSize);
        ad.reo_dtheta_fit_r{j} = ad.reo_dtheta_fit_r{j} + r*deg2rad(sno.dtAngleBinSize);
        ad.reo_dtheta_fit_l{j} = ad.reo_dtheta_fit_l{j} + l*deg2rad(sno.dtAngleBinSize);
    end
    ad.reo_dtheta_fit{j} = ad.reo_dtheta_fit_r{j} + ad.reo_dtheta_fit_l{j};
end
ad.dtxf = rad2deg(dtxf);
ad.reobasedirections = [mod(sno.preferredDirection + 180, 360) - 180, mod(sno.preferredDirection, 360) - 180, mod(sno.preferredDirection + 270, 360) - 180, mod(sno.preferredDirection + 90, 360) - 180];

%calculate mean reorientation magnitude, direction, and number of head sweeps vs. previous heading

pd = reo_prevdir;
[~,my,se] = meanyvsx_slidingwindow(pd(nhs >= sno.minHS), reodir(nhs >= sno.minHS).^2,  txrad, deg2rad(sno.reoBinSize), 'step', true);
ad.reoMag_step = my;
ad.reoMag_step_eb = se; 
[~,my,se] = meanyvsx_slidingwindow(pd(nhs >= sno.minHS), reodir(nhs >= sno.minHS).^2,  txrad, deg2rad(sno.reoBinSize), 'gaussian', true);
ad.reoMag_gauss = my;
ad.reoMag_gauss_eb = se; 


[~,my,se] = meanyvsx_slidingwindow(pd(nhs >= sno.minHS), reodir(nhs >= sno.minHS),  txrad, deg2rad(sno.reoBinSize), 'step', true);
ad.reoDir_step = my;
ad.reoDir_step_eb = se; 
[~,my,se] = meanyvsx_slidingwindow(pd(nhs >= sno.minHS), reodir(nhs >= sno.minHS),  txrad, deg2rad(sno.reoBinSize), 'gaussian', true);
ad.reoDir_gauss = my;
ad.reoDir_gauss_eb = se; 

[~,my,se] = meanyvsx_slidingwindow(pd(nhs >= sno.minHS), nhs(nhs >= sno.minHS),  txrad, deg2rad(sno.reoBinSize), 'step', true);
ad.meanNumHS_step = my;
ad.meanNumHS_step_eb = se; 
[~,my,se] = meanyvsx_slidingwindow(pd(nhs >= sno.minHS), nhs(nhs >= sno.minHS),  txrad, deg2rad(sno.reoBinSize), 'gaussian', true);
ad.meanNumHS_gauss = my;
ad.meanNumHS_gauss_eb = se; 


%on to head sweeps

ad.hstx = (-180:sno.hsBinSpacing:(180-sno.hsBinSpacing)) + sno.preferredDirection;
ad.hstxc = (-180:sno.hsBinSpacing:(180)) + sno.preferredDirection;
hstx = ad.hstx;
%hstxc = ad.hstxc;
hstxe = ad.hstxc - sno.hsBinSpacing/2;
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
    
ad.allHSDir = my;
ad.allHSDir_eb = se;
for k = 1:length(hstx)
 
    ad.headSwingAcceptanceRateRight(k) = mean(hsacc(bin == k & toright));
    ad.headSwingAcceptanceRateRight_eb(k) = sqrt(ad.headSwingAcceptanceRateRight(k) * (1-ad.headSwingAcceptanceRateRight(k))) / sqrt(nnz(bin ==k & toright));

    ad.headSwingAcceptanceRateLeft(k) = mean(hsacc(bin == k & toleft));
    ad.headSwingAcceptanceRateLeft_eb(k) = sqrt(ad.headSwingAcceptanceRateLeft(k) * (1-ad.headSwingAcceptanceRateLeft(k))) / sqrt(nnz(bin ==k & toleft));

    ad.headSwingRejectionBiasRight(k) = nnz(~hsacc(bin == k & toright))/nnz(~hsacc(bin == k)) / mean(toright(bin == k));
    ad.headSwingRejectionBiasRight_eb(k) = sqrt(nnz(~hsacc(bin == k & toright)))/nnz(~hsacc(bin == k)) / mean(toright(bin == k)); %this is a kludge for now
    
    ad.headSwingRejectionBiasLeft(k) = nnz(~hsacc(bin == k & toleft))/nnz(~hsacc(bin == k)) / mean(toleft(bin == k));
    ad.headSwingRejectionBiasLeft_eb(k) = sqrt(nnz(~hsacc(bin == k & toleft)))/nnz(~hsacc(bin == k)) / mean(toleft(bin == k)); %this is a kludge for now
    
    ad.meanRejectedHeadSwingDir(k) = mean(hsmt(bin == k & ~hsacc));
    ad.meanRejectedHeadSwingDir_eb(k) = std(hsmt(bin == k & ~hsacc))/sqrt(nnz(bin == k & ~hsacc));

    ad.meanAcceptedHeadSwingDir(k) = mean(hsmt(bin == k & hsacc));
    ad.meanAcceptedHeadSwingDir_eb(k) = std(hsmt(bin == k & hsacc))/sqrt(nnz(bin == k & hsacc));
end

perp = abs(sind(sno.preferredDirection - tdeg)) > cosd(hbw);

ad.headSwingAcceptanceRateTowards = mean(hsacc(perp & towards));
% std(hsacc(perp&towards))/sqrt(nnz(perp&towards));
% sqrt(headSwingAcceptanceRateTowards * (1-headSwingAcceptanceRateTowards)) / sqrt(nnz(perp & towards));
ad.headSwingAcceptanceRateTowards_eb = sqrt(ad.headSwingAcceptanceRateTowards * (1-ad.headSwingAcceptanceRateTowards)) / sqrt(nnz(perp & towards));

ad.headSwingAcceptanceRateAway = mean(hsacc(perp & ~towards));
ad.headSwingAcceptanceRateAway_eb = sqrt(ad.headSwingAcceptanceRateAway * (1-ad.headSwingAcceptanceRateAway)) / sqrt(nnz(perp & ~towards));

ad.headSwingAcceptanceRateIndex = (ad.headSwingAcceptanceRateTowards - ad.headSwingAcceptanceRateAway)/(ad.headSwingAcceptanceRateTowards + ad.headSwingAcceptanceRateAway);
ad.headSwingAcceptanceRateIndex_eb = sqrt(var(hsacc(perp & towards)) * nnz(perp & towards) + var(hsacc(perp & ~towards)) * nnz(perp & ~towards))/nnz(perp)/(ad.headSwingAcceptanceRateTowards + ad.headSwingAcceptanceRateAway);


%distribution of head swing angles
hsdtx = sort([(sno.minHSTheta - sno.hsdtAngleBinSize/2):(-sno.hsdtAngleBinSize):0 (sno.minHSTheta + sno.hsdtAngleBinSize/2):(sno.hsdtAngleBinSize):180]);
hsdtx = sort([-hsdtx hsdtx]);
ad.hsdtx = hsdtx;
topref = cos(hstail - deg2rad(sno.preferredDirection)) > 1/sqrt(2);
frompref = cos(hstail - deg2rad(sno.preferredDirection)) < -1/sqrt(2);
leftofpref = sin(hstail - deg2rad(sno.preferredDirection)) > 1/sqrt(2);
rightofpref = sin(hstail - deg2rad(sno.preferredDirection)) < -1/sqrt(2);
hsmt = adjustForPolarHistogram(hsmt, deg2rad(hsdtx));
ad.hs_all_dist{1} = hist(hsmt(topref), deg2rad(hsdtx));
ad.hs_all_dist{2} = hist(hsmt(frompref), deg2rad(hsdtx));
ad.hs_all_dist{3} = hist(hsmt(leftofpref), deg2rad(hsdtx));
ad.hs_all_dist{4} = hist(hsmt(rightofpref), deg2rad(hsdtx));

ad.hs_acc_dist{1} = hist(hsmt(topref & hsacc), deg2rad(hsdtx));
ad.hs_acc_dist{2} = hist(hsmt(frompref & hsacc), deg2rad(hsdtx));
ad.hs_acc_dist{3} = hist(hsmt(leftofpref & hsacc), deg2rad(hsdtx));
ad.hs_acc_dist{4} = hist(hsmt(rightofpref & hsacc), deg2rad(hsdtx));


ad.hs_rej_dist{1} = hist(hsmt(topref & ~hsacc), deg2rad(hsdtx));
ad.hs_rej_dist{2} = hist(hsmt(frompref & ~hsacc), deg2rad(hsdtx));
ad.hs_rej_dist{3} = hist(hsmt(leftofpref & ~hsacc), deg2rad(hsdtx));
ad.hs_rej_dist{4} = hist(hsmt(rightofpref & ~hsacc), deg2rad(hsdtx));


%repeat head sweep analysis, but only for first head sweeps

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
ad.firstHSBias = my;
ad.firstHSBias_eb =se;
[~,my,se] = meanyvsx(td, toleft, deg2rad(hstxe));
ad.firstHSDir = my;
ad.firstHSDir_eb = se;

[~,my,se] = meanyvsx(td, hsmt, deg2rad(hstxe));
ad.firstHSMeanDir = my;
ad.firstHSMeanDir_eb = se;

    
hsacc = eset.gatherSubField('firsths', 'accepted');
hsacc = hsacc(htv);

topref = cos(hstail - deg2rad(sno.preferredDirection)) > 1/sqrt(2);
frompref = cos(hstail - deg2rad(sno.preferredDirection)) < -1/sqrt(2);
leftofpref = sin(hstail - deg2rad(sno.preferredDirection)) > 1/sqrt(2);
rightofpref = sin(hstail - deg2rad(sno.preferredDirection)) < -1/sqrt(2);
hsmt = adjustForPolarHistogram(hsmt, deg2rad(hsdtx));
ad.firsths_all_dist{1} = hist(hsmt(topref), deg2rad(hsdtx));
ad.firsths_all_dist{2} = hist(hsmt(frompref), deg2rad(hsdtx));
ad.firsths_all_dist{3} = hist(hsmt(leftofpref), deg2rad(hsdtx));
ad.firsths_all_dist{4} = hist(hsmt(rightofpref), deg2rad(hsdtx));

ad.firsths_acc_dist{1} = hist(hsmt(topref & hsacc), deg2rad(hsdtx));
ad.firsths_acc_dist{2} = hist(hsmt(frompref & hsacc), deg2rad(hsdtx));
ad.firsths_acc_dist{3} = hist(hsmt(leftofpref & hsacc), deg2rad(hsdtx));
ad.firsths_acc_dist{4} = hist(hsmt(rightofpref & hsacc), deg2rad(hsdtx));


ad.firsths_rej_dist{1} = hist(hsmt(topref & ~hsacc), deg2rad(hsdtx));
ad.firsths_rej_dist{2} = hist(hsmt(frompref & ~hsacc), deg2rad(hsdtx));
ad.firsths_rej_dist{3} = hist(hsmt(leftofpref & ~hsacc), deg2rad(hsdtx));
ad.firsths_rej_dist{4} = hist(hsmt(rightofpref & ~hsacc), deg2rad(hsdtx));


[ad.autocorr_theta_withinruns, ad.autocorr_np_withinruns,ad.autocorr_tx_withinruns] = eset.autocorrelate('theta', 'isangle', true, 'withinRuns', true,'timerange',sno.autocorr_timerange);

[ad.autocorr_theta_inruns, ad.autocorr_np_inruns,ad.autocorr_tx_inruns] = eset.autocorrelate('theta', 'isangle', true, 'inRuns', true,'timerange',sno.autocorr_timerange);

fn = {'meanrunchange','instantaneousdthetavstheta', 'instantaneousdthetavstheta', 'speedVsDir'};
tx = ad.tx;
for j = 1:length(fn)
    yd = ad.([fn{j} '_step']);
    eb = ad.([fn{j} '_step_eb']);
    inds = zeros(size(tx));
    for k = 1:length(tx)
        [~,I] = min(abs(ad.txf - tx(k)));
        inds(k) = I;
    end
    ad.(fn{j}) = yd(inds);
    ad.([fn{j} '_eb']) = eb(inds);
end
fn = {'reoMag','reoDir','meanNumHS'};
tx = ad.reotx;
for j = 1:length(fn)
    yd = ad.([fn{j} '_step']);
    eb = ad.([fn{j} '_step_eb']);
    inds = zeros(size(tx));
    for k = 1:length(tx)
        [~,I] = min(abs(ad.txf - tx(k)));
        inds(k) = I;
    end
    ad.(fn{j}) = yd(inds);
    ad.([fn{j} '_eb']) = eb(inds);
end

% 
% fields = {'tx', 'txc', 'txe', 'hstx', 'hstxc','hstxe', 'dtx', 'dtxc', 'hsdtx', 'hsdtx', 'thetahist', 'reohist', 'runStartDirectionHist', 'headSwingAcceptanceRateTowards',...
%     'headSwingAcceptanceRateAway', 'reoMag', 'reoMag_eb', 'reoStd', 'meanNumHS', 'speedVsDir','firstHSBias', 'firstHSDir',... %reoStd is deprecated
%     'headSwingAcceptanceRateLeft', 'headSwingAcceptanceRateRight','meanrunchange','instantaneousdthetavstheta',...
%     'thetahist_eb', 'reohist_eb', 'runStartDirectionHist_eb', 'headSwingAcceptanceRateTowards_eb',...
%     'headSwingAcceptanceRateAway_eb', 'speedVsDir_eb','firstHSBias_eb', 'firstHSDir_eb',...
%     'headSwingAcceptanceRateLeft_eb', 'headSwingAcceptanceRateRight_eb','meanrunchange_eb','instantaneousdthetavstheta_eb',...
%     'allHSDir', 'allHSDir_eb', 'reoDir', 'reoDir_eb', 'reotx', 'reotxe', 'reotxc','headSwingRejectionBiasRight','headSwingRejectionBiasRight_eb','headSwingRejectionBiasLeft','headSwingRejectionBiasLeft_eb',...
%     'meanRejectedHeadSwingDir', 'meanRejectedHeadSwingDir_eb','meanAcceptedHeadSwingDir','meanAcceptedHeadSwingDir_eb', 'firstHSMeanDir', 'firstHSMeanDir_eb',...
%     'reoMagIndex', 'reoMagIndex_eb','reoMagIndex_ci','meanReoMag', 'reoDirIndex', 'reoDirIndex_eb','reoDirIndex_ci','reoRateIndex','reoRateIndex_ci','reoRateIndex_eb','reoRateIndexPerp','reoRateIndexPerp_eb','headSwingAcceptanceRateIndex','headSwingAcceptanceRateIndex_eb',...
%     'navind', 'navind_eb','reo_prevdir', 'reo_dtheta','reo_rate_model','reo_dir_model','reoMagIndex_old', 'reoMagIndex_old_eb','reoDirIndex_old', 'reoDirIndex_old_eb',...
%     'reobasedirections','reo_dtheta_dist_hs', 'reo_dtheta_dist_hs_eb', 'reo_dtheta_dist_nohs', 'reo_dtheta_dist_nohs_eb','reo_dtheta_fit', 'reo_dtheta_fit_m','reo_dtheta_fit_r','reo_dtheta_fit_l','dtxf', ...
%     'hs_all_dist', 'hs_acc_dist', 'hs_rej_dist','firsths_all_dist', 'firsths_acc_dist', 'firsths_rej_dist', 'navind_expt', 'navind_expt_eb', ...
%     'autocorr_tx_withinruns', 'autocorr_theta_withinruns', 'autocorr_np_withinruns', 'autocorr_tx_inruns', 'autocorr_theta_inruns', 'autocorr_np_inruns','reo_numHS'};%,...
%  %   'reo_dir_model_fixed_theta_0', 'reoMagIndex_fixed_theta_0', 'reoMagIndex_fixed_theta_0_eb', 'reoDirIndex_fixed_theta_0', 'reoDirIndex_fixed_theta_0_eb'};
%     %'reo_dtheta_dist_topref','reo_dtheta_dist_frompref','reo_dtheta_dist_leftofpref','reo_dtheta_dist_rightofpref', 'reo_dtheta_fit_topref','reo_dtheta_fit_frompref','reo_dtheta_fit_leftofpref','reo_dtheta_fit_rightofpref',...
%     %'reo_dtheta_dist_topref_eb','reo_dtheta_dist_frompref_eb','reo_dtheta_
%     %dist_leftofpref_eb','reo_dtheta_dist_rightofpref_eb','navind_ci' 
%  
% for j = 1:length(fields)
%     try
%         eval(['analyzedData.' fields{j} ' = ' fields{j} ';']);
%     catch me
%         disp(me.getReport);
%     end
% end
% 
% 
% 
