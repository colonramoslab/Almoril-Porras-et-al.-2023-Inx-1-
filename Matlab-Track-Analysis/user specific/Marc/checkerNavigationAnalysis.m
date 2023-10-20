function ad = checkerNavigationAnalysis (esets, checker_navigation_options, checkercalcs)
%function ad = checkerNavigationAnalysis (esets, checker_navigation_options, checkercalcs)
% checkner navigation options
% cno.angleBinSize = 30; % in degrees
% cno.preferredDirection = 0;
% cno.hsBinSize = 90;
% cno.minHS = 1;
% cno.distBin = 0.05;
if (nargin > 0 && length(esets) > 1)
    existsAndDefault('checker_navigation_options', []);
    for j = 1:length(esets)
        if (exists('checkercalcs', 'var') && length(checkercalcs) >= j)
            ad(j) = checkerNavigationAnalysis(esets(j), checker_navigation_options, checkercalcs(j)); %#ok<AGROW>
        else
            ad(j) = checkerNavigationAnalysis(esets(j), checker_navigation_options); %#ok<AGROW>
        end
    end
    return;
end

cno.angleBinSize = 30; % in degrees
cno.preferredDirection = 0;
cno.hsBinSize = 90;
cno.minHS = 1;
cno.distBin = 0.05;

if (nargin == 0)
    ad = cno;
    return;
end

eset = esets(1);

if (nargin > 1 && isstruct(checker_navigation_options))
    fn = fieldnames(checker_navigation_options);
    for j = 1:length(fn)
        cno.(fn{j}) = checker_navigation_options.(fn{j});
    end
end

tx = (-180:cno.angleBinSize:(180-cno.angleBinSize)) + cno.preferredDirection;
txc = (-180:cno.angleBinSize:(180)) + cno.preferredDirection;
txe = txc - cno.angleBinSize/2;

ad.tx = tx;
ad.txc = txc;

if (existsAndDefault('checkercalcs', []))
    cc = checkercalcs;
else
    cc = checkerCalculations('eset');
end

fn = fieldnames(cc);
for j = 1:length(fn)
    eval([fn{j} ' = cc.' fn{j} ';']);
end

reo_ttb = adjustForPolarHistogram (reo_ttb, deg2rad(tx)); 
run_ttb = adjustForPolarHistogram (run_ttb, deg2rad(tx));

minpp = (60/eset.expt(1).track(1).dr.interpTime);

%r = eset.gatherField('reorientation');
h1 = hist(reo_ttb(reo_onb & reo_numHS >= cno.minHS), deg2rad(tx));
h2 = hist(run_ttb(run_onb), deg2rad(tx));
ad.reorate_thetabound = (h1./h2) * minpp;
ad.reorate_thetabound_eb = sqrt(h1)./h2* minpp;
ad.reorate_thetabound = ad.reorate_thetabound([1:end 1]);
ad.reorate_thetabound_eb = ad.reorate_thetabound_eb([1:end 1]);

h1 = hist(reo_ttb(reo_onb & reo_numHS == 0), deg2rad(tx));
h2 = hist(run_ttb(run_onb), deg2rad(tx));
ad.pauserate_thetabound = (h1./h2) * minpp;
ad.pauserate_thetabound_eb = sqrt(h1)./h2* minpp;
ad.pauserate_thetabound = ad.pauserate_thetabound([1:end 1]);
ad.pauserate_thetabound_eb = ad.pauserate_thetabound_eb([1:end 1]);

distx = (ceil(min(run_dtb)/cno.distBin):floor(max(run_dtb)/cno.distBin))*cno.distBin;

h1 = hist(reo_dtb(reo_numHS >= cno.minHS), distx);
h2 = hist(run_dtb, distx);
ad.reorate_disttobound = (h1./h2) * minpp;
ad.reorate_disttobound_eb = sqrt(h1)./h2* minpp;
ad.distx = distx;

h1 = hist(reo_dtb(abs(reo_ttb ) < pi/3 & reo_numHS >= cno.minHS), distx);
h2 = hist(run_dtb(abs(run_ttb ) < pi/3), distx);
ad.reorate_disttobound_tolight = h1./h2 * minpp;
ad.reorate_disttobound_tolight_eb = sqrt(h1)./h2* minpp;
ad.distx = distx;

h1 = hist(reo_dtb(abs(reo_ttb ) > 2*pi/3 & reo_numHS  >= cno.minHS), distx);
h2 = hist(run_dtb(abs(run_ttb) > 2*pi/3), distx);
ad.reorate_disttobound_todark = h1./h2 * minpp;
ad.reorate_disttobound_todark_eb = sqrt(h1)./h2* minpp;
ad.distx = distx;

h1 = hist(reo_dtb(abs(reo_ttb ) < pi/3 & reo_numHS  == 0), distx);
h2 = hist(run_dtb(abs(run_ttb ) < pi/3), distx);
ad.pauserate_disttobound_tolight = h1./h2 * minpp;
ad.pauserate_disttobound_tolight_eb = sqrt(h1)./h2* minpp;
ad.distx = distx;

h1 = hist(reo_dtb(abs(reo_ttb ) > 2*pi/3 & reo_numHS  == 0), distx);
h2 = hist(run_dtb(abs(run_ttb) > 2*pi/3), distx);
ad.pauserate_disttobound_todark = h1./h2 * minpp;
ad.pauserate_disttobound_todark_eb = sqrt(h1)./h2* minpp;
ad.distx = distx;

[~,ad.reosize_vs_ttb,ad.reosize_vs_ttb_eb] = meanyvsx(reo_ttb(reo_onb), abs(rad2deg(reo_dtheta(reo_onb))), deg2rad(txe));
ad.reosize_vs_ttb = ad.reosize_vs_ttb([1:end 1]);
ad.reosize_vs_ttb_eb = ad.reosize_vs_ttb_eb([1:end 1]);

[~,ad.reodir_vs_ttb,ad.reodir_vs_ttb_eb] = meanyvsx(reo_ttb(reo_onb), (rad2deg(reo_dtheta(reo_onb))), deg2rad(txe));
ad.reodir_vs_ttb = ad.reodir_vs_ttb([1:end 1]);
ad.reodir_vs_ttb_eb = ad.reodir_vs_ttb_eb([1:end 1]);

[~,ad.numhs_vs_ttb,ad.numhs_vs_ttb_eb] = meanyvsx(reo_ttb(reo_onb), reo_numHS(reo_onb), deg2rad(txe));
ad.numhs_vs_ttb = ad.numhs_vs_ttb([1:end 1]);
ad.numhs_vs_ttb_eb = ad.numhs_vs_ttb_eb([1:end 1]);


%headsweeps
%hs = eset.gatherField('headSwing');
hstx45 = -135:90:135;
hstxe45 = -180:90:180;
hstx0 = -180:90:90;
hstxe0 = -225:90:135;

ad.hstx0 = hstx0;
ad.hstx45 = hstx45;

hsvalidinds = hs_valid & hs_onb;

hs_ttb = adjustForPolarHistogram (hs_ttb, deg2rad(hstx45)); 
[~,ad.hs_acctoleft45, ad.hs_acctoleft45_eb] = meanyvsx(hs_ttb(hsvalidinds & (hs_dir > 0)), hs_accepted(hsvalidinds & (hs_dir > 0)),deg2rad(hstxe45));
[~,ad.hs_acctoright45, ad.hs_acctoright45_eb] = meanyvsx(hs_ttb(hsvalidinds & (hs_dir < 0)), hs_accepted(hsvalidinds & (hs_dir < 0)),deg2rad(hstxe45));
[~,ad.hs_acctolight45, ad.hs_acctolight45_eb] = meanyvsx(hs_ttb(hsvalidinds & (hs_delta_dist > 0)), hs_accepted(hsvalidinds & (hs_delta_dist > 0)),deg2rad(hstxe45));
[~,ad.hs_acctodark45, ad.hs_acctodark45_eb] = meanyvsx(hs_ttb(hsvalidinds & (hs_delta_dist < 0)), hs_accepted(hsvalidinds & (hs_delta_dist < 0)),deg2rad(hstxe45));

[~,ad.hs_dir45, ad.hs_dir45_eb] = meanyvsx(hs_ttb(hsvalidinds ), (hs_dir(hsvalidinds) > 0),deg2rad(hstxe45));

hs_ttb = adjustForPolarHistogram (hs_ttb, deg2rad(hstx0)); 
[~,ad.hs_acctoleft0, ad.hs_acctoleft0_eb] = meanyvsx(hs_ttb(hsvalidinds & (hs_dir > 0)), hs_accepted(hsvalidinds & (hs_dir > 0)),deg2rad(hstxe0));
[~,ad.hs_acctoright0, ad.hs_acctoright0_eb] = meanyvsx(hs_ttb(hsvalidinds & (hs_dir < 0)), hs_accepted(hsvalidinds & (hs_dir < 0)),deg2rad(hstxe0));
[~,ad.hs_dir0, ad.hs_dir0_eb] = meanyvsx(hs_ttb(hsvalidinds ), (hs_dir(hsvalidinds) > 0),deg2rad(hstxe0));
[~,ad.hs_acctolight0, ad.hs_acctolight0_eb] = meanyvsx(hs_ttb(hsvalidinds & (hs_delta_dist > 0)), hs_accepted(hsvalidinds & (hs_delta_dist > 0)),deg2rad(hstxe0));
[~,ad.hs_acctodark0, ad.hs_acctodark0_eb] = meanyvsx(hs_ttb(hsvalidinds & (hs_delta_dist < 0)), hs_accepted(hsvalidinds & (hs_delta_dist < 0)),deg2rad(hstxe0));

firsthsvalidinds = firsths_valid & firsths_onb;


firsths_ttb = adjustForPolarHistogram (firsths_ttb, deg2rad(hstx45)); 
[~,ad.firsths_acctoleft45, ad.firsths_acctoleft45_eb] = meanyvsx(firsths_ttb(firsthsvalidinds & (firsths_dir > 0)), firsths_accepted(firsthsvalidinds & (firsths_dir > 0)),deg2rad(hstxe45));
[~,ad.firsths_acctoright45, ad.firsths_acctoright45_eb] = meanyvsx(firsths_ttb(firsthsvalidinds & (firsths_dir < 0)), firsths_accepted(firsthsvalidinds & (firsths_dir < 0)),deg2rad(hstxe45));
[~,ad.firsths_dir45,ad.firsths_dir45_eb] = meanyvsx(firsths_ttb(firsthsvalidinds ), (firsths_dir(firsthsvalidinds) > 0),deg2rad(hstxe45));

firsths_ttb = adjustForPolarHistogram (firsths_ttb, deg2rad(hstx0)); 
[~,ad.firsths_acctoleft0, ad.firsths_acctoleft0_eb] = meanyvsx(firsths_ttb(firsthsvalidinds & (firsths_dir > 0)), firsths_accepted(firsthsvalidinds & (firsths_dir > 0)),deg2rad(hstxe0));
[~,ad.firsths_acctoright0, ad.firsths_acctoright0_eb] = meanyvsx(firsths_ttb(firsthsvalidinds & (firsths_dir < 0)), firsths_accepted(firsthsvalidinds & (firsths_dir < 0)),deg2rad(hstxe0));
[~,ad.firsths_dir0,ad.firsths_dir0_eb] = meanyvsx(firsths_ttb(firsthsvalidinds ), (firsths_dir(firsthsvalidinds) > 0),deg2rad(hstxe0));
