%working on chi-squared for reorientation model for ea_spatial 2 ppm/cm

%load lady gaga spatial ethyl acetate calcs, set ad = ea_spatial_ad(6);
if (~exist('ad', 'var'))
    ladyGagaSpatialEthylAcetate
    close all
    ad = ea_spatial_ad(6);
end

pd = ad.reo_prevdir;
dt = ad.reo_dtheta;
nhs = ad.reo_numHS;

pd = pd(nhs > 0);
dt = dt(nhs > 0);

nreo = length(pd);
model = ad.reo_dir_model;
clear x2;

nbins = 180;
bs = 180/nbins;
bincenters = linspace(bs/2, 180-bs/2, nbins);

snp = @(x, xd) max(skewNormalPDF(xd, x(1), x(2), x(3))+skewNormalPDF(-xd, x(1), x(2), x(3))+skewNormalPDF(2*pi-xd, x(1), x(2), x(3)), eps);
np = @(x,xd) max(normpdf(xd, x(1), x(2)) + normpdf(-xd, x(1), x(2))+ normpdf(2*pi-xd, x(1), x(2)), eps);
skewnormfun = @(x) -sum(log(snp(x,abs(dt))));
normfun = @(x) -sum(log(np(x,abs(dt))));
x = [mean(abs(dt))/2, std(abs(dt)), 2];

skewnormparams = fminsearch(skewnormfun, x);

%this should be unnecessary because taking the mean and standard deviation
%already maximizes the likelihood of the data for a gaussian, but since
%we're adding in the reflection about 0, we'll do it anyway
x = [mean(abs(dt)), std(abs(dt))];
normparams = fminsearch(normfun, x);
%{
nbins = 30:5:150;
for j = 1:length(nbins)
    bs = 180/nbins(j);
    bincenters = linspace(bs/2, 180-bs/2, nbins(j));

   
    chisq_sn_red(j) = sum((o-e_sn).^2./e_sn)/(nbins(j) - 3);
    chisq_n_red(j) = sum((o-e_n).^2./e_sn)/(nbins(j) - 2);
end
plot (nbins, chisq_sn_red, nbins, chisq_n_red);
return;
%}

o = hist(abs(dt), deg2rad(bincenters));
e_sn = snp(skewnormparams, deg2rad(bincenters))*deg2rad(bs)*length(dt);
e_n = np(normparams, deg2rad(bincenters))*deg2rad(bs)*length(dt);
figure(1); clf();
bar(bincenters, o); hold on;

myfun = @(x) sum((snp(x, deg2rad(bincenters))*deg2rad(bs)*length(dt)-o).^2./(snp(x, deg2rad(bincenters))*deg2rad(bs)*length(dt)));
skewnormparams_mincs = fminsearch(myfun, skewnormparams);
e_sn_mincs = snp(skewnormparams_mincs, deg2rad(bincenters))*deg2rad(bs)*length(dt);


plot (bincenters, e_sn, 'r--', bincenters, e_sn_mincs, 'k--', bincenters, e_n , 'g--', 'LineWidth', 3);
chisq_sn_red = sum((o-e_sn).^2./e_sn)/(nbins - 3)
chisq_n_red = sum((o-e_n).^2./e_sn)/(nbins - 2)
chisq_sn_mincs_red = sum((o-e_sn_mincs).^2./e_sn_mincs)/(nbins - 3)


nbins = 30; 
bs = 360/nbins;
bincenters = linspace(-180+bs/2, 180-bs/2, nbins);
if (~exist('pdf', 'var'))
    pdf = zeros(length(bincenters), length(pd));
    
    for k = 1:length(pd)
        pdf(:,k) = model.pdfOfdTheta(model.params, pd(k), deg2rad(bincenters))*deg2rad(bs);
    end
    extra = (1-sum(pdf));
    ratio = sum(pdf(bincenters < 0,:)) + 0.5*sum(pdf(bincenters == 0,:))./sum(pdf);
%    pdf(1,:) = pdf(1,:) + extra .* ratio;
 %   pdf(end,:) = pdf(end,:) + extra.*(1 - ratio);
end

directions = [0 180 90 -90];
o = zeros(length(bincenters), length(directions));
e = o;
for j = 1:length(directions)
    inds = cos(pd - deg2rad(directions(j))) > 1/sqrt(2);
    o(:,j)= hist(dt(inds), deg2rad(bincenters));
    e(:,j) = sum(pdf(:,inds),2);
end
figure(2); clf();
for j = 1:4
    subplot(2,2,j); bar(bincenters, o(:,j)); hold on; plot (bincenters, e(:,j), 'r-'); hold off;
    title (num2str(directions(j)));
end
x2 = sum((o(:) - e(:)).^2./(e(:)));
x2reduced = x2 / (length(o(:)) - 7)

nixedParams = {{},{'C'}, {'B', 'E'}, {'D', 'E'}, {'B','C','E'}};
altmodeldescriptions = {'orig model', 'no left/right bias', 'no size bias', 'no skew', 'no bias'};
if (~exist('altmodel', 'var'))
    for j = 1:length(nixedParams)
        clear startValues fixedValues
        for k = 1:length(model.params)
            startValues.(model.paramkey{k}) = model.params(k);
        end
        fixedValues = [];
        for k = 1:length(nixedParams{j})
            startValues.(nixedParams{j}{k}) = 0;
            fixedValues.(nixedParams{j}{k}) = 0;
        end
        altmodel(j) = fitReorientationAngleDistribution(pd, dt, fixedValues, startValues);
    end
end
if (~exist('altpdf', 'var'))
    for j = 1:length(altmodel)
        altpdf{j} = zeros(length(bincenters), length(pd));
        
        for k = 1:length(pd)
            altpdf{j}(:,k) = altmodel(j).pdfOfdTheta(altmodel(j).params, pd(k), deg2rad(bincenters))*deg2rad(bs);
        end
        extra = (1-sum(altpdf{j}));
        ratio = sum(altpdf{j}(bincenters < 0,:)) + 0.5*sum(altpdf{j}(bincenters == 0,:))./sum(altpdf{j});
     %   altpdf{j}(1,:) = altpdf{j}(1,:) + extra .* ratio;
      %  altpdf{j}(end,:) = altpdf{j}(end,:) + extra.*(1 - ratio);
    end
end

for j = 1:length(directions)
    for k = 1:length(altmodel)
        inds = cos(pd - deg2rad(directions(j))) > 1/sqrt(2);
        alte{k}(:,j) = sum(altpdf{k}(:,inds),2);
    end
end

ls = {'c-','r--', 'g--', 'm--', 'y--'};

figure(2); clf();
for j = 1:4
    subplot(2,2,j); bar(bincenters, o(:,j)); hold on; plot (bincenters, e(:,j), 'k-', 'LineWidth', 2); 
    for k = 1:length(altmodel)
         plot (bincenters, alte{k}(:,j), ls{k}, 'LineWidth', 2); 
    end
    hold off;
    title (num2str(directions(j)));
end
clear altx2;
for k = 1:length(altmodel)
    altx2(k) = sum((o(:) - alte{k}(:)).^2./(e(:)));
    altx2reduced(k) = altx2(k) / (length(o(:)) - 7 + length(nixedParams{k}));
end
