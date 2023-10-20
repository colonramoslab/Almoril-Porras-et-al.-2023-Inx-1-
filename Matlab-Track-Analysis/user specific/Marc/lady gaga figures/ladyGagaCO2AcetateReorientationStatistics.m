%statistics for figure 3(d)
figdir = 'C:\Users\Marc\Documents\figures\lady gaga paper\co2 spatial';

if (~exist('co2_spatial_ad', 'var'))
    load (fullfile(figdir, 'co2_spatial_calcs'));
end
descriptions = {'CS CO$_2$ 2500 ppm/cm', 'Gr63a CO$_2$ 2500 ppm/cm', 'CS CO$_2$ 500 ppm/cm', 'CS CO$_2$ 250 ppm/cm', ...
    'CS clean air', 'CS no air flow', 'CS CO$_2$ 1000 ppm/cm'};


ad = co2_spatial_ad(1);
welcht = (ad.headSwingAcceptanceRateTowards - ad.headSwingAcceptanceRateAway)/sqrt(ad.headSwingAcceptanceRateTowards_eb^2 + ad.headSwingAcceptanceRateAway_eb.^2);

dir = [0, 180, 90, -90];
pleft_data = zeros(size(dir));
pleft_data_eb = pleft_data;
rsize_data = zeros(size(dir));
rsize_data_eb = rsize_data;
rdir_data = rsize_data;
rdir_data_eb = rdir_data;
for j = 1:4
    inds = cos(ad.reo_prevdir - deg2rad(dir(j))) > 1/sqrt(2) & ad.reo_numHS > 0;
    nreo(j) = nnz(inds);
    pleft_data(j) = mean(ad.reo_dtheta(inds) > 0);
    pleft_data_eb(j) = sqrt(pleft_data(j)*(1 - pleft_data(j))/nnz(inds));
    
    rsize_data(j) = rad2deg(sqrt(mean(ad.reo_dtheta(inds).^2)));
    rsm = mean(ad.reo_dtheta(inds).^2);
    rsv = std(ad.reo_dtheta(inds).^2);
    rsize_data_eb(j) = rad2deg(std(mean(rsm + rsv*randn(nnz(inds),1000),1)));
    
    rdir_data(j) = rad2deg(mean(ad.reo_dtheta(inds)));
    rdir_data_eb(j) = rad2deg(std(ad.reo_dtheta(inds)))/nnz(inds);
    
    [rdm, rsm, rpl] = ad.reo_dir_model.reorientationMeansVsTheta(ad.reo_dir_model.params, ad.reo_prevdir(inds), ad.reo_dir_model.hessian, 0.95);
    rdir_model_exact(j) = rad2deg(mean(rdm));
    rsize_model_exact(j) = rad2deg(sqrt(mean(rsm)));
    pleft_model_exact(j) = mean(rpl);
    
    fun = @(x) mean(ad.reo_dir_model.reorientationMeansVsTheta(x, ad.reo_prevdir(inds), ad.reo_dir_model.hessian, 0.95));
    minrd = Inf;
    maxrd = -Inf;
    for k = 1:10
        [mnr, mxr] = confidenceRange(ad.reo_dir_model.params, ad.reo_dir_model.hessian, fun, 0.95, 'numtrials', 1E3);
        minrd = min(minrd,mnr);
        maxrd = max(maxrd,mxr);
    end
    rdir_model_exact_ci(j,:) = rad2deg([minrd, maxrd]);
    
    fun = @(x) mean(getSomeArgsOut(3, ad.reo_dir_model.reorientationMeansVsTheta, x, ad.reo_prevdir(inds), ad.reo_dir_model.hessian, 0.95));
    minpl = Inf;
    maxpl = -Inf;
    for k = 1:10
        [mnp, maxp] = confidenceRange(ad.reo_dir_model.params, ad.reo_dir_model.hessian, fun, 0.95, 'numtrials', 1E3);
        minpl = min(minpl,mnp);
        maxpl = max(maxpl,maxp);
    end
    pleft_model_exact_ci(j,:) = [minpl, maxpl];
    
    minrs = Inf;
    maxrs = -Inf;
    for k = 1:10
        fun = @(x) sqrt(mean(getSomeArgsOut(2, ad.reo_dir_model.reorientationMeansVsTheta, x, ad.reo_prevdir(inds), ad.reo_dir_model.hessian, 0.95)));
        [mns, mxs] = confidenceRange(ad.reo_dir_model.params, ad.reo_dir_model.hessian, fun, 0.95, 'numtrials', 1E3);
        minrs = min(minrs,mns);
        maxrs = max(maxrs,mxs);
    end
    rsize_model_exact_ci(j,:) = rad2deg([minrs, maxrs]);
    
    
    
end

[rdir_model, rsize_model, pleft_model] = ad.reo_dir_model.reorientationMeansVsTheta(ad.reo_dir_model.params, deg2rad(dir), ad.reo_dir_model.hessian, 0.95);
rdir_model = rad2deg(rdir_model);
rsize_model = rad2deg(sqrt(rsize_model));

clear report;
ln = 1; report{ln} = [descriptions{1} ' - reorientation statistics'];
for j = 1:4
    ln = ln+1; report{ln} = ['direction - ' num2str(dir(j))];
    ln = ln+1; report{ln} = ['number of reorientations: ' num2str(nreo(j))];
    ln = ln+1; report{ln} = ['p(left) = (data) ' num2str(100*pleft_data(j),3) ' +/- ' num2str(100*pleft_data_eb(j),3)];
    ln = ln+1; report{ln} = ['p(left) = (model) ' num2str(100*pleft_model_exact(j),3) '{ 95% conf region  ' num2str(100*pleft_model_exact_ci(j,1),3) ' - ' num2str(100*pleft_model_exact_ci(j,2),3) ' }'];
    ln = ln+1; report{ln} = ['meandir = (data) ' num2str(rdir_data(j),3) ' +/- ' num2str(rdir_data_eb(j),3)];
    ln = ln+1; report{ln} = ['meandir = (model) ' num2str(rdir_model_exact(j),3) '{ 95% conf region  ' num2str(rdir_model_exact_ci(j,1),3) ' - ' num2str(rdir_model_exact_ci(j,2),3) ' }'];
    ln = ln+1; report{ln} = ['rms dtheta = (data) ' num2str(rsize_data(j),3) ' +/- ' num2str(rsize_data_eb(j),3)];
    ln = ln+1; report{ln} = ['rms dtheta = (model) ' num2str(rsize_model_exact(j),3) '{ 95% conf region  ' num2str(rsize_model_exact_ci(j,1),3) ' - ' num2str(rsize_model_exact_ci(j,2),3) ' }'];
end

ss = {'condition', 'direction', 'num reorientations', 'p(left) meas', 'std. err p(left)', 'p(left) model', 'p(left) model 95% ci lb', 'p(left) model 95% ci ub', 'rms dtheta meas', 'std. err rms dtheta', 'rms dtheta model', 'rms dt ci lb', 'rms ct ci ub'};
for j = 1:4
    ss(j+1,:) = {'co2', dir(j), nreo(j), pleft_data(j)*100, pleft_data_eb(j)*100, pleft_model_exact(j)*100, pleft_model_exact_ci(j,1)*100, pleft_model_exact_ci(j,2)*100, rsize_data(j), rsize_data_eb(j), rsize_model_exact(j), rsize_model_exact_ci(j,1), rsize_model_exact_ci(j,2)};
end

%report'
xlswrite(fullfile(figdir, 'etac spatial reo stats for fig 3d.xls'), ss);



fid = fopen (fullfile(figdir, 'co2 spatial reo nums for fig 3.txt'),'w');
for j = 1:length(report)
    fprintf(fid, '%s\n', report{j});
end
fclose(fid);
clear fid;
    

