basedirs = {'E:\larvalco2\Extracted\Ethyl Acetate 233ppm CS\50mL air in 1.5 L air\CS', 'E:\larvalco2\Extracted\Ethyl Acetate 2.25ppm CS\40mL air in 2 L air\CS', ...
    'E:\larvalco2\Extracted\Ethyl Acetate masking experiment\25ppm masking 25ppm in vlaves\CS', 'E:\larvalco2\Extracted\Ethyl Acetate 331ppm Gr63a\50mL air in 1 L air\Gr63a',...
    'E:\larvalco2\Extracted\Ethyl Acetate pure Spatial Gr63a\50mL co2 in 2 L air\Gr63a', 'E:\larvalco2\Extracted\ethyl acetate 2000X diluted Spatial\50mL air through 2000X diluted ethyl acetate\CS',...
    'E:\larvalco2\Extracted\Ethyl Acetate Or83b1\23ppm at output spatial\Or83b1', 'E:\larvalco2\Extracted\ethyl acetate Or83b2 spatial\23ppm middle in 2L air\Or83b2'}; 

esetnames = {'cs_ea_0_460', 'cs_ea_0_5', 'cs_ea_25_50', 'gr63a_ea_0_660', 'gr63a_ea_50_2', 'cs_ea_0_48', 'or83b_ea_0_46','Or83b2_ea_0_46'};
figdir = 'C:\Users\Marc\Documents\figures\lady gaga paper\ethyl acetate spatial\';
for j = 1:length(basedirs)
    d = dir(fullfile (basedirs{j}, 'matfiles', [esetnames{j} '_experiment*.mat']));
    if (isempty(d))
        error (['problem with ' num2str(j)']);
    end
end

%undiluted = 2200 ppm at 50 / 2000, which would say vapor pressure of ethyl
%acetate is 0.088 atm = 
%<math>\scriptstyle \log_{10} P_{mmHg} = 7.09808 - \frac {1238.71}
%{217.0+T}</math> = 82 mm mercury at 22 C = 0.108 atm -- close

%2000 x dilution at 50 / 2000 approx 20 ppm mean

descriptions = {'CS EA 20 ppm/cm', 'CS EA 0.2 ppm/cm', 'CS EA 2 ppm/cm, masked by 25 ppm EA', 'GR63a EA 30 ppm/cm', ...
    'GR63a EA 200 ppm/cm', 'CS EA 2 ppm/cm', 'OR83b1 EA 2 ppm/cm', 'OR83b2 EA 2 ppm/cm'};

toload = [6 1 2 7 8];

basedirs = basedirs(toload);
esetnames = esetnames(toload);
descriptions = descriptions(toload);

if (~exist('eset', 'var'))
    ecl = ESetCleaner();
    ecl.rpmCut = 1;
    ecl.askFirst = false;
    ecl.showFigsInReport = false;
    for j = 1:length(basedirs)
        eset(j) = ExperimentSet.fromMatFiles(fullfile (basedirs{j}, 'matfiles', esetnames{j}));
        ecl.getReport(eset(j));
        ecl.clean(eset(j));
    end
    resegment = true;
end
if (exist('eset', 'var') && resegment)
     for j = 1:length(eset)
        eset(j).executeTrackFunction('setSegmentSpeeds');
        eset(j).executeTrackFunction('segmentTrack');
     end
     resegment = false;
end

timerange = [120 1020];
if (~exist('etac_spatial_statistics', 'var'))
    etac_spatial_statistics = calculateStatisticsOfEset(eset, 'timerange', timerange);
    save (fullfile(figdir, 'etac spatial statistics'), 'etac_spatial_statistics');

end

es = etac_spatial_statistics;
report = {'ethyl acetate spatial statistics'};
ln = 1;
ln = ln + 1;report{ln} = [num2str(mean([es.headSwingAccepted])) ' of headswings are accepted.'];
ln = ln + 1;report{ln} = [num2str(mean([es.headSwingValid])) ' of headswings are valid.'];
hsa = [es.headSwingAccepted];
ln = ln + 1;report{ln} = [num2str(mean(hsa(logical([es.headSwingValid])))) ' of valid headswings are accepted.'];
ln = ln + 1;report{ln} = [num2str(mean([es.lastHeadSwingAccepted])) ' of last headswings are accepted.'];
ln = ln + 1;report{ln} = [num2str(mean([es.lastHeadSwingValid])) ' of last headswings are valid.'];
hsa = [es.lastHeadSwingAccepted];
ln = ln + 1;report{ln} = [num2str(mean(hsa((logical([es.lastHeadSwingValid]))))) ' of valid last headswings are accepted.'];

ln = ln + 1;report{ln} = 'FIGURE 2';


for j = 1:length(es)
    
    ln = ln + 1;report{ln} =  ([descriptions{j} ': ' num2str(etac_spatial_statistics(j).numExpts) ' expts, ' num2str(etac_spatial_statistics(j).numAnimals) ' animals, ' num2str(etac_spatial_statistics(j).animalTime/3600, 2) ' hrs of observation']);
end
ln = ln + 1; report{ln} = ['---- experiment: ' descriptions{1} ' --------'];
ln = ln + 1;report{ln} = ('FIGURE 3');
ln = ln + 1;report{ln} = (['3c:  runs up etac gradient to ' num2str(es(1).directions(1)) ': ' num2str(es(1).numRunsInDirection(1))]);
ln = ln + 1;report{ln} = (['3c:  runs down etac gradient to ' num2str(es(1).directions(2)) ': ' num2str(es(1).numRunsInDirection(2))]);
for j = 1:4
    ln = ln + 1;report{ln} = (['3e: num reos with 1 or more HS from ' num2str(es(1).directions(j)) ': ' num2str(es(1).numReosWithHSFromDirection(j))]);
end

ln = ln + 1;report{ln} =('FIGURE 4');
ln = ln + 1;report{ln} = ['4a,b,c) etac total reorientations = ' num2str(sum(es(1).numReosWithHSFromDirection)) ' should = ' num2str(nnz(es(1).reoNumHS > 0 & es(1).reoStartTime >= timerange(1) & es(1).reoStartTime <= timerange(2)))];
ln = ln + 1;report{ln} = ['4d: etac num reorientations perp to gradient = ' num2str(sum(es(1).numReosWithHSFromDirection(3:4)))];
ln = ln + 1;report{ln} = ['4e: etac num valid hs perp to gradient = ' num2str(sum(es(1).numHSFromDirection(3:4)))];
ln = ln + 1;report{ln} = ['4d: etac num valid first hs perp to gradient = ' num2str(nnz((abs(cos(es(1).firstHeadSwingPrevDir)) < 1/sqrt(2)) & es(1).firstHeadSwingValid & es(1).firstHeadSwingTime >= timerange(1) & es(1).firstHeadSwingTime <= timerange(2)))];
ln = ln + 1;report{ln} = ['4d: etac prob valid first hs perp to gradient accepted = ' num2str(mean(es(1).firstHeadSwingAccepted((abs(cos(es(1).firstHeadSwingPrevDir)) < 1/sqrt(2)) & es(1).firstHeadSwingValid & es(1).firstHeadSwingTime >= timerange(1) & es(1).firstHeadSwingTime <= timerange(2))))];

disp(report')

fid = fopen (fullfile(figdir, 'etac spatial nums for fig.txt'),'w');
for j = 1:length(report)
    fprintf(fid, '%s\n', report{j});
end
fclose(fid);
clear fid;