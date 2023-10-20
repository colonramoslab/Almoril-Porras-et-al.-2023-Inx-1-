basedirs = {'E:\larvalco2\Extracted\50 mL CO2 in 2 L air\CS', 'E:\larvalco2\Extracted\50 mL CO2 in 2 L air\gr63a',...
    'E:\larvalco2\Extracted\10mL CO2 2L air\CS', 'E:\larvalco2\Extracted\5p4 mL co2 in 2 L air\CS',...
    'E:\larvalco2\Extracted\control air 50 mL in 2L air\CS', 'E:\larvalco2\Extracted\no air flow\CS',...
    'E:\larvalco2\Extracted\20p4mL CO2 2L air\CS'};
esetnames = {'cs_50', 'gr63a_50', 'cs_10', 'cs_5', 'cs_0', 'cs_noair', 'cs_20'};
labelnames = {'CS 0-5\%', 'Gr63a 0-5\%', 'CS 0-1\%', 'CS 0-0.5\%', 'CS air control', 'CS no air control', 'CS 0-2\%'};
figdir = 'C:\Users\Marc\Documents\figures\lady gaga paper\co2 spatial';

for j = 1:length(basedirs)
    d = dir(fullfile (basedirs{j}, 'matfiles', [esetnames{j} '_experiment*.mat']));
    if (isempty(d))
        error (['problem with ' num2str(j)']);
    end
end

descriptions = {'CS CO$_2$ 2500 ppm/cm', 'Gr63a CO$_2$ 2500 ppm/cm', 'CS CO$_2$ 500 ppm/cm', 'CS CO$_2$ 250 ppm/cm', ...
    'CS clean air', 'CS no air flow', 'CS CO$_2$ 1000 ppm/cm'};

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
if (~exist('co2_spatial_statistics'))
    co2_spatial_statistics = calculateStatisticsOfEset(eset, 'timerange', timerange);
    save (fullfile(figdir, 'co2 spatial statistics'), 'co2_spatial_statistics');

end

es = co2_spatial_statistics;
report = {'co2 spatial statistics'};
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
    
    ln = ln + 1;report{ln} =  ([labelnames{j} ': ' num2str(co2_spatial_statistics(j).numExpts) ' expts, ' num2str(co2_spatial_statistics(j).numAnimals) ' animals, ' num2str(co2_spatial_statistics(j).animalTime/3600, 2) ' hrs of observation']);
end
ln = ln + 1; report{ln} = ['---- experiment: ' labelnames{1} ' --------'];
ln = ln + 1;report{ln} = ('FIGURE 3');
ln = ln + 1;report{ln} = (['3c:  runs up co2 gradient to ' num2str(es(1).directions(1)) ': ' num2str(es(1).numRunsInDirection(1))]);
ln = ln + 1;report{ln} = (['3c:  runs down co2 gradient to ' num2str(es(1).directions(2)) ': ' num2str(es(1).numRunsInDirection(2))]);
for j = 1:4
    ln = ln + 1;report{ln} = (['3e: num reos from ' num2str(es(1).directions(j)) ': ' num2str(es(1).numReosWithHSFromDirection(j))]);
end

ln = ln + 1;report{ln} =('FIGURE 4');
ln = ln + 1;report{ln} = ['4a,b,c) co2 total reorientations = ' num2str(sum(es(1).numReosWithHSFromDirection)) ' should = ' num2str(nnz(es(1).reoNumHS > 0))];
ln = ln + 1;report{ln} = ['4d: co2 num reorientations perp to gradient = ' num2str(sum(es(1).numReosWithHSFromDirection(3:4)))];
ln = ln + 1;report{ln} = ['4e: co2 num hs perp to gradient = ' num2str(sum(es(1).numHSFromDirection(3:4)))];
ln = ln + 1;report{ln} = ['4d: co2 num first hs perp to gradient = ' num2str(nnz((abs(cos(es(1).firstHeadSwingPrevDir)) < 1/sqrt(2)) & es(1).firstHeadSwingValid))];

disp(report')

fid = fopen (fullfile(figdir, 'co2 spatial nums for fig.txt'),'w');
for j = 1:length(report)
    fprintf(fid, '%s\n', report{j});
end
fclose(fid);
clear fid;