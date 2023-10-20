basedirs = {'E:\larvalco2\Extracted\Ethyl Acetate Temporal\sq 4min 0-5p2mL no cond\CS', 'E:\larvalco2\Extracted\Ethyl Acetate Temporal\sq 4min 0-10mL no cond\CS', 'E:\larvalco2\Extracted\Ethyl Acetate Temporal\sq 4min 0-50mL no cond\CS'};
    
esetnames = {'cs_ea_0_220_sq_240','cs_ea_0_440_sq_240' ,'cs_ea_0_2200_sq_240'};

if (~exist('eset', 'var'))
    for j = 1:3
        eset(j) = ExperimentSet.fromMatFiles(fullfile(basedirs{j}, 'matfiles', esetnames{j}));
        ecl = ESetCleaner;
        ecl.rpmCut = 1;
        ecl.askFirst = false;
        ecl.showFigsInReport = false;
        ecl.getReport(eset(j));
        ecl.clean(eset(j));
    end

    for j = 1:3
        eset(j).executeTrackFunction('segmentTrack');
    end
end

tno = temporalNavigationAnalysis;
tno.fieldname = 'vocppm';
tno.rampType = 'square';
if (~exist('tc', 'var'))
    tc = temporalCalculations(eset, tno);
end
tno.timeBinSize = 20;
tno.period = 240;
if (~exist('ad', 'var'))
    ad = temporalNavigationAnalysis(eset, 'vocppm', tno, tc);
end