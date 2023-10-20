setEthylAcetateFilesTemporal

existsAndDefault('timerange',  []);
existsAndDefault('trimrect', []);
existsAndDefault ('redogqs', false);
existsAndDefault('retrim', false(size(basedirs)));
existsAndDefault('resegment', false(size(basedirs)));

if (~exist('ea_eset', 'var'))
   % if (matlabpool('size') == 0)
   %     matlabpool
   % end
   
    for j = 1:length(basedirs)
        ea_eset(j) = ExperimentSet.fromMatFiles(fullfile(basedirs{j}, 'matfiles',esetnames{j}),'all',false); %#ok<SAGROW>
        resegment(j) = true;
       % retrim(j) = true;
    end
    %matlabpool close
end
%
%labelnames = {'Ethyl Acetate (conditioned) 0-1100 ppm, 220 ppm/min', 'Ethyl Acetate (unconditioned) 0-2200 ppm, 440 ppm/min', 'Ethyl Acetate (conditioned) 0-2200 ppm, 440 ppm/min',...
%    'Ethyl Acetate (conditioned) 1100-2200 ppm, 220 ppm/min', 'Ethyl Acetate (unconditioned) 0-2200 ppm exponential', 'Ethyl Acetate (unconditioned) 0-22 ppm exponential', 'Ethyl Acetate (unconditioned) 0-0.2 ppm exponential',...
%    'Ethyl Acetate square 0-220 ppm',...
%    'Ethyl Acetate square 0-440 ppm', 'Ethyl Acetate square 0-2200 ppm','Ethyl Acetate Exponential 0 - 2 ppm'};
%%
ecl = ESetCleaner;
ecl.rpmCut = 1;
ecl.showFigsInReport = false;
ecl.askFirst = false;
if (redogqs || isempty(ea_eset(end).expt(end).globalQuantity))
    %set the two exponential files by set point instead of meter
    
    
    for j = 1:length(ea_eset)
        ecl.getReport(ea_eset(j));
        ecl.clean(ea_eset(j));
        [ea_eset(j).expt.globalQuantity] = deal([]);
        ea_eset(j).executeExperimentFunction('addMetaDataFields');
    end
    
    esetnums = [7 11];
    vocampl = [0.2 2]/50;
    for j = 1:length(esetnums)
        for k = 1:length(ea_eset(esetnums(j)).expt)
            expt = ea_eset(esetnums(j)).expt(k);
            expt.addGlobalQuantity(vocppmFromSetPoint(expt.globalQuantity(1), vocampl(j)));
        end
    end
    
    good = {[3,4], [1,2,3], [1], [2,3], 2:4, 1:5, [], 1, 1, 1, []};
    bad = {[1,2], [], [2,3,4], [1], [1], [], [], 2, 2, 2};

    for j = 1:length(ea_eset);
        
           
        if (~isempty(good{j}) && ~isempty(bad{j}))
            expt = ea_eset(j).expt;
            gind = good{j}(1);
            for k = 1:length(bad{j})
                badind = bad{j}(k);
                g1 = expt(badind).globalQuantity(1).yData;
                g2 = expt(gind).globalQuantity(1).yData;
                g1 = g1(1:min(length(g1), length(g2)));
                g2 = g2(1:min(length(g1), length(g2)));
                [xc,lags] = xcorr(g1,g2);
                [~,I] = max(xc);
                u = (1:length(g2)) - lags(I);
                g2new = interp1(g2, u, 'linear', 'extrap');
                u = (1:length(expt(badind).globalQuantity(3).xData.et)) - lags(I);
                expt(badind).globalQuantity(3).yData = interp1(expt(gind).globalQuantity(3).yData, u, 'linear', 'extrap');
                expt(badind).assignGlobalQuantities();
            end
       end
       ea_eset(j).executeExperimentFunction('assignGlobalQuantities');
       gq = GlobalQuantity;
       gq.fieldname = 'dlogvocppm';
       gq.xField = {'vocppm', 'dvocppm'};
       gq.xData = []; gq.yData = [];
       gq.derivationMethod = @(xin, xdata, ydata) xin.dvocppm./xin.vocppm;
       for k = 1:length(ea_eset(j).expt)
           ea_eset(j).expt(k).addGlobalQuantity(gq);
       end
       ea_eset(j).expt.addTonToff('vocppm', ramptype{j});
       dqs = {'curv', 'sspineTheta', 'vel_dp', 'speed'};
       for k = 1:length(dqs)
           ea_eset(j).gatherField(dqs{k});
       end
       ea_eset(j).toMatFiles(fullfile(basedirs{j}, 'matfiles', esetnames{j}));
       
    end
    redogqs = false;
    
end


%{
for j = 1:length(ea_eset)
    expt = ea_eset(j).expt;
    for k = 1:length(expt)
        etg = expt(k).globalQuantity(1).xData.et;
        gassp = expt(k).globalQuantity(1).yData;
        etv = expt(k).globalQuantity(3).xData.et;
        vcoppm = expt(k).globalQuantity(3).yData;
        plot (etg, gassp, etv, vcoppm);
        title ([num2str(j) ' , ' num2str(k)]);
        pause
    end
end
  %}      

for j = 1:length(ea_eset)
    if (retrim(j))
        disp(['trimming - ' num2str(j)']);
        ea_eset(j).executeExperimentFunction('trimTracks', timerange, trimrect);
        resegment(j) = true; %#ok<SAGROW>
        retrim(j) = false;
        clear tc;
    end
    if resegment(j)
        disp(['segmenting - ' num2str(j)]);
        ea_eset(j).executeTrackFunction('setSegmentSpeeds');
        ea_eset(j).executeTrackFunction('segmentTrack');
        resegment(j) = false; %#ok<SAGROW>
        clear tc;
    end
    po(j).legendEntry = labelnames{j}; %#ok<SAGROW>
    po(j).onColor = [0.6 0 0.6];
    po(j).reverse = true;
    
end
return;
period = [600 600 600 600 600 600 600 240 240 240 600];
for j = 1:length(period)
    tno(j).period = period(j);
    tno(j).timeBinSize = tno(j).period/20;
    tno(j).rampType = ramptype{j};
    tno.fieldname = 'vocppm';
    tno.timerange = [180 1800];
end    

if (~exist ('tc', 'var'))
    for j = 1:length(tno)
        tc{j} = temporalCalculations(ea_eset(j), tno(j));
    end
    clear ad;
end
if (~exist ('ad', 'var'))
    for j = 1:length(tc)
        ad(j) = temporalNavigationAnalysis(ea_eset(j), tno(j).fieldname, tno(j), tc{j});
    end
end


temporalNavigationFigures(ad, po);
