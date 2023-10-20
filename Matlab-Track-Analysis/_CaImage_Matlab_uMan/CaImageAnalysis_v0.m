% CaImageAnalysis_v0.m
% converts uManTime to stimTime

%% Prior to this script:
% Quantify ROIs from video file in imageJ
% 1. With first ROI as background, Create ROIs in ROI manager for all neurons
%    - analyze-> tools-> ROI manager
%    - save ROIs within ROI manager, more-> save
% 2. Set measurement conditions 
%   - analyze-> set measurements-> mean gray value
% 3. Meaure with ROI measure, - more-> multimeasure
%   - measure all slices & one row per slice should be checked
% 4. Set save conditions within measurement table
%   - Results-> Options, uncheck all Results Table Options values
% 5. Save as default name, 'Results.csv', within same directory as video.



%% get video file & align to stimulus information
[imFile,imDir]=uigetfile('*.tif', 'Select video file for analysis');

[ tsF, tsT ] = produceTimeSeries(imFile,imDir);






% dFtoF with dT/dt & dF/dt on time axis
for i=1:size(tsF.data,2)
    smDiff=diff(tsT.data);
    threshDiff= 10*median(abs(smDiff));
    figure();hold on; plot(tsF.data(:,i)); 
    yyaxis right; plot(smDiff); plot(diff(tsF.data(:,i)),'color',[0 0 1]);
    line(1:size(tsT.data,1),repmat(threshDiff,[size(tsT.data,1),1]));
    text(.5,.9,strcat('worm #',num2str(i)),'units', 'normalized');
    posInfo=get(gcf,'position');
    set(gcf,'position',[posInfo(1)/4, posInfo(2:4)]);
end

% for i=1:size(tsF.data,2)
%     figure();plot(tsT.data,tsF.data(:,i),'linestyle', 'none','marker','*','color',rand(1,3))
%     text(.5,.9,strcat('worm #',num2str(i)),'units', 'normalized');
%     posInfo=get(gcf,'position');
%     set(gcf,'position',[posInfo(1)/4, posInfo(2:4)]);
% end


for i=1:size(tsF.data,2)
    figure(); plot(diff(tsT.data),diff(tsF.data(:,i)),'linestyle', 'none','marker','*','color',rand(1,3));
    text(.5,.9,strcat('worm #',num2str(i)),'units', 'normalized');
    posInfo=get(gcf,'position');
    set(gcf,'position',[(posInfo(1)+3*posInfo(1)/4), posInfo(2:4)]);   
end

