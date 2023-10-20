function [tsF, tsT] = CaImageAnalysis(varargin)
% CaImageAnalysis:
% converts uManTime to stimTime

% Input: pointer to image file for analysis.
%       Empty=> gui select image file. OK
%       file name as string. OK
%       Folder pick image file (first). NOTE ISSUE: HAVE TO PICK stimFile!


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

callVis=0; % set to 1 to see all traces individually with dT/dt
imFile=[];
imDir=[];

if length(varargin)==1
    [imDir, imFile, ext]=fileparts(varargin{1});
    if isempty(ext)
        imDir=varargin{1};
        retDir=pwd;
        cd(imDir);
        imFs=dir('*.tif');
        imFile=imFs(1).name;
        cd(retDir);
    end
    imFile=strcat(imFile,ext);
end

varargin=assignApplicable(varargin);


if isempty(imFile)
    %% get video file & align to stimulus information
    [imFile,imDir]=uigetfile('*.tif', 'Select video file for analysis');
end

if ~strcmp('\',imDir(end))
    imDir=strcat(imDir,'\');
end

[ tsF, tsT ] = produceTimeSeries(imFile,imDir);


if callVis==1
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
end

end

