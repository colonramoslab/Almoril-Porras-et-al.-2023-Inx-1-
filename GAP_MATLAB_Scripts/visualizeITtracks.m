function [fPath, hRat, ITscore] = visualizeITtracks(fPath, varargin)
%analyzeIT Extract IT behavior as function of gradient position
%   Detailed explanation goes here

% Outputs
% fPath, file path for analyzed data
% hRat, fraction of time spent IT for x-dimension in nBins
% ITscore, MAKE THIS PER TRACK METRIC!!

% optional inputs
%   travThresh, objects with less total travel are excluded (schmutz filter)
%   smW, smoothing window, default 50 frames
%   lThresh, length of continuity required to be an IT
%   nBins, hMax, hMin: histogram is from hMin:hMax in nBins steps


% initialize optional inputs
travThresh=200; % length of total travel to be considered as a worm...
smW=1; % window for smoothing to find speed.50
contThresh=50; % super arbitrary atm
prefix='';
trimX=5; % trim values from the edges for final calls.
trimY=5; % trim values from the edges for final calls.
hMax=900; hMin=0; % max & min of histogram bin edges
nBins=10; % number of bins for histogram

varargin=assignApplicable(varargin);
% assignApplicable seems to choke on directory assignment...?
while ~isempty(varargin)
    eval([varargin{1},'=','''',varargin{2},''';']);
    varargin(1:2)=[];
end

%% Load a trial data set.

% if inappropriate input, ask for file.
if ~exist('fPath','var')||isempty(fPath)
    [fN, dDir]=uigetfile('.mat','Please select .mat file with eset');
    fPath=fullfile(dDir,fN);
end

% trim 'experiment_1', which is required for fromMatFiles
bors=findstr(fPath,'experiment_1');
if ~isempty(bors)
    fPath=fPath(1:(bors-2)); %
end

eset=ExperimentSet.fromMatFiles(fPath); % load eset for analysis
% % Give tracks values... Only required if Track-wise analysis used.
% eset.expt(1).segmentTracks();

%% Establish saving directory
% Figures & saving
if  ~exist('saveDir','var')||isempty(saveDir)
    saveDir=fileparts(fPath);
    % create unique save directory for each.
    saveDir=fullfile(fPath,'plots');
end
if ~exist(saveDir,'dir')
    mkdir(saveDir);
end

%% Extract migration information
% Get values for travel in isothermal direction, Yiso, and...
% travel in orthothermal direction, Ytemp.
% Use plotXt function to extract & visualize data.
% Ytemp is temperature dimension
[ ~, ~, ~, ~, Ytemp] = plotXt( eset.expt,...
    'dim',1,'ylab','orthothermal','yts',{});
minX=min(min(Ytemp)); maxX=max(max(Ytemp));
set(gca, 'ylim', [minX, maxX]);
% Yiso is isothermal dimension
[ ~, ~, ~, ~, Yiso ] = plotXt( eset.expt,...
    'dim',2,'ylab','isothermal','yts',{});  
minY=min(min(Yiso)); maxY=max(max(Yiso)); % Edges of data. Useful to trim borders
% NOTE: trim borders later
set(gca, 'ylim', [minY, maxY]);

% Get rid of the immobile bits using travThresh
% travThresh=200; % make non-arbitrary and dynamic based on distrubition?
travLrs=nansum(abs(diff(Yiso)))+nansum(abs(diff(Ytemp)))>travThresh;
y=Yiso(:,travLrs);
x=Ytemp(:,travLrs);

%% Quantify IT based on Y-dimension speed component
% Assay IT by proportion of speed in Iso dimension, fracITspeed.
% find speed in each direction, but use smoothing to make more reliable
% over smoothing window of:
% smW=50;
smt1=smoothdata(y,1,'gaussian',smW);
smt2=smoothdata(x,1,'gaussian',smW);
itDiff=diff(smt1);
orthoDiff=diff(smt2);
speed=sqrt(itDiff.^2+orthoDiff.^2);
fracITspeed=abs(itDiff)./speed;


%% Simplify for figures into linear arrays without nans
% trim X&Y to match speed for plotting purposes
% turn into single vectors for ease of plotting
ySpeed=fracITspeed(:);
speed=speed(:);
xDisp=x(1:end-1,:);
xDisp=xDisp(:);
yDisp=y(1:end-1,:);
yDisp=yDisp(:);

% cut out nan values for position.
inVals=~isnan(yDisp); % can use later with other features
ySpeed=ySpeed(inVals);
speed=speed(inVals);
xDisp=xDisp(inVals);
yDisp=yDisp(inVals);

outThresh=nanmean(speed)+5*nanstd(speed);
speed(speed>outThresh)=nan; % remove a few outliers (17/53k)


% Extract fracITspeed as function of temperature position
% [ITvalues,indexArray] = indexAbyB(ySpeed,xDisp);

%% Call IT based on criteria
% 1. Simple fracITspeed threshold. % Do data give me a clear threshold?
% 2. Continuity above threshold. Maintained for XX length.
% 3. Maintain minimum speed?


% Combine speed, ySpeed, and continuity to ID IT
speedThresh=0.1;
yspeedThresh=0.9;
contThresh=50; % Could instead filter by physical length of continuity. 100

% first filter
passVals=and((speed>speedThresh),(ySpeed>yspeedThresh));

% calculate continuity of passing first filter
% forgiveness... eliminate very short failures.
giveThresh=5; % stretches of less than this many failures are stitched together
[lFail, ~] = stretchCounter(~passVals);
passVals(lFail<giveThresh)=1;
% encode continuity after stitching
[lPass, ~] = stretchCounter(passVals);
%
passVals2=and(passVals,lPass>contThresh);

% trim values from edges:
xtrim= and( (minX+trimX) <xDisp, xDisp< (maxX-trimX) );
ytrim= and( (minY+trimY) <yDisp, yDisp< (maxY-trimY) );
xytrim= and(xtrim,ytrim);

isIT=and(passVals2,xytrim);
notIT=and(~passVals2,xytrim);

% Extract IT calls as function of temperature position, ITcount
[ITcount,~] = indexAbyB(isIT,xDisp);
[notITcount,~] = indexAbyB(notIT,xDisp);
[ITscore,~] = indexAbyB(ySpeed,xDisp);


h1v=histcounts(xDisp(isIT)-minX,hMin:(hMax-hMin)/10:hMax);
h2v=histcounts(xDisp(notIT)-minX,hMin:(hMax-hMin)/10:hMax);
hRat=h1v./(h1v+h2v);
Hist=cat(1,h1v,h2v,hRat);

% Create IT score as a function of track.
% first assign track number to each value.
trackNum=repmat(1:size(x,2),[size(x,1)-1,1]);
trackPos=repmat(1:size(x,1)-1,[size(x,2),1])';
% linearize
trackNum=trackNum(:);
trackPos=trackPos(:);
% remove positional nans, from above
trackNum=trackNum(inVals);
trackPos=trackPos(inVals);
% assign track values of isIT & notIT, back to a matrix... 
ITmat=nan(size(x));
notITmat=nan(size(x));
for ii=1:size(x,2)
    TOI=(trackNum==ii);
    POI=trackPos(TOI);
    ITmat(POI,ii)=isIT(TOI);
    notITmat(POI,ii)=notIT(TOI);
end

ITtrack=nansum(ITmat,1);
notITtrack=nansum(notITmat,1);
ITscore=ITtrack./(ITtrack+notITtrack);
%ITscore(isnan(ITscore))=0;
ITtrackscore=cat(1,ITtrack,notITtrack,ITscore);

%% Figures
% QC plot to make sure all is well on re-conversion
ITindex=ITmat;
ITindex(isnan(ITindex))=0;
ITindex=logical(ITindex);
xITindex=x.*ITindex;
yITindex=y.*ITindex;
nITindex=notITmat;
nITindex(isnan(nITindex))=0;
nITindex=logical(nITindex);
xnITindex=x.*nITindex;
ynITindex=y.*nITindex;
for k=1:size(xnITindex,2)
    saveName=strcat(char(prefix),'_ITcalls_',string(k),'.png');
    axis([800 1700 400 1400]);
    figure(); hold on;
    plot(xITindex(:,k),yITindex(:,k),'.r')
    plot(xnITindex(:,k),ynITindex(:,k),'.k')
    vectorSave(gcf, fullfile(saveDir,saveName));
    close(figure);
end

close all;

end