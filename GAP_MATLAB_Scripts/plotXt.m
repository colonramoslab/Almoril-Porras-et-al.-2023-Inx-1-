function [ fh, Xbar, Ybar, X, Y ] = plotXt( eset, varargin )
%plotXt Summary of this function goes here
%   Detailed explanation goes here


%key parameters
%   eset.track.startFrame
%   eset.track.endFrame
%  -where i=worm/track over 7199 frames & j=frame within track
%   eset.track(i).pt(j).ind %frame of point, aka time
%   eset.track(i).pt(j).loc % position of point: [x;y]

% diagnostic: collect all starting locations, very little variation
%for i=1:length(eset.track)
%    X(i)=eset.track(i).pt(1).loc(1); %position in X (Temperature
%    Y(i)=eset.track(i).pt(1).loc(2); % position with respect to Y, irrelevant
%end

%
field='loc'; % position of point [x,y];
dim=1; % x dimension or temperature in 'loc' field
duRun=7200;
ylab='gradient temperature';
yts={'18C', '19C', '20C', '21C', '22C'};

varargin=assignApplicable(varargin);

fh=[];

X=nan(duRun,length(eset(1).track));
Y=nan(duRun,length(eset(1).track));

for foo=1:numel(eset)
    for i=1:length(eset(foo).track)
        for j=1:length(eset(foo).track(i).pt);
            % create two matrices with Y & X values
            % each column is a single track, X(j,i) & Y(j,i)
            ptInd=eset(foo).track(i).pt(j).ind+1;
            X(ptInd,i)=ptInd;
            Y(ptInd,i)=eset(foo).track(i).pt(j).(field)(dim);
        end
    end
end

fh=figure();
plot(X,Y);
xmin=0;xmax=duRun;ymin=400;ymax=2200;
axis([xmin,xmax,ymin,ymax]);
ystep=450;
ticks=400:ystep:2200;
set(gca,'ytick',ticks)

ylabel(ylab);
set(gca,'ytick',ticks)
set(gca,'yticklabel',yts)

xstep= 1200; %15min
ticks=0:xstep:duRun;
labels=0:10:60;
set(gca,'xtick',ticks)
set(gca,'xticklabel',labels)
xlabel('migration time (min)');

%minimum is ~400 for cryophilic runs

Xbar=nanmean(X,2);
Ybar=nanmean(Y,2);

end