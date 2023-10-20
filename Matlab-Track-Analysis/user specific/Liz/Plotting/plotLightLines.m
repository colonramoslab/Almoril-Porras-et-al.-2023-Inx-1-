function [ ] = plotLightLines(period, minY, maxY)
%function [ ] = plotLightLines(period)
%plots lines onto phototaxis plots to show when light is rising/falling

firstFrame=0;
lastFrame=1400;

hold on

%plots light rising lines in yellow (first frame is yellow since all videos
%start in the dark)
for j=firstFrame:period:lastFrame
    line([j j],[minY maxY],'LineWidth',2,'LineStyle','--','Color',[0 1 0]);
end

%plots light falling lines in grey 
for k=firstFrame+period/2:period:lastFrame
    line([k k],[minY maxY],'LineWidth',2,'LineStyle','--','Color',...
        [0.501960784313725 0.501960784313725 0.501960784313725]);
end

hold off

end

