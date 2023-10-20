for i = 1:length(eset)
    for j = 1:length(eset(i).expt);
            %%Finding the max and min values of all tracks
            fignum=figure();
            clf(fignum);  %create a figure and clear it
            hold on;
            axis([0 2500 0 2000]);
            yL = get(gca,'YLim');
            minXPosition = eset(i).expt(j).track(1,1).pt(1,1).imOffset(1);
            minYPosition = eset(i).expt(j).track(1,1).pt(1,1).imOffset(2);
            maxXPosition = eset(i).expt(j).track(1,1).pt(1,1).imOffset(1);
            maxYPosition = eset(i).expt(j).track(1,1).pt(1,1).imOffset(2);
            for k = 1:length(eset(i).expt(j).track)
                for l = 1:length(eset(i).expt(j).track(1,k).pt)
                    if eset(i).expt(j).track(1,k).pt(1,l).imOffset(1)<minXPosition
                        minXPosition = eset(i).expt(j).track(1,k).pt(1,l).imOffset(1);
                    end
                    if eset(i).expt(j).track(1,k).pt(1,l).imOffset(2)<minYPosition
                        minYPosition = eset(i).expt(j).track(1,k).pt(1,l).imOffset(2);
                    end
                    if eset(i).expt(j).track(1,k).pt(1,l).imOffset(1)>maxXPosition
                        maxXPosition = eset(i).expt(j).track(1,k).pt(1,l).imOffset(1);
                    end
                    if eset(i).expt(j).track(1,k).pt(1,l).imOffset(2)>maxYPosition
                        maxYPosition = eset(i).expt(j).track(1,k).pt(1,l).imOffset(2);
                    end
                end
            end
            areacheck = (maxXPosition-minXPosition)*(maxYPosition-minYPosition);
            if areacheck > 828100 %910px x 910px area
                    warning(['Tracks found in an area larger than estimated arena']);
            elseif areacheck < 720801 %849px x 849px area
                    warning(['Tracks found in an area much smaller than estimated arena']);
            end
            

        for k = 1:length(eset(i).expt(j).track);
            
            %Get distance traveled (x displacement)
            startingXPosition = eset(i).expt(j).track(1,k).pt(1,1).imOffset(1);
            finalXPosition = eset(i).expt(j).track(1,k).pt(1,length(eset(i).expt(j).track(1,k).pt)).imOffset(1);
            trackLength = 0;
            for l = 1:length(eset(i).expt(j).track(1,k).run)
                trackLength = trackLength + eset(i).expt(j).track(1,k).run(1,l).pathLength;
            end
                if trackLength<100 %arbitrary
                    plotPath(eset(i).expt(j).track(k),'sloc','k-')
                end
                if startingXPosition < (minXPosition+20) %20px from minX
                    plotPath(eset(i).expt(j).track(k))
                end
                if startingXPosition > (maxXPosition-20) %20px from maxX
                    plotPath(eset(i).expt(j).track(k),'sloc','r-')
                end
            totalRunLength = 0;
            posXDisp = 0;
            negXDisp = 0;
            posRunLength = 0;
            negRunLength = 0;
            posRuns = 0;
            negRuns = 0;
            turnOmega=0;
            turnBlip=0;
            turnReversal=0;
            turnSecondRev=0;
            posBiasCount = 0;   
            negBiasCount = 0; 
            posRunTime = 0;            
            negRunTime = 0;
            posVelocity = 0;
            negVelocity = 0;
            posCurveIndex = 0;
            negCurveIndex = 0;
            WVIndex = 0;
            posWVIndex = 0;
            negWVIndex = 0;
        end
    end
end
