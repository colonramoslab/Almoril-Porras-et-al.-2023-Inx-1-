function [] = GusAutoAnalyze(filename, varargin)

outputdirectory = '';
binfile = '';
fileprefix = '';
printOpt=1;
fullOpt=1;
varargin=assignApplicable(varargin);
segmentOptions=WormSegmentOptions;

funPath='C:\Users\jshha\Dropbox\ColonRamosLab\Matlab\thermotaxisAnalysis\';
if exist(funPath)
    addpath(genpath(funPath));
else
    warning([funPath ': not found. May need to find Matlab-Track-Analysis Folder']);
end


[fid, msg] = fopen(filename);
if fid == (-1)
    error(message('JoshAutoAnalyze:cannotOpenFile', filename, msg));
end


createTimingInfoForFileSets(filename);

tline = fgetl(fid);
while ischar(tline) %loop until eof
    if(strncmp(tline,'    output file:',16))   %look for "output file" line
        remain = tline;
        
        %parse output file for directory and bin file
        [str, remain] = strtok(remain,':');
        [str, remain] = strtok(remain,' :');
        outputdirectory = strcat(outputdirectory,str);
        while true
            [str, remain] = strtok(remain,'\');
            if isempty(str), break; end
            if isempty(strfind(str,'.bin')) ~= 1
                binfile = strcat(outputdirectory,str);
                fileprefix = strtok(str,'.');
                break;
            end
            outputdirectory = strcat(outputdirectory,str);
            outputdirectory = strcat(outputdirectory,'\');
        end
        
        
        eset=ExperimentSet.fromFiles(binfile);
        
        
        %stitch tracks
        nathanAutoStitch
        if printOpt
            %plot tracks and save image
            [fignum]=ExpPlotTracks2(eset.expt); %
            set(fignum,'Units','Inches');
            pos=get(fignum,'Position');
            set(fignum,'PaperPositionMode','Auto','PaperUnits','Inches',...
                'PaperSize',[pos(3),pos(4)]);
            outName=strcat(outputdirectory,fileprefix,'_initial.pdf');
            print(fignum,outName,'-dpdf','-r150','-painters');
            close(fignum);
        end
        %save as Mat file
        eset.toMatFiles(strcat(outputdirectory,fileprefix,'_poststitch'));
        
        %prune tracks from off the edge
        for j = 1:length(eset.expt)
            minXPosition = eset.expt(j).track(1,1).dq.sloc(1,1);
            minYPosition = eset.expt(j).track(1,1).dq.sloc(2,1);
            maxXPosition = eset.expt(j).track(1,1).dq.sloc(1,1);
            maxYPosition = eset.expt(j).track(1,1).dq.sloc(2,1);
            for k = 1:length(eset.expt(j).track)
                for l = 1:length(eset.expt(j).track(1,k).dq.sloc)
                    if eset.expt(j).track(1,k).dq.sloc(1,l)<minXPosition
                        minXPosition = eset.expt(j).track(1,k).dq.sloc(1,l);
                    end
                    if eset.expt(j).track(1,k).dq.sloc(2,l)<minYPosition
                        minYPosition = eset.expt(j).track(1,k).dq.sloc(2,l);
                    end
                    if eset.expt(j).track(1,k).dq.sloc(1,l)>maxXPosition
                        maxXPosition = eset.expt(j).track(1,k).dq.sloc(1,l);
                    end
                    if eset.expt(j).track(1,k).dq.sloc(2,l)>maxYPosition
                        maxYPosition = eset.expt(j).track(1,k).dq.sloc(2,l);
                    end
                end
            end
            %areacheck = (maxXPosition-minXPosition)*(maxYPosition-minYPosition);
            %if areacheck > 828100 %910px x 910px area
            %        warning('Tracks found in an area larger than estimated arena');
            %elseif areacheck < 720801 %849px x 849px area
            %        warning('Tracks found in an area much smaller than estimated arena');
            %end
            eset.executeExperimentFunction('pruneTracks',[],[(minXPosition+100) (minYPosition+100) (maxXPosition-100) (maxYPosition-100)]);
        end
        
        %Calculate minPts value for cleaning tracks
        %Assume 22 worms
        ptscutoff = 100;
        %         numpts = [];
        %         for i = 1:length(eset)
        %             for j = 1:length(eset(i).expt);
        %                 for k = 1:length(eset(i).expt(j).track);
        %                     numpts = [numpts eset(i).expt(j).track(k).npts];
        %                 end
        %             end
        %         end
        %
        %         %sort numpts vector
        %         numpts = sort(numpts,'descend');
        %
        %
        %         if length(numpts) > 23
        %             ptscutoff = numpts(24);
        %         else
        %             ptscutoff = 0;
        %         end
        %
        
        %Clean tracks
        ecl=ESetCleaner();
        ecl.askFirst=false;
        ecl.minSpeed=0.25; % formerly 0.1 & 0.0001
        ecl.minPts=ptscutoff;
        ecl.minDist=50;
        ecl.clean(eset);
        
        %save as Mat file
        eset.toMatFiles(strcat(outputdirectory,fileprefix,'_postclean'));
        if fullOpt
            if printOpt
                %plot tracks and save image
                [fignum]=ExpPlotTracks2(eset.expt); %
                set(fignum,'Units','Inches');
                pos=get(fignum,'Position');
                axis([minXPosition maxXPosition minYPosition maxYPosition]);
                set(fignum,'PaperPositionMode','Auto','PaperUnits','Inches',...
                    'PaperSize',[pos(3),pos(4)]);
                outName=strcat(outputdirectory,fileprefix,'_final.pdf');
                print(fignum,outName,'-dpdf','-r150','-painters');
                close(fignum);
                
                %plot superimposed tracks and save image
                nathanPlotSuperimposed
                set(fignum,'Units','Inches');
                pos=get(fignum,'Position');
                set(fignum,'PaperPositionMode','Auto','PaperUnits','Inches',...
                    'PaperSize',[pos(3),pos(4)]);
                outName=strcat(outputdirectory,fileprefix,'_super.pdf');
                print(fignum,outName,'-dpdf','-r150','-painters');
                close(fignum);
            end
            %segment tracks
            eset.expt(1).segmentTracks(segmentOptions);
            
            %get data
            [textfile, msg] = fopen(strcat(outputdirectory,fileprefix,'_data.csv'),'w');
            if textfile == (-1)
                error(message('JoshAutoAnalyze:cannotOpenFile', filename, msg));
            end
            JoshGetData_WVI
            fclose(textfile);
            
            %get theta data
            [textfile, msg] = fopen(strcat(outputdirectory,fileprefix,'_theta_data.csv'),'w');
            if textfile == (-1)
                error(message('JoshAutoAnalyze:cannotOpenFile', filename, msg));
            end
            
            %[turnfile, msg] = fopen(strcat(outputdirectory,fileprefix,'_turn_data.csv'),'w');
            %if turnfile == (-1)
            %    error(message('JoshAutoAnalyze:cannotOpenFile', filename, msg));
            %end
            GusGetTheta
            %JoshGetThetaTurn
            fclose(textfile);
            %fclose(turnfile);
            
            
            if printOpt
                % plot timing of runs
                [fh,~,~,~,~]=plotXt(eset.expt);
                set(fh,'Units','Inches');
                pos=get(fh,'Position');
                set(fh,'PaperPositionMode','Auto','PaperUnits','Inches',...
                    'PaperSize',[pos(3),pos(4)]);
                outName=strcat(outputdirectory,fileprefix,'_timeTemp.pdf');
                print(fh,outName,'-dpdf','-r150','-painters');
                close(fh);
                
                % plot timing of runs
                [fh,~,~,~,~]=plotXt(eset.expt,'dim',2,'ylab','isothermal','yts',{});
                set(gca,'ylim',[0,2000]);
                set(fh,'Units','Inches');
                pos=get(fh,'Position');
                set(fh,'PaperPositionMode','Auto','PaperUnits','Inches',...
                    'PaperSize',[pos(3),pos(4)]);
                outName=strcat(outputdirectory,fileprefix,'_timeIso.pdf');
                print(fh,outName,'-dpdf','-r150','-painters');
                close(fh);
            end
        end
    end
    
    clear eset;
    outputdirectory = '';
    binfile = '';
    fileprefix = '';
    tline = fgetl(fid);
end

fclose(fid);



