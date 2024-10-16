function [] = JoshAutoAnalyzePDF_time(filename, varargin)

outputdirectory = '';
binfile = '';
fileprefix = '';
printOpt=1;
fullOpt=1;
varargin=assignApplicable(varargin);

funPath='C:\Users\dgmdi\OneDrive - Yale University\Documents\BasicTTXanalysis\'; %% before was C:\Users\jshha\Dropbox\ColonRamosLab\Matlab\thermotaxisAnalysis\
if exist(funPath)
    addpath(genpath(funPath));
else
    warning([funPath ': not found. May need to find Matlab-Track-Analysis Folder']);
end


[fid, msg] = fopen(filename);
if fid == (-1)
    error(message('JoshAutoAnalyze:cannotOpenFile', filename, msg));
end

%disp(filename)
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
            cntN=0;
            while exist(outName)
                cntN=cntN+1;
                outName=[outputdirectory fileprefix num2str(cntN) '_initial.pdf'];
            end
            print(fignum,outName,'-dpdf','-r150','-painters');
            close(fignum);
        end
        %save as Mat file
        eset.toMatFiles(strcat(outputdirectory,fileprefix,'_poststitch'));
        
        %prune tracks from off the edge
        eset.executeExperimentFunction('pruneTracks',[],[600 100 2100 1700]);
        
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
        ecl.minSpeed=0.0001; % formerly 0.1
        ecl.minPts=ptscutoff;
        ecl.minDist=0;
        ecl.clean(eset);
        
        %save as Mat file
        eset.toMatFiles(strcat(outputdirectory,fileprefix,'_postclean'));
        if fullOpt
            if printOpt
                %plot tracks and save image
                [fignum]=ExpPlotTracks2(eset.expt); %
                set(fignum,'Units','Inches');
                pos=get(fignum,'Position');
                set(fignum,'PaperPositionMode','Auto','PaperUnits','Inches',...
                    'PaperSize',[pos(3),pos(4)]);
                outName=strcat(outputdirectory,fileprefix,'_final.pdf');
                cntN=0;
                while exist(outName)
                    cntN=cntN+1;
                    outName=[outputdirectory fileprefix num2str(cntN) '_initial.pdf'];
                end
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
            eset.expt(1).segmentTracks();
            
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
            
            [turnfile, msg] = fopen(strcat(outputdirectory,fileprefix,'_turn_data.csv'),'w');
            if turnfile == (-1)
                error(message('JoshAutoAnalyze:cannotOpenFile', filename, msg));
            end
            JoshGetTheta
            JoshGetThetaTurn
            fclose(textfile);
            fclose(turnfile);
            
            
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
                cntN=0;
                while exist(outName)
                    cntN=cntN+1;
                    outName=[outputdirectory fileprefix num2str(cntN) '_initial.pdf'];
                end
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
print('Finished');


