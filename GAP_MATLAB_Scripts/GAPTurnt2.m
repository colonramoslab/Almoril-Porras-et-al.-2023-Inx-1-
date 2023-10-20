function []= GAPTurnt2(dDir)

fileL=dir([dDir, '/*.mat']);


for i=1:length(fileL)
% Load & segment file
load(fullfile(fileL(i).folder,fileL(i).name));
segmentTracks(experiment_1);
expt=experiment_1;

    [turnfile, msg] = fopen(strcat(fullfile(fileL(i).folder,fileL(i).name),'_turn_time_data.csv'),'w');
            if turnfile == (-1)
                error(message('GAPTurnt:cannotOpenFile', filename, msg));
            end
turnOmega=0;
turnBlip=0;
turnReversal=0;
turnSecondRev=0;
fracturnOmega=zeros(12,1);%GAP added on 4/19 for time analysis of turns
fracturnBlip=zeros(12,1);%GAP added on 4/19 for time analysis of turns
fracturnReversal=zeros(12,1);%GAP added on 4/19 for time analysis of turns
fracturnSecondRev=zeros(12,1);%GAP added on 4/19 for time analysis of turns
l=0;
%outputFile = input('Type address of output file:', 's');

if exist('turnfile','var')
else
    turnfile=1;
end


    for j = 1:length(expt)
        % fprintf(textfile,'%s\n',expt(j).fname); MOVED to each
        % line to improve compiling multiple files with uniquely
        % identifying info on each line
        fprintf(turnfile,'Expt,Track,Total time,# Sharp Turns,Omega Turns,Blip Turns,Reversals,Second Reversals,Omega Turns (1-6min),Omega Turns (6-11min),Omega Turns (11-16min),Omega Turns (16-21min),Omega Turns (21-26min),Omega Turns (26-31min),Omega Turns (31-36min),Omega Turns (36-41min),Omega Turns (41-46min),Omega Turns (46-51min),Omega Turns (51-56min),Omega Turns (56-60min),Blip Turns (1-6min),Blip Turns (6-11min),Blip Turns (11-16min),Blip Turns (16-21min),Blip Turns (21-26min),Blip Turns (26-31min),Blip Turns (31-36min),Blip Turns (36-41min),Blip Turns (41-46min),Blip Turns (46-51min),Blip Turns (51-56min),Blip Turns (56-60min),Reversals (1-6min),Reversals (6-11min),Reversals (11-16min),Reversals (16-21min),Reversals (21-26min),Reversals (26-31min),Reversals (31-36min),Reversals (36-41min),Reversals (41-46min),Reversals (46-51min),Reversals (51-56min),Reversals (56-60min),Second Reversals (1-6min),Second Reversals (6-11min),Second Reversals (11-16min),Second Reversals (16-21min),Second Reversals (21-26min),Second Reversals (26-31min),Second Reversals (31-36min),Second Reversals (36-41min),Second Reversals (41-46min),Second Reversals (46-51min),Second Reversals (51-56min),Second Reversals (56-60min)\n');

        for k = 1:length(expt(j).track)
            %Get # of sharp turns
            numSharpTurns = length(expt(j).track(1,k).sharpTurn);
            Totaleti = expt(j).track(1,k).dq.eti(1,end)-(expt(j).track(1,k).dq.eti(1,1));
            
            for m = 1:length(expt(j).track(1,k).sharpTurn)
                if expt(j).track(1,k).sharpTurn(1,m).typeCode == -1
                    turnOmega = turnOmega + 1;
                elseif expt(j).track(1,k).sharpTurn(1,m).typeCode == 0
                    turnBlip = turnBlip + 1;
                elseif expt(j).track(1,k).sharpTurn(1,m).typeCode == 1
                    turnReversal = turnReversal + 1;
                elseif expt(j).track(1,k).sharpTurn(1,m).typeCode == 2
                    turnSecondRev = turnSecondRev +1;
                end
            end
            for m = 1:length(expt(j).track(1,k).sharpTurn)
                for n=1:12
                    if expt(j).track(1,k).sharpTurn(1,m).startInd>(120+120*5*(n-1)) && expt(j).track(1,k).sharpTurn(1,m).startInd<=(120+(120*5*n))
                        if expt(j).track(1,k).sharpTurn(1,m).typeCode == -1
                            fracturnOmega(n) = fracturnOmega(n) + 1;
                        elseif expt(j).track(1,k).sharpTurn(1,m).typeCode == 0
                            fracturnBlip(n) = fracturnBlip(n) + 1;
                        elseif expt(j).track(1,k).sharpTurn(1,m).typeCode == 1
                            fracturnReversal(n) = fracturnReversal(n) + 1;
                        elseif expt(j).track(1,k).sharpTurn(1,m).typeCode == 2
                            fracturnSecondRev(n) = fracturnSecondRev(n) +1;
                        end
                    end
                end
            end
            fprintf(turnfile,'%s,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i\n',expt(j).fname,k,Totaleti,numSharpTurns,turnOmega,turnBlip,turnReversal,turnSecondRev,fracturnOmega(1),fracturnOmega(2),fracturnOmega(3),fracturnOmega(4),fracturnOmega(5),fracturnOmega(6),fracturnOmega(7),fracturnOmega(8),fracturnOmega(9),fracturnOmega(10),fracturnOmega(11),fracturnOmega(12),fracturnBlip(1),fracturnBlip(2),fracturnBlip(3),fracturnBlip(4),fracturnBlip(5),fracturnBlip(6),fracturnBlip(7),fracturnBlip(8),fracturnBlip(9),fracturnBlip(10),fracturnBlip(11),fracturnBlip(12),fracturnReversal(1),fracturnReversal(2),fracturnReversal(3),fracturnReversal(4),fracturnReversal(5),fracturnReversal(6),fracturnReversal(7),fracturnReversal(8),fracturnReversal(9),fracturnReversal(10),fracturnReversal(11),fracturnReversal(12),fracturnSecondRev(1),fracturnSecondRev(2),fracturnSecondRev(3),fracturnSecondRev(4),fracturnSecondRev(5),fracturnSecondRev(6),fracturnSecondRev(7),fracturnSecondRev(8),fracturnSecondRev(9),fracturnSecondRev(10),fracturnSecondRev(11),fracturnSecondRev(12));
            % added "posTimeAvg, negTimeAvg, posVelocityAvg, negVelocityAvg, posCurveIndexAvg, negCurveIndexAvg, posWVIAvg, negWVIAvg,StartCountBias,StartLengthBias" before ttxi
            % also included expt name in each line
            turnOmega=0;
            turnBlip=0;
            turnReversal=0;
            turnSecondRev=0;
            fracturnOmega=zeros(12,1);%GAP added on 4/19 for time analysis of turns
            fracturnBlip=zeros(12,1);%GAP added on 4/19 for time analysis of turns
            fracturnReversal=zeros(12,1);%GAP added on 4/19 for time analysis of turns
            fracturnSecondRev=zeros(12,1);%GAP added on 4/19 for time analysis of turns
        end
    end
            fclose(turnfile);
end
end

