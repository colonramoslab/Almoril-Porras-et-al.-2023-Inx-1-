function []= GAPTurnt3(dDir,initialX)

fileL=dir([dDir, '/*.mat']);


for i=1:length(fileL)
% Load & segment file
load(fullfile(fileL(i).folder,fileL(i).name));
segmentTracks(experiment_1);
expt=experiment_1;

    [turnfile, msg] = fopen(strcat(fullfile(fileL(i).folder,fileL(i).name),'_binned_turn_rate_data.csv'),'w');
            if turnfile == (-1)
                error(message('GAPTurnt:cannotOpenFile', filename, msg));
            end
Reotheta=0;
turnOmega=0;
Omegatheta=0;
turnBlip=0;
Bliptheta=0;
turnReversal=0;
Reversaltheta=0;
turnSecondRev=0;
SecondRevtheta=0;
TbinturnOmega=zeros(10,1);%GAP added on 5/12 for spatial analysis of turns
TbinturnBlip=zeros(10,1);%GAP added on 5/12 for spatial analysis of turns
TbinturnReversal=zeros(10,1);%GAP added on 5/12 for spatial analysis of turns
TbinturnSecondRev=zeros(10,1);%GAP added on 5/12 for spatial analysis of turns
bineti=NaN(10,1);%GAP added on 5/15 for spatial analysis of turn rate
framepass = 0;%GAP added on 5/16 for spatial analysis of turn rate
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
        fprintf(turnfile,'Expt,Track,Total time,Timecheck,# Reorientations,# Sharp Turns,Omega Turns,Blip Turns,Reversals,Second Reversals,AVG dTheta in Reorientations,AVG dTheta in Omega Turns,AVG dTheta in Blip Turns,AVG dTheta in Reversals,AVG dTheta in Second Reversals,Mins in bin(1),Mins in bin(2),Mins in bin(3),Mins in bin(4),Mins in bin(5),Mins in bin(6),Mins in bin(7),Mins in bin(8),Mins in bin(9),Mins in bin(10),Omega Turns (1),Omega Turns (2),Omega Turns (3),Omega Turns (4),Omega Turns (5),Omega Turns (6),Omega Turns (7),Omega Turns (8),Omega Turns (9),Omega Turns (10),Blip Turns (1),Blip Turns (2),Blip Turns (3),Blip Turns (4),Blip Turns (5),Blip Turns (6),Blip Turns (7),Blip Turns (8),Blip Turns (9),Blip Turns (10),Reversals (1),Reversals (2),Reversals (3),Reversals (4),Reversals (5),Reversals (6),Reversals (7),Reversals (8),Reversals (9),Reversals (10),Second Reversals (1),Second Reversals (2),Second Reversals (3),Second Reversals (4),Second Reversals (5),Second Reversals (6),Second Reversals (7),Second Reversals (8),Second Reversals (9),Second Reversals (10)\n');

        for k = 1:length(expt(j).track)
            %Get # of sharp turns
            numSharpTurns = length(expt(j).track(1,k).sharpTurn);
            numReo = length(expt(j).track(1,k).reorientation);
            Totaleti = expt(j).track(1,k).dq.eti(1,end)-(expt(j).track(1,k).dq.eti(1,1));
            
            for m = 1:length(expt(j).track(1,k).reorientation)
                Reotheta = Reotheta + abs(expt(j).track(1,k).reorientation(1,m).dTheta);
            end
            ReothetaAVG = Reotheta/numReo;
            
            for m = 1:length(expt(j).track(1,k).sharpTurn)
                if expt(j).track(1,k).sharpTurn(1,m).typeCode == -1
                    turnOmega = turnOmega + 1;
                    Omegatheta = Omegatheta + abs(expt(j).track(1,k).sharpTurn(1,m).dTheta);
                elseif expt(j).track(1,k).sharpTurn(1,m).typeCode == 0
                    turnBlip = turnBlip + 1;
                    Bliptheta = Bliptheta + abs(expt(j).track(1,k).sharpTurn(1,m).dTheta);
                elseif expt(j).track(1,k).sharpTurn(1,m).typeCode == 1
                    turnReversal = turnReversal + 1;
                    Reversaltheta = Reversaltheta + abs(expt(j).track(1,k).sharpTurn(1,m).dTheta);
                elseif expt(j).track(1,k).sharpTurn(1,m).typeCode == 2
                    turnSecondRev = turnSecondRev +1;
                    SecondRevtheta = SecondRevtheta + abs(expt(j).track(1,k).sharpTurn(1,m).dTheta);
                end
            end
            OmegathetaAVG = Omegatheta/turnOmega;
            BlipthetaAVG = Bliptheta/turnBlip;
            ReversalthetaAVG = Reversaltheta/turnReversal;
            SecondRevthetaAVG = SecondRevtheta/turnSecondRev;
            
            for m = 1:length(expt(j).track(1,k).sharpTurn)
                for n=1:10
                    if expt(j).track(1,k).sharpTurn(1,m).loc(1,1)>(initialX+90*(n-1)) && expt(j).track(1,k).sharpTurn(1,m).loc(1,1)<=(initialX+90*(n))
                        if expt(j).track(1,k).sharpTurn(1,m).typeCode == -1
                            TbinturnOmega(n) = TbinturnOmega(n) + 1;
                        elseif expt(j).track(1,k).sharpTurn(1,m).typeCode == 0
                            TbinturnBlip(n) = TbinturnBlip(n) + 1;
                        elseif expt(j).track(1,k).sharpTurn(1,m).typeCode == 1
                            TbinturnReversal(n) = TbinturnReversal(n) + 1;
                        elseif expt(j).track(1,k).sharpTurn(1,m).typeCode == 2
                            TbinturnSecondRev(n) = TbinturnSecondRev(n) +1;
                        end
                    end
                end
            end
            for n=1:10
                for m = 1:length(expt(j).track(1,k).dq.sloc)
                    if expt(j).track(1,k).dq.sloc(1,m)>(initialX+90*(n-1)) && expt(j).track(1,k).dq.sloc(1,m)<=(initialX+90*(n))
                        framepass = framepass + 1;
                    end
                end
                bineti(n,1) = framepass/120;
                framepass = 0;
            end
            totalbineticheck = sum(bineti)*60;
            TbinturnOmega = TbinturnOmega./bineti;
            TbinturnBlip = TbinturnBlip./bineti;
            TbinturnReversal = TbinturnReversal./bineti;
            TbinturnSecondRev = TbinturnSecondRev./bineti;
            fprintf(turnfile,'%s,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i\n',expt(j).fname,k,Totaleti,totalbineticheck,numReo,numSharpTurns,turnOmega,turnBlip,turnReversal,turnSecondRev,ReothetaAVG,OmegathetaAVG,BlipthetaAVG,ReversalthetaAVG,SecondRevthetaAVG,bineti(1),bineti(2),bineti(3),bineti(4),bineti(5),bineti(6),bineti(7),bineti(8),bineti(9),bineti(10),TbinturnOmega(1),TbinturnOmega(2),TbinturnOmega(3),TbinturnOmega(4),TbinturnOmega(5),TbinturnOmega(6),TbinturnOmega(7),TbinturnOmega(8),TbinturnOmega(9),TbinturnOmega(10),TbinturnBlip(1),TbinturnBlip(2),TbinturnBlip(3),TbinturnBlip(4),TbinturnBlip(5),TbinturnBlip(6),TbinturnBlip(7),TbinturnBlip(8),TbinturnBlip(9),TbinturnBlip(10),TbinturnReversal(1),TbinturnReversal(2),TbinturnReversal(3),TbinturnReversal(4),TbinturnReversal(5),TbinturnReversal(6),TbinturnReversal(7),TbinturnReversal(8),TbinturnReversal(9),TbinturnReversal(10),TbinturnSecondRev(1),TbinturnSecondRev(2),TbinturnSecondRev(3),TbinturnSecondRev(4),TbinturnSecondRev(5),TbinturnSecondRev(6),TbinturnSecondRev(7),TbinturnSecondRev(8),TbinturnSecondRev(9),TbinturnSecondRev(10));
            % added "posTimeAvg, negTimeAvg, posVelocityAvg, negVelocityAvg, posCurveIndexAvg, negCurveIndexAvg, posWVIAvg, negWVIAvg,StartCountBias,StartLengthBias" before ttxi
            % also included expt name in each line
            Reotheta=0;
            turnOmega=0;
            Omegatheta=0;
            turnBlip=0;
            Bliptheta=0;
            turnReversal=0;
            Reversaltheta=0;
            turnSecondRev=0;
            SecondRevtheta=0;
            ReothetaAVG=0;
            OmegathetaAVG = 0;
            BlipthetaAVG = 0;
            ReversalthetaAVG = 0;
            SecondRevthetaAVG = 0;
            
            
            TbinturnOmega=zeros(10,1);%GAP added on 5/12 for spatial analysis of turns
            TbinturnBlip=zeros(10,1);%GAP added on 5/12 for spatial analysis of turns
            TbinturnReversal=zeros(10,1);%GAP added on 5/12 for spatial analysis of turns
            TbinturnSecondRev=zeros(10,1);%GAP added on 5/12 for spatial analysis of turns
            bineti=NaN(10,1);%GAP added on 5/15 for spatial analysis of turns
        end
    end
            fclose(turnfile);
end
end

