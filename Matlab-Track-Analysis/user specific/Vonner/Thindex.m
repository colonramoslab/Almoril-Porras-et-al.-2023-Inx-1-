%% For Canton S flies (L1, L2, and L3)
%actual x1=[17.30,20.00,22.70,24.90,27.60,30.05]; %midpoint of actual temp range for each eset for L1s: 14.8:19.8, 17.6:22.4, 20.2:25.2, 22.5:27.3, 25.2:30.0, 27.4:32.7
%actual x2=[17.50,20.05,22.50,25.15,27.50,29.90]; %midpoint of actual temp range for each eset for L2s: 14.9:20.1, 17.5:22.6, 20.1:24.9, 22.6:27.7, 25.1:29.9, 27.5:32.3
%actual x3=[17.65,19.95,22.70,24.90,27.50,30.10]; %midpoint of actual temp range for each eset for L3s: 15.1:20.2, 17.4:22.5, 20.2:25.2, 22.4:27.4, 25.0:30.0, 27.6:32.6
x1=[17.5,20,22.5,25,27.5,30,32.5,35,37.5];
x2=[17.5,20,22.5,25,27.5,30,32.5,35,37.5];
x3=[17.5,20,22.5,25,27.5,30,32.5,35,37.5];
y1=[0.5134,.7612,.6899,.2542,.1918,.0681,-.1492,-.5123,-.6731]; %thindex for each expt on L1s; %input matrix (for 062510, thindex=+0.0488)
y2=[.5134,.4250,.3093,-0.0795,-.3227,-0.0304,0.0048,-.2915,-.6779]; %thindex for each expt on L2s; %input matrix (for 062510, thindex=+0.0934)
y3=[0.3159,-0.0451,-0.4262,-0.3107,-0.3486,-0.3521,-0.5103,-0.3605,-0.3157]; %thindex for each expt on L3s; %input matrix (for 062510, thindex=-0.3018)
subplot(1,3,1), plot(x1,y1);
axis([15,40,-1,1]);
xlabel('start temperature (degrees C)');
ylabel('thermotaxis index');
title('CS L1');
subplot(1,3,2), plot(x2,y2);
axis([15,40,-1,1]);
xlabel('start temperature (degrees C)');
ylabel('thermotaxis index');
title('CS L2');
subplot(1,3,3), plot(x3,y3);
axis([15,40,-1,1]);
xlabel('start temperature (degrees C)');
ylabel('thermotaxis index');
title('CS L3');


%% For TRPL[302] and dTRPA1[ins] mutants (L1 only)
figure(2);
x1=[20,22.5,25,27.5,30,32.5,35,37.5];
y1=[0.6519,-.2085,-.8219,-0.6190,-0.5180,-.4320,-.3171,0.0398]; %for TrpL[302] L1
x2=[17.5,20,22.5,25,27.5,30,32.5,35,37.5];
y2=[.7842,.4789,.8158,.7444,.3967,.4589,-.3004,-.3233,-.4681]; %for dTrpA1[ins]
subplot(1,2,1), plot(x1,y1);
axis([15,40,-1,1]);
xlabel('start temperature (degrees C)');
ylabel('thermotaxis index');
title('TrpL[302] L1');
subplot(1,2,2), plot(x2,y2);
axis([15,40,-1,1]);
xlabel('start temperature (degrees C)');
ylabel('thermotaxis index');
title('dTrpA1[ins] L1');