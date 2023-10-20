% calculate the speed over time
timeAxis=0:120:3600;
runsp = eset.gatherFromSubField('run', 'speed', 'trimpts', 5); % pixle per second
runlocation= eset.gatherFromSubField('run', 'sloc', 'trimpts', 5);
runx=runlocation(1,:); % x location;
runt=(runx-200)*4/2200+18; %temperature 18C at 200 and 22C at 2400
runc=(runx-200)*50/2200; %salt concentration 50mM at 200 and 100mM at 2400
runsp = runsp/8.836; % mm per second

runtime = eset.gatherFromSubField('run', 'eti', 'trimpts', 5);
[newtime, meanspeed, stderrorspeed,stdspeed, sumy] = meanyvsx(runtime, runsp, timeAxis);
sizey=sumy./meanspeed;
err=stdspeed./sqrt(sizey/100); % independent n= total time/2 decorr. time=total frame/4*25
figure;
errorbar(newtime, meanspeed, err);


tempAxis=18:0.5:24;
cheAxis=0:2.5:50;
[newtemp, meanspeed, stderrorspeed,stdspeed, sumy] = meanyvsx(runt, runsp, tempAxis);
[newchem, meanspeed2, stderrorspeed2,stdspeed2, sumy2] = meanyvsx(runc, runsp, cheAxis);
sizey=sumy./meanspeed;
err=stdspeed./sqrt(sizey/100); % independent n= total time/2 decorr. time=total frame/4*25
figure;
errorbar(newtemp, meanspeed, err);

sizey2=sumy2./meanspeed2;
err2=stdspeed2./sqrt(sizey2/100);
figure;
errorbar(newchem, meanspeed2, err2);


