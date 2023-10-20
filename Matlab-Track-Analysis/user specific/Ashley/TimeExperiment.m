function [x1,y1]=TimeExperiment(expt, dt)

t=expt.gatherField('timeon');
v=expt.gatherField('speed');

[x1,y1]=meanyvsx(t,v,0:2:(2*dt));

h=expt.gatherField('hsstart');



end