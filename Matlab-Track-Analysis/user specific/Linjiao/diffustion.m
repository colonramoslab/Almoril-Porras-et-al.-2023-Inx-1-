xb=[];
yb=[];
rb=[];
tb=[];
for i=1:length(eset.expt)
    for j=1:length(eset.expt(i).track)
        %if eset.expt(i).track(j).npts>1800 %&& eset.expt(i).track(j).startFrame<120
        if eset.expt(i).track(j).startFrame<120
           location=eset.expt(i).track(j).getDerivedQuantity('iloc');
           xk=location(1,:);
           yk=location(2,:);
           tk=eset.expt(i).track(j).getDerivedQuantity('eti');
           xbb=xk-xk(1);
           ybb=yk-yk(1);
           rbb=sqrt(xbb.^2+ybb.^2);
           tbb=tk-tk(1);
           xb=[xb, xbb];
           yb=[yb, ybb];
           rb=[rb, rbb];
           tb=[tb, tbb];
        end
    end
end
xb=xb./88.36; %pix to cm
yb=yb./88.36;
rb=rb./88.36;
xsq = xb.^2;
ysq = yb.^2;
rsq = rb.^2;
%timeAxis=0:60:1800;
timeAxis=0:30:900;

[et,mxs]=meanyvsx(tb,xsq,timeAxis);
figure;
plot(et,mxs,'.');

hold off;
title('mean x square over time');

[et,mx, stderrx, stdx]=meanyvsx(tb,xb,timeAxis);
figure;
%plot(et,mx,'.');
errorbar(et,mx,stdx,'.');
hold on;
% errorbar(et,mx,stderrx,'r.');
% title('mean x  over time');

[et,mrs]=meanyvsx(tb,rsq,timeAxis);
figure;
plot(et,mrs,'.');
title('mean r square over time');

[et,mys]=meanyvsx(tb,ysq,timeAxis);
figure;
plot(et,mys,'.');
hold off;
title('mean y square over time');


xAxis=-13:0.5:13;
ind0=find(tb >0 & tb < 60);
a0 = xb(ind0);
[y0,x0] = hist(a0,xAxis);
y0=y0./sum(y0);
ind15=find(tb >840 & tb < 900);
a15=xb(ind15);
[y15,x15]=hist(a15,xAxis);
y15=y15./sum(y15);
ind30=find(tb >1740 & tb < 1800);
a30=xb(ind30);
[y30,x30]=hist(a30,xAxis);
y30=y30./sum(y30);


% yAxis=-13:0.5:13;
% ind0=find(tb >0 & tb < 60);
% a0 = yb(ind0);
% [y0,x0] = hist(a0,yAxis);
% y0=y0./sum(y0);
% ind15=find(tb >840 & tb < 900);
% a15=yb(ind15);
% [y15,x15]=hist(a15,yAxis);
% y15=y15./sum(y15);
% ind30=find(tb >1740 & tb < 1800);
% a30=yb(ind30);
% [y30,x30]=hist(a30,yAxis);
% y30=y30./sum(y30);

%len=length(xb);
% xc=zeros(1,len); %calculate the <(x-x0)^2-<x-x0>^2>
% for i=1:length(xb)
%     meanx=mean(xb(tb==tb(i)));
%     xcc=(xb(i))^2-meanx^2;
%     xc(i)=xcc;
% end
% [et,rmxs]=meanyvsx(tb,xc,timeAxis);
% figure(5)
% plot(et,rmxs,'.');
% title('root mean x square over time');


% yc=zeros(1,len); %calculate the <(y-y0)^2-<y-y0>^2>
% for i=1:length(yb)
%     meany=mean(yb(tb==tb(i)));
%     ycc=(yb(i))^2-meany^2;
%     yc(i)=ycc;
% end
% [et,rmys]=meanyvsx(tb,yc,timeAxis);
% figure(5)
% plot(et,rmys,'.');
% title('root mean y square over time');
    
    
        