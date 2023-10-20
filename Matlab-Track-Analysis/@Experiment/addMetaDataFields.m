function addMetaDataFields(expt, timfname, varargin)
%parses the metadata file (.mdat) for global quantities

if (length(expt) > 1)
    existsAndDefault('timfname', []);
    for j = 1:length(expt)
        expt(j).addMetaDataFields(timfname, varargin);
    end
    return;
end
        

existsAndDefault('timfname', expt.timfname);

header = {'GAS_CH1_SP',  'co2filt_0', 'co2filt_1', 'PID_PPM', 'LED_ON'};
timheader = {'GAS_CH1_SP_time', 'co2fast_0_time', 'co2fast_1_time', 'PID_PPM_time', 'LED_ON_time'};
gname = {'gassp', 'co2ppm0', 'co2ppm1', 'vocppm', 'ledOn'};
smoothTime = [0, 2, 2, 7, 0];
derivTime = [0.1, 1, 1, 1, 0];
isgas = [true, true, true, true, false];
gasvel = 1.15; %changed from 1.5 by mhg 4/12
gasdir = [1;0]; 
gasorigin = [-29;0]; %change from -15 by mhg4/12
fnum = 1:length(expt.elapsedTime);
et = expt.elapsedTime;

try
    datastruct = importdata(timfname);
catch me
    disp(me);
    return;
end

for j = 1:length(header)
    col = find(strcmpi(datastruct.colheaders, header{j}), 1);
    if (~isempty(col))
        try
            ydata = datastruct.data(:,col);
            col = find(strcmpi(datastruct.colheaders, timheader{j}), 1);
            if (~isempty(col))
                timdata = datastruct.data(:,col)/1000;
            else
                timdata = et;
            end
            inds = find(isfinite(ydata)&isfinite(timdata));
            ydat = ydata(inds);
            timdata = timdata(inds);
            [~,I] = min(timdata);
            ydat = ydat(I:end); %get rid of any early bad time values before reset
            timdata = timdata(I:end);
            [timdata,I] = unique(timdata);
            ydat = ydat(I);
            p = polyfit(timdata, et(inds(I)),1);
            et2 = p(1)*timdata + p(2);
            eti = linspace(min(et2), max(et2), (max(et2)-min(et2))/expt.dr.interpTime);
            xdat = eti;
            ydat = interp1(et2, ydat, eti);
            if (smoothTime(j) > 0)
                ydat = lowpass1D(ydat, smoothTime(j)/expt.dr.interpTime, 'padType', 'linear');
            end
            gq = GlobalQuantity();
            gq.yData = ydat;
            gq.fieldname = gname{j};
            if (isgas(j))
                gq.xField = {'eti', 'sloc'};
                xData.et = xdat;
                xData.origin = gasorigin;
                xData.flowspeed = gasvel;
                xData.flowdir = gasdir;
                gq.xData = xData;
                gq.derivationMethod = @GlobalQuantity.timeVaryingGasDerivation;
            else
                gq.xField = 'eti';
                gq.xData = xdat;
            end
        
            
        
            expt.addGlobalQuantity(gq);
            if (derivTime(j) > 0)
                gq.xData = [];
                gq.yData = [];
                gq.xField = gname{j};
                gq.fieldname = ['d' gname{j}];
                gq.derivationMethod = @(xin, xData, yData) deriv(xin, derivTime(j)/expt.dr.interpTime);
                expt.addGlobalQuantity(gq);
            end
        catch me
            disp (['problem adding metadata field ' gname{j}]);
            disp (me.getReport);
        end
            
    end
end
