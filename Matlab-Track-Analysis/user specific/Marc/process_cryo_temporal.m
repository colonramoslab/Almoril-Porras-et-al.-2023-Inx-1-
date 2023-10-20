basedir = 'D:\Marc Processed\cryo temporal n2\';
dates = {'1029','1023','1030'};

maxToLoad = 1;
if ~exist('expt', 'var')
    for j = 1:(min(maxToLoad,length(dates)))
        fn = [basedir '2009' dates{j} '_n2_g15_tracks'];
        expt(j) = Experiment.fromFile([fn '.bin'], [fn '.tim'], false, [], 10);
    end
end

if (~exist('fixedt', 'var'))
    for j = 1:length(expt)
        fn = [basedir '2009' dates{j} '_n2_g15_tracks'];
        data = load([fn '.tmp']);
        temperature{j} = data(end,:);
        tt = lowpass1D(temperature{j} - mean(temperature{j}), 400);
        t = 1:length(tt);
        %func = @(x,xdata) x(1)*cos(2*pi*x(2)*xdata + x(3));
        amp = sqrt(2*mean(tt.^2));
        inds = find(diff(sign(tt)) > 0);
        freq = 2*pi/median(diff(inds));
        phase = inds(1)+pi/2;
        phase = atan2(sin(phase), cos(phase));
    
        func = @(x,xdata) amp*cos(freq*xdata + phase);
        phase = lsqcurvefit(func, phase, t, tt);
        func = @(x,xdata) amp*cos(x(1)*xdata + phase);
        op = optimset('lsqcurvefit');
        op.TolX = 1E-9;
    %   op.Algorithm = 'levenberg-marquardt';
        freq = lsqcurvefit(func, freq, t, tt,[],[],op);
        x = [amp, freq*length(t), phase];
    
        lb = [0 0 -pi];
        ub = [Inf length(t) pi];
   
    
        op.Jacobian = 'on';
        %op.PlotFcns = @optimplotresnorm;
        op.DerivativeCheck = 'on';
        x2 = lsqcurvefit(@cosineForFitting, x, t, tt, lb, ub, op);
        figure(1);
        plot (t,tt,t,cosineForFitting(x,t),t,cosineForFitting(x2,t))
        figure(2);
        plot (t, tt-cosineForFitting(x2,t));
        p = polyfit (t,tt-cosineForFitting(x2,t),1);
        x3 = lsqcurvefit(@cosineForFitting, x2, t, tt - polyval(p, t), lb, ub, op);
        figure(3);
        plot (t,temperature{j}, 'b-', t, cosineForFitting(x3,t)+polyval(p,t)+mean(temperature{j}), 'r-');
        fixedt{j} = cosineForFitting(x3,t)+polyval(p,t)+mean(temperature{j});
        interval{j} = mean(diff(data(end-1,:)))/1000;
    end
end

for j = 1:length(expt)
    expt(j).executeTrackFunction('addGlobalQuantity', 'temperature', expt(j).elapsedTime, fixedt{j});
    dtemp = deriv(fixedt{j}, 1)/interval{j};
    expt(j).executeTrackFunction('addGlobalQuantity', 'dtemperature', expt(j).elapsedTime, dtemp);
end