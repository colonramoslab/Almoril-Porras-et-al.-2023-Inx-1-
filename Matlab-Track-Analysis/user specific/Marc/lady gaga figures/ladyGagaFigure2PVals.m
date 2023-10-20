
descriptions = {'CS EA 20 ppm/cm', 'CS EA 0.2 ppm/cm', 'CS EA 2 ppm/cm, masked by 25 ppm EA', 'GR63a EA 30 ppm/cm', ...
    'GR63a EA 200 ppm/cm', 'CS EA 2 ppm/cm', 'OR83b1 EA 2 ppm/cm', 'OR83b2 EA 2 ppm/cm'};
labelnames = {'CS 0-5\%', 'Gr63a 0-5\%', 'CS 0-1\%', 'CS 0-0.5\%', 'CS air control', 'CS no air control', 'CS 0-2\%'};
%these data are created by ladyGagaSpatialEthylAcetate and
%ladyGagaSpatialCO2

load ('C:\Users\Marc\Documents\figures\lady gaga paper\ethyl acetate spatial\ea_spatial_calcs.mat');
load (fullfile('C:\Users\Marc\Documents\figures\lady gaga paper\co2 spatial', 'co2_spatial_calcs'));

pval = 0.01;
cutoff = norminv(1-pval, 0, 1);
ni = [ea_spatial_ad.navind];
ne = [ea_spatial_ad.navind_eb];
eagt01 = (ni - 0.1)./ne > cutoff;
eagt005 = (ni - 0.05)./ne > cutoff;
eagt00 = ni./ne > cutoff;
ealt00 = ni./ne < -cutoff;
ealt01 = (ni + 0.1)./ne < -cutoff;
ealt005 = (ni + 0.05)./ne < -cutoff;

pval = 0.01;
cutoff = norminv(1-pval, 0, 1);
ni = [co2_spatial_ad.navind];
ne = [co2_spatial_ad.navind_eb];
cogt01 = (ni - 0.1)./ne > cutoff;
cogt005 = (ni - 0.05)./ne > cutoff;
cogt00 = ni./ne > cutoff;
colt00 = ni./ne < -cutoff;
colt01 = (ni + 0.1)./ne < -cutoff;
colt005 = (ni + 0.05)./ne < -cutoff;


for j = 1:length(descriptions)
    xmsg = '-blank-';
    ymsg = xmsg;
    if (eagt00(1,j))
        xmsg = '+';
    end
    if (eagt005(1,j))
        xmsg = '++';
    end
    if (eagt01(1,j))
        xmsg = '+++';
    end
    if (eagt00(2,j))
        ymsg = '+';
    end
    if (eagt005(2,j))
        ymsg = '++';
    end
    if (eagt01(2,j))
        ymsg = '+++';
    end
    
    if (ealt00(1,j))
        xmsg = '-';
    end
    if (ealt005(1,j))
        xmsg = '--';
    end
    if (ealt01(1,j))
        xmsg = '---';
    end
    if (ealt00(2,j))
        ymsg = '-';
    end
    if (ealt005(2,j))
        ymsg = '--';
    end
    if (ealt01(2,j))
        ymsg = '---';
    end
    eamsg{j} = [descriptions{j} ' x: ', xmsg, ',  y: ', ymsg]; 
end


for j = 1:length(labelnames)
    xmsg = '-blank-';
    ymsg = xmsg;
    if (cogt00(1,j))
        xmsg = '+';
    end
    if (cogt005(1,j))
        xmsg = '++';
    end
    if (cogt01(1,j))
        xmsg = '+++';
    end
    if (cogt00(2,j))
        ymsg = '+';
    end
    if (cogt005(2,j))
        ymsg = '++';
    end
    if (cogt01(2,j))
        ymsg = '+++';
    end
    
    if (colt00(1,j))
        xmsg = '-';
    end
    if (colt005(1,j))
        xmsg = '--';
    end
    if (colt01(1,j))
        xmsg = '---';
    end
    if (colt00(2,j))
        ymsg = '-';
    end
    if (colt005(2,j))
        ymsg = '--';
    end
    if (colt01(2,j))
        ymsg = '---';
    end
    comsg{j} = [labelnames{j} ' x: ', xmsg, ',  y: ', ymsg]; 
end

msg = [eamsg comsg];
msg'