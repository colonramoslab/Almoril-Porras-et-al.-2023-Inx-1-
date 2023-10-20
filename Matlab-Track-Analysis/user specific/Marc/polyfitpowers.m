function [xp,yp] = polyfitpowers(lp)
    %P = [p00 p10 p01 p20 p11 p02 p30 p21 p12 p03...]
    n=(sqrt(1+8*lp)-3)/2;
    xp = zeros(1,lp);
    yp = zeros(1,lp);
    k = 1;
    for order = 0:n
        for j = 1:(order+1)
            xp(k) = order - j + 1;
            yp(k) = order - xp(k);
            k = k+1;
        end
    end
end
