function v = myvander(x,y,lp)
    v = zeros([length(x) lp]);
    [xp,yp] = polyfitpowers(lp);
    for j = 1:lp
        v(:,j) = x.^(xp(j)).*y.^(yp(j));
    end
end