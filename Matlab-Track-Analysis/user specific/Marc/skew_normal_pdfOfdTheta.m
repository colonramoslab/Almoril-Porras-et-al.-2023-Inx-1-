function p = skew_normal_pdfOfdTheta (x,thetaIn, dTheta)
%function p = skew_normal_pdfOfdTheta (x,thetaIn, dTheta)


[r,l] =  skew_normal_pdfOfdThetaBreakdown (x,thetaIn, dTheta);
p = r + l;
[r2,l2] = skew_normal_pdfOfdThetaBreakdown(x, thetaIn, 2*pi+dTheta);
p = p + r2 + l2;
[r2,l2] = skew_normal_pdfOfdThetaBreakdown(x, thetaIn, -2*pi+dTheta);
p = p + r2 + l2;

if (any (p < 0))
    %disp (['negative probability from x = ' num2str(x','%10.5e\t')]);
  %  pause
end
if (any (~isfinite(p)))
   % disp (['non-finite probability from x = ' num2str(x','%10.5e\t')]);
   % pause
end
p(~isfinite(p)) = eps;
p(p < 0) = eps;
end