function [meandt, meandtsq, pleft, meandt_ci, meandtsq_ci] = skew_normal_model_reorientationMeansVsTheta (params, thetaIn, hess, CL)
% function [meandt, meandtsq, pleft, meandt_ci, meandtsq_ci] =
% skew_normal_model_reorientationMeansVsTheta (params, thetaIn, hess, CL)
% P(dt | theta) = (0.5 + g) * S(U, sigma, skew, dt) + (0.5 - g) * S(U,
% sigma, skew, -dt)
%
% U = A - B*cos(theta - theta_0);
% skew = D - E*cos(theta - theta_0);
% g = C*(sin(theta-theta_0))
%other fit parameters, sigma, theta_0 
%
%model.params = [A, B, C, D, E, sigma, theta_0];theta_0 = x(8);


    x = params;
    theta_0 = x(7,:);
    if (size(thetaIn,2) > 1)
        thetaIn = thetaIn';
    end
    nti = length(thetaIn);
    tt = repmat(thetaIn,1,length(theta_0)) - repmat(theta_0,nti, 1);

    g = repmat(x(3,:),nti,1).*sin(tt);
    u = repmat(x(1,:),nti,1) - repmat(x(2,:),nti,1).*cos(tt);
    skew = repmat(x(4,:),nti,1) - repmat(x(5,:),nti,1).*cos(tt);
    sigma = repmat(x(6,:),nti,1);

    [mu,sig] = skewNormalMean(u(:), sigma(:), skew(:));
    mu = reshape(mu, size(tt));
    sig = reshape(sig, size(tt));
    meandt = -2*g.*mu;
    meandtsq = (mu.^2 + sig.^2);

    %sigma = sqrt(meandtsq - meandt.^2);
    pleft = 0.5 - g;

    if (nargout > 3 && nargin > 2)
        [~,meandt_ci] = reorientationMeanVsTheta(params, thetaIn, hess, CL);
    end
    if (nargout > 4 && nargin > 2)
        [~,meandtsq_ci] = reorientationMeanSqVsTheta(params, thetaIn, hess, CL);
    end
end

function [meandt, meandt_ci] = reorientationMeanVsTheta (params, thetaIn, hess, confidenceLevel)
    x = params;
    theta_0 = x(7,:);
    tt = thetaIn - theta_0;

    g = x(3,:).*sin(tt);
    u = x(1,:) - x(2,:).*cos(tt);
    skew = x(4,:) - x(5,:).*cos(tt);
    sigma = x(6,:);

    mu = skewNormalMean(u, sigma, skew);
    meandt = -2*g.*mu;

    if (nargout > 1 && nargin > 2)
        meandt_ci = zeros(2, length(thetaIn), length(confidenceLevel));
        for j = 1:length(thetaIn);
            fun = @(x) reorientationMeanVsTheta (x, thetaIn(j));
            [minv maxv] = confidenceRange(x, hess, fun, confidenceLevel, 'numtrials', 1E5);
            meandt_ci(1,j,:) = minv;
            meandt_ci(2,j,:) = maxv;
        end
    end

end


function [meandtsq, meandtsq_ci] = reorientationMeanSqVsTheta (params, thetaIn, hess, confidenceLevel)
    x = params;
    theta_0 = x(7,:);
    tt = thetaIn - theta_0;

    %g = x(3,:).*sin(tt);
    u = x(1,:) - x(2,:).*cos(tt);
    skew = x(4,:) - x(5,:).*cos(tt);
    sigma = x(6,:);

    [mu,sig] = skewNormalMean(u, sigma, skew);
    meandtsq = (mu.^2 + sig.^2);

    if (nargout > 1 && nargin > 2)
        meandtsq_ci = zeros(2, length(thetaIn), length(confidenceLevel));
        for j = 1:length(thetaIn);
            fun = @(x) reorientationMeanSqVsTheta (x, thetaIn(j));
            [minv maxv] = confidenceRange(x, hess, fun, confidenceLevel, 'numtrials', 1E5);
            meandtsq_ci(1,j,:) = minv;
            meandtsq_ci(2,j,:) = maxv;
        end
    end

end


