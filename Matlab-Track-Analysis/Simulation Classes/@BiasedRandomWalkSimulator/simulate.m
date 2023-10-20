function simulate(s, varargin)
%function simulate(s, varargin)
%runs the BiasedRandomWalk simulation
%no optional args at this time

for j = 2:s.npts
    step = s.timestep * transpose([s.speed.*cos(s.vdir); s.speed.*sin(s.vdir)]);
    %size(step)
    %size(s.loc)
    %size(s.loc(:,:,j))
    s.loc(:,:,j) = s.loc(:,:,j-1) + step;
    
    updateVdir(s);
end

function updateVdir(s)
%function updateVdir(s)
s.vdir = s.vdir + s.driftCoeff * sqrt (s.timestep) * randn(size(s.vdir));
s.vdir = mod(s.vdir, 2*pi);
reoprob = s.timestep * interp1(s.thetaAxis, s.reoRate, s.vdir, 'linear');
reo = find(rand(size(reoprob)) < reoprob);

reoMag = interp1(s.thetaAxis, s.reoMag, s.vdir(reo), 'linear');
deltaTheta = randn(size(reoMag)).*reoMag;

%randomly flip based on s.reoBias, but only if flipping moves you closer to
%biasDirection
flip = rand(size(reoMag)) < interp1(s.thetaAxis, s.reoBias, s.vdir(reo), 'linear') & ...
    cos(s.vdir(reo) - deltaTheta - s.biasDirection) > cos(s.vdir(reo) + deltaTheta - s.biasDirection);
deltaTheta(flip) = -deltaTheta(flip);


s.vdir(reo) = s.vdir(reo) + deltaTheta;
s.vdir = mod(s.vdir, 2*pi);