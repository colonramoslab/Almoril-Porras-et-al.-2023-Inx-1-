function draw(reo, varargin)
%MaggotReorientation.draw(varargin)
%function draw(reo, varargin)


ih = ishold;

center = mean(reo.track.dq.sloc(:,reo.inds),2);
allpts = [reo.track.dq.sloc(:,reo.inds) reo.track.dq.stail(:,reo.inds) reo.track.dq.shead(:,reo.inds)];
radius = max(sqrt (sum( (allpts - repmat(center,1,length(allpts))).^2)));

t = [0:0.1:2*pi]; plot (center(1) + radius * cos(t), center(2)+radius*sin(t), 'k--'); hold on

for j = 1:length(reo.headSwing)
    reo.headSwing(j).draw();
end

if (~ih)
    hold off
end