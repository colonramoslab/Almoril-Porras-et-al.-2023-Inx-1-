function [ang_new] = AngleAdd (ang1,ang2)

% if indeg==1
%     ang1=deg2rad(ang1);
%     ang2=deg2rad(ang2);
% end

%add angles
ang_new=ang1+ang2;

%now make sure angle is between -pi and pi
Z = exp(i*ang_new);
ang_new=angle(Z);



end