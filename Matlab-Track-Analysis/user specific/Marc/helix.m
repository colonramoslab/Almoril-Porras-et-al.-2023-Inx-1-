z = 1:100;
theta = 4*pi * z / max(z);
r_inner = 1;
r_outer = 2;

xi = r_inner * cos(theta);
yi = r_inner * sin(theta);
xo = r_outer * cos(theta);
yo = r_outer * sin(theta);

plot3 (xi,yi,z, xo, yo, z);