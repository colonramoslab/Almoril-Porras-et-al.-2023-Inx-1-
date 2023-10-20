% copied manually from igor
L3Atemps= [17.3 19.7 22.4];
L3AturnHeat= [2.4 2.71 3.36];
L3AturnCool= [4.03 2.62 2.47];
L3AthetaHeat= [53.2 62.1 65.1];
L3AthetaCool= [58.4 65.6 59.3];
L3AhsHeat= [18.9 13.5 18.5];
L3AhsCool= [24.9 11.2 13.2];

colordef black;
figure;
plot(L3Atemps, L3AturnHeat,'r.-', 'LineWidth', 2, 'MarkerSize', 22);
hold on;
plot(L3Atemps, L3AturnCool,'b.-', 'LineWidth', 2, 'Markersize', 22);
hold off;
% axis labels
xlabel('Temperature (C)');
ylabel('Turning rate (turns/min)');
% legend
legh= legend('Heating','Cooling');
setBlack1Params;
set(legh,'Position',[0.6525    0.7474    0.2170    0.1313]);
% turn this on at the end to get anti-aliasing!
myaa;

figure;
plot(L3Atemps, L3AthetaHeat,'r.-', 'LineWidth', 2, 'MarkerSize', 22);
hold on;
plot(L3Atemps, L3AthetaCool,'b.-', 'LineWidth', 2, 'Markersize', 22);
hold off;
% axis labels
xlabel('Temperature (C)');
ylabel('Reorientations with > 1 head sweep (%)');
% legend
legh= legend('Heating','Cooling');
setBlack1Params;
set(legh,'Position',[0.6547    0.1697    0.2170    0.1313]);
% turn this on at the end to get anti-aliasing!
myaa;

figure;
plot(L3Atemps, L3AhsHeat,'r.-', 'LineWidth', 2, 'MarkerSize', 22);
hold on;
plot(L3Atemps, L3AhsCool,'b.-', 'LineWidth', 2, 'Markersize', 22);
hold off;
% axis labels
xlabel('Temperature (C)');
ylabel('Turning rate index');
% legend
legh= legend('Heating','Cooling');
setBlack1Params;
set(legh,'Position',[0.6525    0.7474    0.2170    0.1313]);
% turn this on at the end to get anti-aliasing!
myaa;
