function [] = SaveImages()

dir_name = uigetdir;

figure(1);
saveas(gca,[dir_name '\run_duration_vs_heading.jpg']);

figure(2);
saveas(gca,[dir_name '\run_speed.jpg']);

figure(3);
saveas(gca,[dir_name '\run_duration_hist.jpg']);

figure(4);
saveas(gca,[dir_name '\run_heading_hist.jpg']);

figure(5);
saveas(gca,[dir_name '\run_steering.jpg']);

figure(6);
saveas(gca,[dir_name '\run_steering_means.jpg']);

figure(7);
saveas(gca,[dir_name '\turn_steering.jpg']);

figure(8);
saveas(gca,[dir_name '\turn_steering_means.jpg']);

figure(10);
saveas(gca,[dir_name '\turn_HS_count.jpg']);

figure(11);
saveas(gca,[dir_name '\turn_HS_steering.jpg']);

figure(12);
saveas(gca,[dir_name '\turn_HS_steering_means.jpg']);

figure(15);
saveas(gca,[dir_name '\turn_HS_initbias.jpg']);

figure(16);
saveas(gca,[dir_name '\turn_HS_acceptbias.jpg']);




end