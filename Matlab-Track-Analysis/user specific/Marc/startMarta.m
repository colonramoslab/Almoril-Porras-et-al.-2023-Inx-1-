close all
clear classes
bd = 'C:\mae\marta analysis explanations\ExampleData\raw data for three animals\time series of spine skeleton etc\';
fn = '20100122_100901@w1118@n@t1@b_a1v_15s1x15s0s#n#n#n@1@.00001.spine';
mt = MartaTrack.fromFile([bd fn]);
mt.so.stop_speed_cut = 0.5;
mt.so.start_speed_cut = 0.575;
mt.segmentTrack;
mt.plotSegmentation;