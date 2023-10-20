function directionPlots(track)
%function directionPlots(track)
eti = track.getDerivedQuantity('eti');
th = unwrap(track.getDerivedQuantity('theta'));
mh = track.getDerivedQuantity('smhdir');
tm = track.getDerivedQuantity('stmdir');

hth = atan2(mh(2,:), mh(1,:));
tth = atan2(tm(2,:), tm(1,:));

bob = unwrap([th;hth;tth], [], 1);
hth = bob(2,:);
tth = bob(3,:);

plot (eti, th, eti, hth, eti, tth); legend ('vel dir', 'head dir', 'tail dir');