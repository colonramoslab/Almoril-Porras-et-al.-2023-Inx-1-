function compareWithExtraction(track, ind)

pt = track.pt(ind);

subplot(1,2,1); pt.drawTrackImage(track.expt.camcalinfo); colormap jet; shading interp

xl = get(gca, 'XLim');
yl = get(gca, 'Ylim');
colormap jet; shading interp
hold(gca, 'on');
track.plotPath('sloc', 'b.-','highlightinds', ind);
xlim(xl);
ylim(yl);
hold(gca, 'off');

itp = ImTrackPoint(pt);
subplot(1,2,2); itp.drawTrackImage(track.expt.camcalinfo); 

xl = get(gca, 'XLim');
yl = get(gca, 'Ylim');
colormap jet; shading interp
hold(gca, 'on');
track.plotPath('sloc', 'b.-','highlightinds', ind);
xlim(xl);
ylim(yl);
hold(gca, 'off');