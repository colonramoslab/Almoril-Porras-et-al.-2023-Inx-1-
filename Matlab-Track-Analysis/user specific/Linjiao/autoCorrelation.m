cryo=eset;
[ac, np, tx] = cryo.autocorrelate('vnorm');
clf
semilogy(tx(ac>0), ac(ac>0)./np(ac>0), 'b.');
xlim([0 600]);
xlabel ('$\tau$ (s)','Interpreter', 'Latex');
ylabel('$\langle\hat{v}(t)\cdot\hat{v}(t + \tau)\rangle$','Interpreter','Latex');
title ('Auto-Correlation of velocity direction');
ylim([0.01 1]);
snapnow;
