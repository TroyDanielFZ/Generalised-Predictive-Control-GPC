% clf;
% subplot(141);
% plot(dat1.time, dat1.signals(1).values, '--r');
% hold on;
% plot(dat1.time, dat1.signals(2).values, '--.g');
% plot(dat1.time, dat1.signals(3).values, '-b');
% plot(xlim, [0,0], '--.c');
% title('Case 1');

% subplot(142);
% plot(dat2.time, dat2.signals(1).values, '--r');
% hold on;
% plot(dat2.time, dat2.signals(2).values, '--.g');
% plot(dat2.time, dat2.signals(3).values, '-b');
% plot(xlim, [0,0], '--.c');
% title('Case 2');

% subplot(143);
% plot(dat3.time, dat3.signals(1).values, '--r');
% hold on;
% plot(dat3.time, dat3.signals(2).values, '--.g');
% plot(dat3.time, dat3.signals(3).values, '-b');
% plot(xlim, [0,0], '--.c');
% title('Case 3');

% subplot(144);
% plot(dat4.time, dat4.signals(1).values, '--r');
% hold on;
% plot(dat4.time, dat4.signals(2).values, '--.g');
% plot(dat4.time, dat4.signals(3).values, '-b');
% plot(xlim, [0,0], '--.c');
% title('Case 4');

% h = legend('$y_r$', '$y$', '$e$', 'location', 'east');
% % h = legend('Adaptive','PID','$x_d$');
% set(h,'FontSize',14, 'interpreter', 'latex')

% ----------------------------------------------------------------------------------
% Hysterisis + Delay
% ----------------------------------------------------------------------------------

fig = figure(1); clf;
fig.Name = 'GPC';
subplot(121);
plot(dat.time + 0.6, dat.signals(1).values, '--r');
hold on;
plot(dat.time, dat.signals(2).values, '-.c');
plot(dat.time, dat.signals(3).values, '-b');
h = legend('$y_r$', '$e$', '$y$');
set(h,'FontSize',14, 'interpreter', 'latex')
xlim([0, max(dat.time)]);
ylim([-15, 15]);

subplot(122);
plot(-dat.signals(4).values(151:end), dat.signals(3).values(151:end))

fig = figure(2); clf;
fig.Name = 'PID';
subplot(121);
plot(dat1.time, dat1.signals(1).values, '--r');
hold on;
plot(dat1.time, dat1.signals(2).values, '-.c');
plot(dat1.time, dat1.signals(3).values, '-b');
h = legend('$y_r$', '$e$', '$y$');
set(h,'FontSize',14, 'interpreter', 'latex')

subplot(122);
plot(-dat1.signals(4).values(151:end), dat1.signals(3).values(151:end))
	
