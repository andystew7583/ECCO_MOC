%%%
%%% plotTaabw.m
%%% 
%%% Example script to plot time series of diagnosed and predicted AABw
%%% transport.
%%%

load Taabw.mat;

%%% Apply 30-day smoothing
Taabw_anom_diagnosed_smooth = smooth(Taabw_anom_diagnosed,30);
Taabw_anom_wind_smooth = smooth(Taabw_anom_wind,30);

figure(1);
plot(tt,Taabw_anom_diagnosed_smooth/1e6);
hold on;
plot(tt,Taabw_anom_wind_smooth/1e6);
hold off;
set(gca,'FontSize',14);
datetick('x')
yhandle = ylabel('AABW flux anomaly $T^\prime_\mathrm{AABW}$ (Sv)','interpreter','latex');
set(gca,'YLim',[-25 25]);
leghandle = legend('Diagnosed','Reconstructed from wind stress','interpreter','latex','Location','NorthEast');
set(leghandle,'Orientation','horizontal');