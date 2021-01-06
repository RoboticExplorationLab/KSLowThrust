

using JLD2
@load "plotting_data_30.jld2" thirty
@load "plotting_data_60.jld2" sixty
@load "plotting_data_100.jld2" hundred

mat"
idx = 5
%(1:idx:end)
figure
hold on
subplot(3,1,1)
hold on
h(4)=plot($thirty.t_hist(1:idx:end),$thirty.a_hist(1:idx:end)/1000,'linewidth',2)
h(3)=plot($sixty.t_hist(1:idx:end),$sixty.a_hist(1:idx:end)/1000,'linewidth',2)
h(2)=plot($hundred.t_hist(1:idx:end),$hundred.a_hist(1:idx:end)/1000,'linewidth',2)
h(1) = plot([0 $hundred.t_hist(end)], 42e3*ones(2,1),'Color','#9e9e9e','linewidth',2);
legend([h(1)],'GEO SMA','Location','SouthEast')
ylabel('Semi-Major Axis (km)')
xlabel('Time (days)')
xlim([0,$hundred.t_hist(end)])
ylim([20e3,46e3])
subplot(3,1,2)
hold on
plot($thirty.t_hist(1:idx:end),$thirty.e_hist(1:idx:end),'linewidth',2)
plot($sixty.t_hist(1:idx:end),$sixty.e_hist(1:idx:end),'linewidth',2)
plot($hundred.t_hist(1:idx:end),$hundred.e_hist(1:idx:end),'linewidth',2)
ylabel('Eccentricity')
xlabel('Time (days)')
xlim([0,$hundred.t_hist(end)])
subplot(3,1,3)
hold on
plot($thirty.t_hist(1:idx:end),$thirty.i_hist(1:idx:end),'linewidth',2)
plot($sixty.t_hist(1:idx:end),$sixty.i_hist(1:idx:end),'linewidth',2)
plot($hundred.t_hist(1:idx:end),$hundred.i_hist(1:idx:end),'linewidth',2)
legend('30 day','60 day','100 day','location','northeast')
ylabel('Inclination (deg)')
xlabel('Time (days)')
xlim([0,$hundred.t_hist(end)])
hold off
%addpath('matlab2tikz')
addpath('/Users/kevintracy/devel/Low-Thrust-TrajOpt/matlab2tikz/src')
disp('pre')
matlab2tikz('mytest.tikz')
disp('post')
"

# mat"
# figure
# hold on
# subplot(2,2,1)
# hold on
# plot($thirty.t_hist,$thirty.a_hist/1000)
# plot($sixty.t_hist,$sixty.a_hist/1000)
# plot($hundred.t_hist,$hundred.a_hist/1000)
# h(1) = plot($hundred.t_hist, 42e3*ones(length($hundred.t_hist),1),'k--')
# legend([h(1)],'GEO SMA','Location','SouthEast')
# ylabel('Semi-Major Axis (km)')
# xlabel('Time (days)')
# xlim([0,$hundred.t_hist(end)])
# ylim([20e3,46e3])
# subplot(2,2,2)
# hold on
# plot($thirty.t_hist,$thirty.e_hist)
# plot($sixty.t_hist,$sixty.e_hist)
# plot($hundred.t_hist,$hundred.e_hist)
# ylabel('Eccentricity')
# xlabel('Time (days)')
# xlim([0,$hundred.t_hist(end)])
# subplot(2,2,3)
# hold on
# h1=plot($thirty.t_hist,$thirty.i_hist)
# h2=plot($sixty.t_hist,$sixty.i_hist)
# h3=plot($hundred.t_hist,$hundred.i_hist)
# %legend('30 day','60 day','100 day','location','southoutside')
# ylabel('Inclination (deg)')
# xlabel('Time (days)')
# xlim([0,$hundred.t_hist(end)])
# hold off
#
# hL = subplot(2,2,4);
# poshL = get(hL,'position');     % Getting its position
#
# lgd = legend(hL,[h1;h2;h3],'RandomPlot1','RandomPlot2','RandomPlot3');
# set(lgd,'position',poshL);      % Adjusting legend's position
# axis(hL,'off');
# set(hL,'position',[0 0 0.1 0.1])
# saveas(gcf,'test.png')
# "
