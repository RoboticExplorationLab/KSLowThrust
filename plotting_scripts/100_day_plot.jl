using JLD2, Attitude, MATLAB, SatelliteDynamics, Dierckx
# include ksfunctions.jl
include(joinpath(dirname(dirname(@__FILE__)),"ks_functions.jl"))

@load "rate_control_transfers/112_day_transfer.jld2" X U

dscale = 1e7 # m
tscale = 20000 # s
uscale = 10000.0
dt = 4e-2


t_days,t_hist = get_time_transfer(X,dt,tscale)

r_eci_hist = dscale*mat_from_vec([x_from_u(X[i][1:4]) for i = 1:length(X)])
v_eci_hist = (dscale/tscale)*mat_from_vec([xdot_from_u(X[i][1:4],X[i][5:8]) for i = 1:length(X)])

# TODO: Spline this in a data efficient way
new_t = discretize_t(t_hist,4)
mat"
$rnew = spline($t_hist,$r_eci_hist,$new_t)/1000;
"
mat"
rnew = spline($t_hist,$r_eci_hist,$new_t)/1000;
figure('Renderer', 'painters', 'Position', [10 10 1300 300])
hold on
[x,y,z] = sphere(40);
imgRGB = imread('/Users/kevintracy/devel/Low-Thrust-TrajOpt/matlab_video/earth.jpg');
warp(6378*x,6378*y,6378*z,rot90(imgRGB,2))
plot3(rnew(1,:),rnew(2,:),rnew(3,:),'linewidth',0.4)
xlabel('ECI X (km)','Interpreter','Latex')
ylabel('ECI Y (km)','Interpreter','Latex')
zlabel('ECI Z (km)','Interpreter','Latex')
set(gca,'TickLabelInterpreter','latex')
set(gca,'Fontsize',28)
%grid on
axis equal
hold off
view(100,5)
saveas(gcf,'test.eps','epsc')
saveas(gcf,'test.png')
"


# top view
mat"
rnew = spline($t_hist,$r_eci_hist,$new_t)/1000;
figure('Renderer', 'painters', 'Position', [10 10 900 900])
hold on
[x,y,z] = sphere(40);
imgRGB = imread('/Users/kevintracy/devel/Low-Thrust-TrajOpt/matlab_video/earth.jpg');
warp(6378*x,6378*y,6378*z,rot90(imgRGB,2))
plot3(rnew(1,:),rnew(2,:),rnew(3,:),'linewidth',1.0)
xlabel('ECI X (km)','Interpreter','Latex')
ylabel('ECI Y (km)','Interpreter','Latex')
zlabel('ECI Z (km)','Interpreter','Latex')
set(gca,'TickLabelInterpreter','latex')
set(gca,'Fontsize',28)
%grid on
axis equal
hold off
view(0,90)
saveas(gcf,'topview.eps','epsc')
saveas(gcf,'topview.png')
"

mat"
close all
close all hidden
rnew = spline($t_hist,$r_eci_hist,$new_t)/1000;
rnew = round(rnew,4,'significant');
figure
hold on
plot(rnew(1,:),rnew(2,:))
xlabel('ECI X (km)')
ylabel('ECI Y (km)')
axis equal
hold off
addpath('/Users/kevintracy/devel/Low-Thrust-TrajOpt/matlab2tikz/src')
%matlab2tikz('topview.txt')
"

mat"
close all
close all hidden
rnew = spline($t_hist,$r_eci_hist,$new_t)/1000;
rnew = round(rnew,4,'significant');
figure
hold on
plot(rnew(1,:),rnew(3,:))
xlabel('ECI X (km)')
ylabel('ECI Y (km)')
axis equal
hold off
%addpath('/Users/kevintracy/devel/Low-Thrust-TrajOpt/matlab2tikz/src')
%matlab2tikz('sideview.txt')
"

# mat"
# close all
# close all hidden
# [X,Y,Z] = sphere(20)
# Re = 6010
# rnew = spline($t_hist,$r_eci_hist,$new_t)/1000;
# rnew = round(rnew,4,'significant');
# figure
# hold on
# plot3(rnew(1,:),rnew(2,:),rnew(3,:),'k')
# h = surf(Re*X,Re*Y,Re*Z)
# axis equal
# view(-69,14)
# hold off
# h.EdgeAlpha = 0
# h.FaceColor = 'blue'
# %get(h)
# addpath('/Users/kevintracy/devel/Low-Thrust-TrajOpt/matlab2tikz/src')
# %matlab2tikz('globe.txt')
# "
# mat"
# close all
# close all hidden
# [X,Y,Z] = sphere(20)
# Re = 6010
# rnew = spline($t_hist,$r_eci_hist,$new_t)/1000;
# rnew = round(rnew,4,'significant');
# figure
# hold on
# %plot3(rnew(1,:),rnew(2,:),rnew(3,:),'k')
# plot(rnew(2,:),rnew(3,:),'k')
# %h = surf(Re*X,Re*Y,Re*Z)
# axis equal
# hold off
# xlabel('ECI X (km)')
# ylabel('ECI Y (km)')
# zlabel('ECI Z (km)')
# %get(h)
# addpath('/Users/kevintracy/devel/Low-Thrust-TrajOpt/matlab2tikz/src')
# %matlab2tikz('globe.txt')
# "
Xm = mat_from_vec(X)
U_real = Xm[12:14,:]

mat"
figure
hold on
title('Thrust Vector')
plot($t_hist,$U_real*10')
ylabel('Thrust Magnitude')
xlabel('Time (days)')
legend('u_x','u_y','u_z','Location','SouthWest')
hold off
xlim([0,$t_hist(end)])
a2 = axes()
a2.Position = [0.6500 0.750 0.2 0.15]; % xlocation, ylocation, xsize, ysize
plot(a2,$t_hist(1:50),10*$U_real(:,1:50)'); axis tight

%saveas(gcf,'long_traj.png')
"

mat"
figure
hold on
subplot(2,2,1:2)
plot($t_hist,$U_real*0.1')
ylabel('Thrust Magnitude (N)')
xlabel('Time (days)')
legend('u_x','u_y','u_z','Location','NorthEast')
hold off
xlim([0,$t_hist(end)])

subplot(2,2,3)
plot($t_hist(1:50),$U_real(:,1:50)*0.1','linewidth',1.2)
ylabel('Thrust Magnitude (N)')
xlabel('Time (days)')

subplot(2,2,4)
plot($t_hist(2000:2050),$U_real(:,2000:2050)*0.1','linewidth',1.2)
xlim([$t_hist(2000),$t_hist(2050)])
ylabel('Thrust Magnitude (N)')
xlabel('Time (days)')
%saveas(gcf,'long_traj.png')
addpath('/Users/kevintracy/devel/Low-Thrust-TrajOpt/matlab2tikz/src')
%matlab2tikz('thrust_100.txt')
"

# get this in RTN

# function thrust_angle_from_rtn(u_rtn)
#     u_rtn = normalize(u_rtn)
#     α = atan(u_rtn[1],u_rtn[2])
#     β = atan(u_rtn[3],norm(u_rtn[1:2]))
#     return [α; β]
# end
# function get_rtn_in_out(U_real,r_eci_hist,v_eci_hist)
#
#     return [thrust_angle_from_rtn(
#                 rECItoRTN([r_eci_hist[:,i];v_eci_hist[:,i]])*U_real[:,i] ) for i = 1:size(U_real,2)]
# end

angle_hist = get_rtn_in_out(U_real,r_eci_hist,v_eci_hist)
angles = mat_from_vec(angle_hist)

mat"
figure
hold on
plot($t_hist,rad2deg($angles'))
legend('In Plane','Out of Plane')
hold off
"
mat"
angles = rad2deg($angles)
figure
hold on
subplot(2,2,1:2)
plot($t_hist,angles')
ylabel('Thrust Angle (deg)')
xlabel('Time (days)')
legend('In Plane','Out of Plane','Location','SouthWest')
hold off
xlim([0,$t_hist(end)])

subplot(2,2,3)
plot($t_hist(1:50),angles(:,1:50)','linewidth',1.2)
ylabel('Thrust Magnitude (N)')
xlabel('Time (days)')

subplot(2,2,4)
plot($t_hist(2000:2050),angles(:,2000:2050)','linewidth',1.2)
xlim([$t_hist(2000),$t_hist(2050)])
ylabel('Thrust Magnitude (N)')
xlabel('Time (days)')
%saveas(gcf,'long_traj.png')
addpath('/Users/kevintracy/devel/Low-Thrust-TrajOpt/matlab2tikz/src')
%matlab2tikz('thrust_100_angles.txt')
"


U_real_vec = [X[i][12:14] for i = 1:length(X)]
Unorm = [norm(U_real_vec[i]) for i = 1:length(X)]

mat"
figure
hold on
plot($t_hist,$Unorm/$uscale)
xlim([$t_hist(4),$t_hist(end)])
xlabel('Time (days)')
ylabel('Thrust Magnitude (m/s^2)')
hold off
"


println("Time in days:  $t_days")


# now let's get orbital elements
const SD = SatelliteDynamics
oe_hist = [SD.sCARTtoOSC([r_eci_hist[:,i];v_eci_hist[:,i]]) for i = 1:size(r_eci_hist,2)]

a_hist = [oe_hist[i][1] for i = 1:length(oe_hist)]
e_hist = [oe_hist[i][2] for i = 1:length(oe_hist)]
i_hist = [rad2deg(oe_hist[i][3]) for i = 1:length(oe_hist)]

mat"
figure
hold on
subplot(3,1,1)
hold on
plot($t_hist,$a_hist/1000)
h = plot($t_hist, 42e3*ones(length($t_hist),1),'r');
legend([h(1)],'GEO SMA','Location','SouthEast')
ylabel('Semi-Major Axis (km)')
xlabel('Time (days)')
xlim([0,$t_hist(end)])
ylim([20e3,46e3])
subplot(3,1,2)
hold on
plot($t_hist,$e_hist)
ylabel('Eccentricity')
xlabel('Time (days)')
xlim([0,$t_hist(end)])
subplot(3,1,3)
hold on
plot($t_hist,$i_hist)
ylabel('Inclination (deg)')
xlabel('Time (days)')
xlim([0,$t_hist(end)])
hold off
"

# matlab animation stuff
output_to_matlab = (new_time_vec = new_t, x_hist_t_correct = rnew)
file = matopen("matlab_video_data_100.mat", "w")
write(file, "output_to_matlab", output_to_matlab)
close(file)

hundred = (t_hist = t_hist, a_hist = a_hist, e_hist = e_hist, i_hist = i_hist, Unorm = Unorm,angles = angles)

using JLD2
@save "plotting_data/plotting_data_100.jld2" hundred
