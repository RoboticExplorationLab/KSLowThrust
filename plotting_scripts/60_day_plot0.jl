using JLD2, Attitude, MATLAB, SatelliteDynamics, Dierckx
# include ksfunctions.jl
include(joinpath(dirname(dirname(@__FILE__)),"ks_functions.jl"))
@load "rate_control_transfers/60_day_v2_transfer.jld2" X U t_hist
dscale = 1e7 # m
tscale = 20000 # s
uscale = 10000.0
dt = 3e-2

r_eci_hist = dscale*mat_from_vec([x_from_u(X[i][1:4]) for i = 1:length(X)])
v_eci_hist = (dscale/tscale)*mat_from_vec([xdot_from_u(X[i][1:4],X[i][5:8]) for i = 1:length(X)])

# TODO: Spline this in a data efficient way
mat"
figure
hold on
plot3($r_eci_hist(1,:),$r_eci_hist(2,:),$r_eci_hist(3,:))
axis equal
hold off
"

Xm = mat_from_vec(X)
U_real = Xm[12:14,:]

mat"
figure
hold on
title('Thrust Vector')
plot($t_hist,$U_real')
ylabel('Thrust Acceleration (0.0001m/s^2)')
xlabel('Time (days)')
hold off
xlim([0,$t_hist(end)])
legend('u_x','u_y','u_z')
%saveas(gcf,'long_traj.png')
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
