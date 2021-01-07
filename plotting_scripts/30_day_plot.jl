using JLD2, Attitude, MATLAB, SatelliteDynamics, Dierckx
# include ksfunctions.jl
include(joinpath(dirname(dirname(@__FILE__)),"ks_functions.jl"))
@load "rate_control_transfers/30_v3_day_transfer.jld2" X U
dscale = 1e7 # m
tscale = 20000 # s
uscale = 10000.0
dt = 3e-2


t_days,t_hist = get_time_transfer(X,dt,tscale)

r_eci_hist = dscale*mat_from_vec([x_from_u(X[i][1:4]) for i = 1:length(X)])
v_eci_hist = (dscale/tscale)*mat_from_vec([xdot_from_u(X[i][1:4],X[i][5:8]) for i = 1:length(X)])

# function discretize_t(t_hist,n)
#     # new time vec
#     new_t = [0.0]
#     for i = 1:(length(t_hist)-1)
#
#         # time spacing betwen two knot points
#         δt = t_hist[i+1] - t_hist[i]
#
#         # add intermediate points between the knot points
#         for j = 1:n
#             push!(new_t,δt*j/n + t_hist[i])
#         end
#     end
#     return new_t
# end
new_t = discretize_t(t_hist,4)
mat"
rnew = spline($t_hist,$r_eci_hist,$new_t);
figure
hold on
plot3(rnew(1,:),rnew(2,:),rnew(3,:))
%axis equal
hold off
"
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

angle_hist = get_rtn_in_out(U_real,r_eci_hist,v_eci_hist)

angles = mat_from_vec(angle_hist)

mat"
figure
hold on
plot(round($t_hist,4,'Significant'),round(rad2deg($angles'),4,'Significant'),'linewidth',1.2)
xlabel('Time (days)')
xlim([0 $t_hist(end)])
ylabel('Angle (degrees)')
legend('In Plane','Out of Plane','Location','southwest')
hold off
addpath('/Users/kevintracy/devel/Low-Thrust-TrajOpt/matlab2tikz/src')
%matlab2tikz('thrust_angles_30.tex')
"

# picture in picture
mat"
figure
hold on
plot(round($t_hist,4,'Significant'),round(rad2deg($angles'),4,'Significant'),'linewidth',1.2)
xlabel('Time (days)')
xlim([0 $t_hist(end)])
ylabel('Angle (degrees)')
legend('In Plane','Out of Plane','Location','southwest')
a2 = axes()
a2.Position = [0.1200 0.2600 0.2 0.2]; % xlocation, ylocation, xsize, ysize
plot(a2,$t_hist(1:500),$angles(:,1:500)'); axis tight
hold off
addpath('/Users/kevintracy/devel/Low-Thrust-TrajOpt/matlab2tikz/src')
%matlab2tikz('thrust_angles_30.tex')
"


# now let's get angle from velocity vector
angle_from_v = [acos(dot(normalize(v_eci_hist[:,i]),normalize(U_real[:,i]))) for i = 1:size(U_real,2)]

mat"
figure
hold on
plot(rad2deg($angle_from_v))
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

thirty = (t_hist = t_hist, a_hist = a_hist, e_hist = e_hist, i_hist = i_hist, Unorm = Unorm)

using JLD2
@save "plotting_data_30.jld2" thirty
