using JLD2
@load "rate_control_transfers/30_v2_day_transfer.jld2" X U
dscale = 1e7 # m
tscale = 20000 # s
uscale = 10000.0
dt = 3e-2
# get time stuff
function get_time_transfer(X,dt,tscale)
    """get the time in days for a run. This allows us to remove t from state"""
    # X = states(solver)
    t = 0.0
    t_hist = zeros(length(X))
    for i = 1:(length(X)-1)

        p0 = X[i][1:4]
        R0 = dot(p0,p0)

        p1 = X[i+1][1:4]
        R1 = dot(p1,p1)

        t += (R0 + R1)*dt/2
        t_hist[i+1] = copy(t)
    end
    t_days = t*tscale/(24*3600)
    t_hist = t_hist*tscale/(24*3600)
    return t_days, t_hist
end

t_days,t_hist = get_time_transfer(X,dt,tscale)
