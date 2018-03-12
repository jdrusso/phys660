#!/usr/bin/python
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

END = -1
COLLISION = 0
BOUNCE = 1

x1, x2 = 1.0, 3.0
v1, v2 = 0.0, 0.0

m1, m2 = 1., 1.
m = m1 + m2

g = -1

Ei = m1*x1 + m2*x2 + 0.5*m1*v1**2 + 0.5*m2*v2**2

### Normalize Xs and Vs
x1 /= Ei/m
x2 /= Ei/m

v1 /= np.sqrt(Ei/m)
v2 /= np.sqrt(Ei/m)


E = 0.5*(m1/m)* v1**2 + 0.5*(m2/m)*v2**2 + (m1/m)*x1 + (m2/m)*x2
print("Initial energy is %f, normalized to %f" % (Ei, E))


tmax = 2000
dt = 1.E-2
steps = tmax/dt


# Normalize positions and velocities


x = np.zeros([2,1])
v = np.zeros([2,1])
t = np.zeros([1])

x[0,0], x[1,0] = x1, x2
v[0,0], v[1,0] = v1, v2


# General algorithm:
#   1. Determine whether next event will be a bounce or a collision.
#   2. Determine the time, velocities and positions of that event.
#   3. If time > t_max, then stop.
#   4. Using the time, positions, and velocities of the balls at that event,
#       calculate the time of the next event.

# Find what the next event will be for a given set of intial conditions, and in
#   how much time.
def next_event(xs, vs, t):

    v1, v2 = vs[:,-1]
    x1, x2 = xs[:,-1]

    # Do this by comparing bounce times, using the initial conditions from the
    #   last elements in v and x.
    bounce_time1 = [x for x in np.roots([0.5*g, v1, x1]) if x > 0][0]
    bounce_time2 = [x for x in np.roots([0.5*g, v2, x2]) if x > 0][0]

    # If the second ball would bounce before the first, then a collision will occur
    collide = bounce_time2 < bounce_time1

    # Find the time of the event
    event_time = (x1-x2) / (v2-v1) if collide else bounce_time1

    # Update v and x using kinematics equations
    #   I.e., find where the balls are at the time of the event
    x1 += v1*event_time + 0.5 * g * event_time**2
    v1 += g*event_time

    x2 += v2*event_time + 0.5 * g * event_time**2
    v2 += g*event_time


    # If it's a collision, also reverse the top velocity
    if collide:
        _v2 = (2*m1)/(m1+m2)*v1 - (m1-m2)/(m1+m2)*v2
        _v1 = (m1-m2)/(m1+m2)*v1 + (2*m2)/(m1+m2)*v2

        v1 = _v1
        v2 = _v2

        #x2 += .0000001

    # If it's a bounce, set x1 to absolute value so it doesn't go through the
    #   floor.
    elif not collide:
        x1 = abs(x1)
        v1 *= (-1)


    # Update the entries in the list
    xs = np.append(xs, [[x1],[x2]], axis=1)
    vs = np.append(vs, [[v1],[v2]], axis=1)
    t = np.append(t, t[-1]+event_time)

    #  Return the updated lists
    return xs, vs, t


# Find all the events, record the times, positions, and velocities at them,
#   and then interpolate functions between.
while t[-1] <= tmax:

    # Come up with the list of events
    x, v, t = next_event(x, v, t)

print("Found %d events" % (len(t)-1))


##################### Poincare section ##################
# At this point, there's a list of positions, velocities, and times at all the
#   collisions. This is what's needed for part 3, plotting the Poincare section
#   of v2 vs x2 at collision times.

print("Generating Poincare sections")

# For the Poincare section, I want to filter out these events to ONLY collisions.
colls = x[0] == x[1]
print("Found %d collisions" % len([c for c in colls if c == True]))


plt.subplot(311)
# Plot only the collision elements by indexing the arrays with a Boolean mask
plt.plot(v[1, colls], x[1, colls], '.')
plt.xlabel("v_2")
plt.ylabel("x_2")
#########################################################




print("Generating trajectory functions")
# Create functions describing trajectories between events
# Build a list of each trajectory function, and the domain over
#   which it applies.
traj1s = []
traj2s = []

for i in range(1,len(t)):
    start = t[i-1]
    end = t[i]

    # These lambda functions calculate the position and velocity of the ball at
    #   some time. This is like a piecewise trajectory function, where each
    #   trajectory is defined within some domain by the initial conditions at
    #   the start of that.
    # The alternative to this would be storing the initial conditions over each
    #   interval and then just calling a trajectory function and passing it those,
    #   but I started with this because I wanted to use np.piecewise() to evaluate
    #   all the times.
    # However, it turned out that doing it that way I had to build a matrix of
    #   booleans that would define at what time each part of the piecewise
    #   function is called. That matrix was of dimension
    #   (number of sampled times * number of events), which rapidly got unmanageably
    #   large. For example, at t=20,000 that was using more than 16GB of memory.
    # As some bonus content on how Python handles lambda functions, variables
    #   in a lambda function are evaluated when the function is called, not when
    #   it's created. Here, that'll cause issues, since I need the initial conditions
    #   to be stored in the lambda. I can get around this by passing them as
    #   default arguments, which are stored when the function is created.
    # It's a little janky, but it works.
    # TODO: There's a better way to do this without having these two mirrored
    #   lines, but whatever.

    # Define position functions
    traj1 = lambda _t, _x=x[0][i-1], _v=v[0][i-1], _start=start: \
        _x + _v*(_t-_start) + 0.5*g*(_t-_start)**2

    traj2 = lambda _t, _x=x[1][i-1], _v=v[1][i-1], _start=start: \
        _x + _v*(_t-_start) + 0.5*g*(_t-_start)**2

    # Define velocity functions
    vel1 = lambda _t, _v=v[0][i-1], _start=start: \
        _v + 0.5*g*(_t-_start)
    vel2 = lambda _t, _v=v[1][i-1], _start=start: \
        _v + 0.5*g*(_t-_start)

    traj1s.append([start, end, traj1, vel1])
    traj2s.append([start, end, traj2, vel2])


# Sample the trajectories at each of the points
print("Sampling points")
x = [[],[]]
v = [[],[]]
sample_times = np.arange(0,tmax,dt)
for time in sample_times:
    while time > traj1s[0][1]:
        traj1s.pop(0)
        traj2s.pop(0)

    x[0].append(traj1s[0][2](time))
    x[1].append(traj2s[0][2](time))

    v[0].append(traj1s[0][3](time))
    v[1].append(traj2s[0][3](time))



print("Plotting trajectories")
################### Old piecewise method ######################
# This evaluates

#print("Generating boolean list")
# Build list of boolean conditions defining where each
#   function applies.
# TODO: This is far and away the slowest part of the program. Find a more efficient algorithm.
#   -> Maybe: Rather than doing all this list indexing
#conds = np.zeros([len(traj1s), len(sample_times)])
#for i in range(len(traj1s)):
#    start = traj1s[i][0]
#    end = traj1s[i][1]
#    for j in range(len(sample_times)):
#        if sample_times[j] >= start and sample_times[j] < end:
#            conds[i][j] = True

#print("Calculating piecewise")
#data_points1 = np.piecewise(sample_times, conds, map(list, zip(*traj1s))[2])
#data_points2 = np.piecewise(sample_times, conds, map(list, zip(*traj2s))[2])
###############################################################
plt.subplot(312)
plt.plot(sample_times, x[0], sample_times, x[1])
plt.xlabel('t')
plt.ylabel('x')


################## Autocorrelation #########################

print("Plotting correlation")
plt.subplot(313)
# plt.acorr(x[0][::10], usevlines=False, normed=True, maxlags=100, linestyle='-')
correlated = np.correlate(x[0][::100], x[0][::100], mode=2)
correlated = correlated[(correlated.size-1)//2 :]
plt.plot(correlated / max(correlated), 'r-')
plt.show()
