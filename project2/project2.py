import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

END = -1
COLLISION = 0
BOUNCE = 1

x1, x2 = 1.0, 1.5
v1, v2 = 0.0, -1.0

m1, m2 = 1., 1.

tmax = 10
dt = 1.E-3
steps = tmax/dt

g = -1

x = np.zeros([2,1])
v = np.zeros([2,1])
t = np.zeros([1])
# Want to add new values with np.append(x, [[x1],[x2]], axis=1)


x[0,0], x[1,0] = x1, x2
v[0,0], v[1,0] = v1, v2

Ei = 0.5*(m1/(m1+m2))*v1**2 + 0.5*(m2/(m1+m2))*v2**2 + m1/(m1+m2)*x1 + m2/(m1+m2)*x2

eps = 2.2204E-16

print("Initial energy is %f" % Ei)



# General algorithm:
#   1. Determine whether next event will be a bounce or a collision.
#   2. Determine the time, velocities and positions of that event.
#   3. If time > t_max, then stop.
#   4. Using the time, positions, and velocities of the balls at that event,
#       calculate the time of the next event.

# Find what the next event will be, and in how much time
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
    print(event_time)

    # Update v and x using kinematics equations
    #   I.e., find where the balls are at the time of the event
    x1 += v1*event_time + 0.5 * g * event_time**2
    v1 += g*event_time

    x2 += v2*event_time + 0.5 * g * event_time**2
    v2 += g*event_time


    # If it's a collision, also reverse the top velocity
    if collide:
        _v2 = (2*m1)/(m1+m2)*v1 - (m1-m2)/(m1+m2)*v2
        _v1 = (m1-m2)/(m1+m2)*v1 - (2*m2)/(m1+m2)*v2

        v1 = _v1
        v2 = _v2

    # If it's a bounce, set x1 to absolute value so it doesn't go through the
    #   floor.
    elif not collide:
        print("Bounce")
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
print("Found %d events" % len(t))

print(x)
print(v)
print(t)
    
# Create functions describing trajectories between events
# Build a list of each trajectory function, and the domain over
#   which it applies.
funcs1 = []
funcs2 = []
for i in range(1,len(t)):
    start = t[i-1]
    end = t[i]
    # Define the lambda with default arguments which are executed when the lambda is created, not when it's evaluated.
    #   Grade A jank right here.
    func1 = lambda _t, _x=x[0][i-1], _v=v[0][i-1], _start=start: _x + _v*(_t-_start) + 0.5*g*(_t-_start)**2
    #func1 = lambda _t, _x=x[0][i-1], _v=v[0][i-1], _start=start:  ("x[0] = %f, v[0] = %f, start = %f, end = %f" % (_x, _v,_start, end))
    func2 = lambda _t, _x=x[1][i-1], _v=v[1][i-1], _start=start: _x + _v*(_t-_start) + 0.5*g*(_t-_start)**2
    funcs1.append([start, end, func1])
    funcs2.append([start, end, func2])

sample_times = np.arange(0,tmax,dt)

# Build list of boolean conditions defining where each
#   function applies.
conds = np.zeros([len(funcs1), len(sample_times)])
for i in range(len(funcs1)):
    for j in range(len(sample_times)):
        if sample_times[j] < funcs1[i][1] and sample_times[j] >= funcs1[i][0]:
            conds[i][j] = True

print(conds[0][0:10])
print(conds[2][0:10])

data_points1 = np.piecewise(sample_times, conds, map(list, zip(*funcs1))[2])
data_points2 = np.piecewise(sample_times, conds, map(list, zip(*funcs2))[2])
#print(data_points)

print(data_points1[0:10])

plt.plot(sample_times, data_points1, sample_times, data_points2)
plt.show()  
