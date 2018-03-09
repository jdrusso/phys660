import numpy as np

END = -1
COLLISION = 0
BOUNCE = 1

x1, x2 = 1.0, 1.5
v1, v2 = -1.0, -1.0

m1, m2 = 1., 2.

tmax = 10
dt = 0.1
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

    # Update v and x using kinematics equations
    x1 += v1*event_time + 0.5 * g * event_time**2
    v1 += g*event_time

    x2 += v2*event_time + 0.5 * g * event_time**2
    v2 += g*event_time

    v1 *= (-1)

    # If it's a collision, also reverse the top velocity
    if collide:
        v2 *= (-1)

    # If it's a bounce, set x1 to absolute value so it doesn't go through the
    #   floor.
    elif not collide:
        x1 = abs(x1)

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
    print("Found event")
    
print(x)
print(v)
print(t)


# for i in range(1,steps):
#     t += dt
#
#
#     # Update velocities and positions
#
#     v[0,i] = v[0,i-1] - (g*dt)
#     v[1,i] = v[1,i-1] - (g*dt)
#
#     x[0,i] = x[0,i-1] + dt*((v[0,i] + v[0,i-1])/2)
#     x[1,i] = x[1,i-1] + dt*((v[1,i] + v[1,i-1])/2)
#
#     # Check for a bounce/collision
#     distance = x[1,i] - x[0,i]
#     if distance < 0 or abs(distance<)
