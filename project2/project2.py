import numpy as np

x1 = 1.0
x2 = 1.5
v1 = -1.0
v2 = -1.0

m1 = 1.
m2 = 2.

tmax = 10
dt = 0.1
steps = tmax/dt

g = -1

x = np.zeros([2,steps])
v = np.zeros([2,steps])

x[0,0], x[1,0] = x1, x2
v[0,0], v[1,0] = v1, v2

Ei = 0.5*(m1/(m1+m2))*v1^2 + 0.5*(m2/(m1+m2))*v2^2 + m1/(m1+m2)*x1 + m2/(m1+m2)*x2

eps = 2.2204E-16

print("Initial energy is %f" % Ei)


# Check if the first event will be a bounce or a collision
# Time for the first ball to hit the ground
bounce_time1 = [x for x in np.roots[0.5*g, v[0,0], x[0,0]] if x > 0][0]
bounce_time2 = [x for x in np.roots[0.5*g, v[0,0], x[1,0]] if x > 0][0]

# Check to see if the first event is a collision
collide = bounce_time2 < bounce_time1

# Find the collision time
if collide:
    event_time = (x[0,0]-x[1,0]) / (v[1,0]-v[0,0])

#
if not collide:
    bounce = True

collided = False
bounced = False
t == 0
for i in range(steps):

    # If a collision happened in the last step, reverse the direction of both balls
    #   Use the new velocity and position to calculate whether it will bounce or
    #   collide next.
    if collided:

        collided = False
        bounce = True
        pass

    # If a bounce happened in the last step, only reverse the direction of the bottom ball.
    #   Also set its position to the absolute value of its position.
    #   The next event WILL be a collision.
    if bounced:

        bounced = False
        collide = True
        pass


    # If there's an impending collision, keep looping until it happens
    if collide:

        # If the collision is going to happen withinin the next timestep, set
        #   collided to True
        continue


    # If there's an impending bounce, keep looping until it happens
    elif bounce:

        # If the bounce will happen within the next timestep, set bounced to
        #   True.

        continue

    else:

        pass






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
