-- Define initial conditions
x1 = 1.0 :: Double
v1 = 1.0 :: Double

x2 = 1.0 :: Double
v2 = 1.0 :: Double

dt = .001 :: Double

-- Start both balls at their initial positions and velocities
-- First event will be a bounce, then collision, then bounce, etc...

-- Need:
--  * Times/positions of collisions
--      What if collisions don't occur at a timestep?
--  * Positions/velocities of arbitrary intervals
--  * Autocorrelation of ^
--  * Lyapunov exponent
-- So, get a time series of x, t for the balls


get_events bounce v1 v2 x1 x2
    | bounce == False = 
        take 
    | bounce == True  =
        take 


main = do

    let tmax = 10
    let steps = tmax/dt

    -- Start them both in the air, travelling down.
    -- Calculate if first event will be a collision or a bounce,
    --  and then proceed accordingly.

    -- Time for ball 1 to hit floor
    let t1 = x1 / v1
    let t2 = x2 / v2
    -- Ball will bounce first if t1 < t2
    let bounce_first = t1 < t2

    
