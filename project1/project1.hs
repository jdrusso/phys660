import Graphics.Gnuplot.Simple
-- import Graphics.EasyPlot

-- import Graphics.Rendering.Chart
-- import Data.Colour
-- import Data.Colour.Names
-- import Data.Default.Class
-- import Graphics.Rendering.Chart.Backend.Cairo
-- import Control.Lens


import Data.List

-- Set some constants
a_0 = 200.0 :: Double
b_0 = 5.0 :: Double
tA = 0.5 :: Double
tB = 1 :: Double


-- This is a pretty sweet iterative way of calculating this
a_pop t_max dt =
  -- Take this many timesteps from the iterator
  take (floor(t_max / dt) + 1) $
  iterate (\a -> a + a*(-dt)) a_0

-- Here, if track both the A and B populations, we can use the same iterative
--  approach.
pops t_max dt =
  take (floor(t_max / dt)) $
  iterate (\n ->
    [n!!0 + n!!0*(-dt) ,
    (n!!1 + (dt) * ((a_0/b_0)*n!!0 - (tA/tB)*n!!1))])
  [a_0, b_0]

main = do

  -- Set some generic attributes for the plots
  let attrib = [XLabel "{/*1.75 Time (arbitrary units)}", YLabel "{/*1.75 Number of Atoms}"]

  let t_max = 10.0
  let dt = 0.01
  let thepops = transpose $ pops t_max dt


  plotPathsStyle (attrib ++ [EPS "logplot.eps",
    Custom "logscale y; set key font \",18\"; set tics font \",18\";" []
    -- ,Title ("{/*2 Radioactive Atom Populations vs Time (tau_A: " ++ show tA ++ ", tau_B: " ++ show tB ++ ")}")
    ])$

    -- Zip a list of times in increments of dt to timestamp the population data
    [(defaultStyle{lineSpec=CustomStyle [LineTitle "N_A", LineWidth 10]}, zip [0,dt..] (thepops!!0)),
    (defaultStyle{lineSpec=CustomStyle [LineType 3, LineTitle "N_B", LineWidth 10]}, zip [0,dt..] (thepops!!1))]
  --
  plotPathsStyle (attrib ++ [EPS "plot.eps",
    Custom "key font \",18\"; set tics font \",18\"" []
    -- ,Title ("Radioactive Atom Populations vs Time (tau_A: " ++ show tA ++ ", tau_B: " ++ show tB ++ ")")
    ])$

    -- Zip a list of times in increments of dt to timestamp the population data
    [(defaultStyle{lineSpec=CustomStyle [LineTitle "N_A", LineWidth 10]}, zip [0,dt..] (thepops!!0)),
    (defaultStyle{lineSpec=CustomStyle [LineType 3, LineTitle "N_B", LineWidth 10]}, zip [0,dt..] (thepops!!1))]


  -- -- Not very accurate
  -- let t_max = 5
  -- let dt = 0.75
  -- plotPath (attrib ++ [PNG "inaccurate.png",
  --   Title ("Inaccurate (Tmax: " ++ show t_max ++ ", dt: " ++ show dt ++ ")")])$
  --   -- Zip a list of times in increments of dt to timestamp the population data
  --   zip [0,dt..] (a_pop t_max dt)
  --
  -- -- Oscillatory
  -- let t_max = 12
  -- let dt = 1.5
  -- plotPath (attrib ++ [PNG "oscillatory.png",
  --   Title ("Oscillatory (Tmax: " ++ show t_max ++ ", dt: " ++ show dt ++ ")")])$
  --   -- Zip a list of times in increments of dt to timestamp the population data
  --   zip [0,dt..] (a_pop t_max dt)
  --
  -- -- Unstable
  -- let t_max = 100
  -- let dt = 5
  -- plotPath (attrib ++ [PNG "unstable.png",
  --   Title ("Unstable (Tmax: " ++ show t_max ++ ", dt: " ++ show dt ++ ")")])$
  --   -- Zip a list of times in increments of dt to timestamp the population data
  --   zip [0,dt..] (a_pop t_max dt)
