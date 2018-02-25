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
ra_0 = 150.0 :: Double
rb_0 = 0.0 :: Double
tA = 1 :: Double
tB = 1 :: Double


-- This is a pretty sweet iterative way of calculating this
a_pop t_max dt =
  -- Take this many timesteps from the iterator
  take (floor(t_max / dt) + 1) $
  iterate (\a -> a + a*(-dt)) a_0

-- Here, if track both the A and B populations, we can use the same iterative
--  approach.
pops t_max dt =
  take (floor(t_max / dt) + 1) $
  iterate (\n ->
    [n!!0 + n!!0*(-dt) ,
    (n!!1 + (dt) * (n!!0 - (tA/tB)*n!!1))])
  [a_0, b_0]


-- Here, we simulate the two populations where each type decays into the other.
crossover t_max dt =
  take (floor(t_max / dt) + 1) $
  iterate (\n ->
    [n!!0 + (0.5*n!!1 - n!!0)*(dt),
    n!!1 + (n!!0 - 0.5*n!!1)*(dt)])
  [ra_0, rb_0]

main = do

  -- Set some generic attributes for the plots
  let attrib = [XLabel "{/*3 Time (arbitrary units)}", YLabel "{/*3 Number of Atoms}"]

  let t_max = 10.0
  let dt = 1
  let thepops = transpose $ pops t_max dt
  -- let thepops = transpose $ crossover t_max dt



  plotPathsStyle (attrib ++ [EPS "logplot.eps",
    Custom ("logscale y; set key font \",36\"; set tics font \",36\";" ++
    "set ylabel offset -7; set xlabel offset 0,-1; set lmargin 16; set bmargin 5; set tmargin 2;"++
    "set format y \"10^{%L}\"; set yrange [0.001:1000]; set ytics (1000,10,.1,.001);"++
    "set xrange [0:10]; set xtics (0, 2, 4, 6, 8, 10)") []
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
