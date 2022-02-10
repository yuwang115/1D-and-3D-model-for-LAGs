-- ## Set two surface level
--UppSurf = 1320   -- steady state & UG_JP
UppSurf = 1500
LowSurf = 0

--UppSurf = 690  -- JP_AM05
--LowSurf = 0

-- ## Set two surface temperature

UppTemp = 248.15 -- -25
--UppTemp = 250.65 -- -22.5
--UppTemp = 251.65 -- -21.5
--LowTemp = 269.45 -- -3.7
--UppTemp = 253.55 -- -19.6
--LowTemp = 270.95 -- -2.2
--LowTemp = 270.23 -- -2.92
LowTemp = 271.85 -- -1.3


datadir = "/projappl/project_2002875/data/antarctic/"
outdir = "./vtuoutputs"
meshdb = "./ColumnMesh"

-- ## Min threshold value for Ice thickness (Real)
MINH=100.0

-- ## levels in the vertical
MLEV=100

-- ## controlling steady state iterations
IMIN=10
IMAX=100

Tol=0.01

-- ## for block preconditioner 
blocktol=0.001

-- ## mesh update
function fs_lower(lower)
  mesh = lower-LowSurf 
  return mesh
 end
 
 function fs_upper(upper)
  mesh = upper-UppSurf 
  return mesh
 end
  
