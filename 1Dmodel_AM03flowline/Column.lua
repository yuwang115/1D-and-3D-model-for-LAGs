-- ## Set two surface level
UppSurf = 3000
LowSurf = 0
--UppSurf = 2500  --for restart
--LowSurf = 0


-- ## Set two surface temperature
--UppTemp = 243.15 -- -30
--LowTemp = 271.15 -- -2
UppTemp = 244.15 -- -29
LowTemp = 270.45 -- -2.7


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
  
