
meshdb  = "../mesh2D_LAGS2"


-- ## Min threshold value for Ice thickness (Real)
MINH=150.0

-- ## levels in the vertical
MLEV=20

-- ## controlling steady state iterations
IMIN=200
IMAX=201

Tol=0.01

DPtol = 0.001

-- ## for block preconditioner 
blocktol          = 0.0001  -- linear system convergence tolerance for individual blocks 1-3 in block preconditioner
blocktol_pressure = 0.00001 -- linear system convergence tolerance for individual block 4 in block preconditioner
blockOutInterval  = 20      -- output interval for blocks in block preconditioner
OuterOutInterval  = 50      -- output interval for outer loop in block preconditioner
OuterMaxIter      = 100000  -- Maximum linear iterations for outer loop in block preconditioner
blockMaxIter      = 2000    -- Maximum linear iterations for blocks in block preconditioner
outerLinTol       = 0.8e-4  -- Linear convergence tolerance for outer loop in block preconditioner

BMB_DATA      = "/projappl/project_2002875/data/antarctic/bmb_roms_tincubic.nc"
BMB_DATA2     = "/projappl/project_2002875/data/antarctic/bmb_susheel.nc"
BMB_DATAS2    = "/projappl/project_2002875/data/antarctic/bmb_susheel2.nc"

-- ## BMB_DATA2     = "/projappl/project_2002875/data/antarctic/antarctic_iceshelf_melt_rate.nc"
