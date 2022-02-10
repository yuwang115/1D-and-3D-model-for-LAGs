
-- ##########################################################
-- ## Material parameters

yearinsec = 365.25 * 24.0 * 60.0 * 60.0
Pa2MPa = 1.0E-06
rhoi = 917.0 / (1.0e6 * yearinsec^2)
gravity = -9.81 * yearinsec^2

-- ## Ice fusion latent heat (0.335×106 J kg−1)
Lf = 0.335e06 * yearinsec^2.0

-- ## Sea level elevation
zsl = 0.0

-- ## ocean water density
rhoo = 1029.0 / (1.0e6 * yearinsec^2.0)

-- ## Specific heat of ocean water (4218 J kg−1 K−1)
cw = 4218.0 * yearinsec^2.0

--  For Glen's flow law (Paterson 2010)
n = 3.0
m = 1.0/n
A1 = 2.89165e-13*yearinsec*1.0e18 
A2 = 2.42736e-02*yearinsec*1.0e18 
Q1 = 60.0e3
Q2 = 115.0e3
Tlim = -10.0

-- ## GLToleranceInit=1.0e-1
GLTolerance=1.0e-3


--  Temperature of the simulation in Celsius
--  with the formula for A works only if T > -10
Tc=-1.0



-- ##########################################################
-- hard coded paths to forcing data

-- VELOCITY_DATA = "/projappl/project_2002875/data/antarctic/antarctica_m2slim.nc"
-- BETA_GUESS    = "/projappl/project_2002875/data/antarctic/aa_v3_e8_l11_beta.nc"

datadir = "/projappl/project_2002875/data/antarctic/"
outdir  = "./vtuoutputs"



-- ##########################################################
-- ## functions for material parameters and conditions



-- ## Set geometry to floatation
function UpperSurface(ShelfThick)
  UpperSurface = ShelfThick * (1.0 - rhoi/rhoo)
  return UpperSurface
end
function LowerSurface(ShelfThick)
  LowerSurface = 0.0 - ShelfThick * (rhoi/rhoo)
  return LowerSurface
end

-- ## identify the grounding line in a 2D domain based on 
-- ## geometry (this is a crude approximation).
function groundingline(thick,bed,surf)
  if ((surf - thick) > (bed + 100.0) or (surf - thick) <= (bed + 1.0)) then
     gl_mask = -1.0
  else
     gl_mask = 1.0
  end
  return gl_mask
end

-- ## for imposing a velocity condition based on temperature
-- ## input is temperature relative to pressure melting
function tempCondition(temp,tempCutoff)
  if (temp < tempCutoff) then
    cond = 1.0
  else
    cond = -1.0
  end
  return cond
end

-- ## Impose basal mass balance as lower surface normal velocity 
-- ## under shelf for steady simulations (e.g. inversions).
-- ## Note on sign:
-- ## Normal velocity is taken to be positive "outward" i.e. 
-- ## approx. downward under the ice shelf.
-- ## bmb is taken to be positive for mass gain and negative for 
-- ## mass loss.
-- ## So negative bmb => positive normal velocity
function bmb_as_vel(bmb,gmask)
  if (gmask < -0.5) then
    vel = -1.0 * bmb
  else
    vel = 0.0
  end
  return vel
end

-- ## assume ice goes from zero to thick in the vertical (coord 3).
-- ## topVel is the upward velocity at the ice upper surface.
-- ## should be called on unit height mesh (before extrusion).
function velocityProfile(Z, topVel)
  vel = topVel * Z
  return vel
end

-- ## 
function bmb_over_time(time)
  bmb = time / 200.0
  return bmb
end

-- ## more accurate identification of grounding line for a body
-- ## force condition if GroundedSolver is present.
-- ## Also checks velocity: allow coarse mesh for low speed GL.
function glCondition(glMask,vel,velCutoff)
  if ( (glMask < 0.5) and (glMask > -0.5) ) then 
    cond = 1.0
  else
    cond = -1.0
  end
  if ( vel < velCutoff ) then
    cond = -1.0
  end
  return cond
end

-- ## set the distance at GL to non-zero values depending on flow speed...
-- ## (experimental; not currently used)
function distBF(vel)
  if vel > refvel then
    dist = 0.0
  else
    dist = 100000.0 * (refvel - vel)/refvel
  end
  return dist
end

-- ## function for setting an upper limit to mesh size based on distance
-- ## (e.g. distance from grounding line)
function refinebydist(distance)
  factor = distance/GLdistlim
  if (factor < 0.0) then
    factor = 0.0
  end	   
  if (factor > 1.0) then
    factor = 1.0
  end	   
  Mmax = Mmaxclose*(1.0-factor) + Mmaxfar*factor
  return Mmax
end

-- ## function for setting an upper limit for mesh size
function setmaxmesh(gldist,bdist,vel,glmask)
  gldistfactor = gldist/GLdistlim
  if (gldistfactor < 0.0) then
    gldistfactor = 0.0
  end	   
  if (gldistfactor > 1.0) then
    gldistfactor = 1.0
  end	   
  bdistfactor = bdist/Bdistlim
  if (bdistfactor < 0.0) then
    bdistfactor = 0.0
  end	   
  if (bdistfactor > 1.0) then
    bdistfactor = 1.0
  end	   
  if (gldistfactor < bdistfactor) then
    distfactor = gldistfactor
  else
    distfactor = bdistfactor
  end
  velfactor = vel/refvel
  if (velfactor > 1.0) then
    velfactor = 1.0
  end
  if velfactor < 0.5 then
    velfactor = 0.5
  end
  Mmax = ( Mmaxclose*(1.0-distfactor) + Mmaxfar*distfactor ) / (velfactor)
  if (glmask < 0.5) then
      if (Mmax > Mmaxshelf) then
      Mmax = Mmaxshelf
    end
  end
  return Mmax
end

-- ## function for setting a lower limit for mesh size
function setminmesh(gldist,bdist,glmask)
  effectivedist = gldist - GLdistlim
  if (effectivedist < 0.0) then
    effectivedist = 0.0
  end
  distfactor = effectivedist / distscale
  if (distfactor > 1.0) then
    distfactor = 1.0
  end
  if (bdist < Bdistlim) then
    distfactor = 0
  end
  Mmin = Mminfine*(1.0-distfactor) + Mmincoarse*distfactor
  if (glmask < 0.5) then
    if (Mmin > (Mmaxshelf - 50.0) ) then
      Mmin = Mmaxshelf - 50.0
    end
  end
  return Mmin
end

-- ## set the lower surface for a given upper surface and thickness
function getlowersurface(upp_surf,thick,bed)
  if (thick < MINH) then
    thick = MINH
  end
  if ((upp_surf - thick) > bed) then
    low_surf = upp_surf - thick
  else
    low_surf = bed
  end
  return low_surf
end  		

-- ## set the upper and lower surfaces to floatation
function floatUpper(thick,bed)
  if (thick < MINH) then
    thick = MINH
  end
  if ( (thick*rhoi/rhoo) >= -(bed) ) then
    upp_surf = bed + thick
  else
    upp_surf = thick - thick*(rhoi/rhoo)
  end
  return upp_surf
end

function floatLower(thick,bed)
  if (thick < MINH) then
    thick = MINH
  end
  if ( (thick*rhoi/rhoo) >= -bed ) then
    low_surf = bed
  else
    low_surf = -thick*rhoi/rhoo
  end
  return low_surf
end  		


-- ## variable timestepping (TODO: dt_init and dt_max and dt_incr should be passed in)
function timeStepCalc(nt)
  dt_init = 0.00001
  dt_max = 0.25
  dt_incr = 1.15
  dt = dt_init * 1.2^nt 
  if ( dt > dt_max ) then
    dt = dt_max
  end
  return dt
end



function sw_pressure(z)
  if (z >  0) then
    p=0.0
  else
    p=-rhoo*gravity*z
  end
  return p
end

function init_velo1(v, g1, g2, zs, zb, z)
  gt = math.sqrt(g1*g1 + g2*g2)
  vin1=-v*(g1/gt)*z
  return vin1
end

function init_velo2(v, g1, g2, z)
  gt = math.sqrt(g1*g1 + g2*g2)
  vin2=-v*(g2/gt)*z
  return vin2
end

-- ## inputs: T is temperature in Kelvin.
-- ## returns temperature in Centigrade.
function K2C(T)
  Th= T - 273.15
  return Th
end  

-- ## inputs: TinC is temperature in Centigrade, p is pressure. 
-- ## Returns temperature relative to pressure melting point.
function relativetemp(TinC,p)
  pe = p
  if (pe < 0.0) then
    pe = 0.0
  end
  Trel = TinC - 9.8e-08*1.0e06*pe
  if (Trel > 0.0) then
    Trel = 0.0
  end
  return Trel
end  

function pressuremelting(p)
  pe = p
  if (pe < 0.0) then
    pe = 0.0
  end  
  Tpmp = 273.15 - 9.8e-08*1.0e06*pe
  return Tpmp
end

-- ## pressure melting temperature approximated as a function of
-- ## depth below ice upper surface.
function pressuremelting_d(depth)
  pe = 0.0 - rhoi * gravity * depth
  if (pe < 0.0) then
    pe = 0.0
  end  
  Tpmp = 273.15 - 9.8e-08*1.0e06*pe
  return Tpmp
end

function initMu(TempRel)
  if (TempRel < Tlim) then
    AF = 3.985e-13 * math.exp( -60.0e03/(8.314 * (273.15 + TempRel)))
  elseif (TempRel > 0) then
    AF = 1.916e03 * math.exp( -139.0e03/(8.314 * (273.15)))
  else
    AF = 1.916e03 * math.exp( -139.0e03/(8.314 * (273.15 + TempRel)))
  end
  glen = (2.0 * AF)^(-1.0/3.0)
  viscosity = glen * yearinsec^(-1.0/3.0) * Pa2MPa
  mu = math.sqrt(viscosity)
  return mu
end

function limitslc(slc)
  slco = slc
  if (slco < 1.0e-06) then
    slco = 1.0e-06
  end
  return slco
end

-- Condition to impose no flux on the lateral side applied if
-- surface slope in the normal direction positive (should result in inflow)
-- and greater than 50m/km
function OutFlow(N1,N2,G1,G2) 
  cond=N1*G1 + N2*G2 - 0.05
  return cond
end

function evalcost(velx,vx,vely,vy)
  if (math.abs(vx)<1.0e06) then
     Cost1=0.5*(velx-tx(1))*(velx-vxy)
  else
     Cost1=0.0
  end
  if (math.abs(vy)<1.0e06) then
     Cost2=0.5*(vely-vy)*(vely-vxy)
  else
     Cost2=0.0
  end   
  return Cost1+Cost2
end

function evalcostder(vel,v)
  if (abs(v) < 1.0e06) then
    return (vel - v)
  else
    return 0.0
  end
end  

function initbeta(slc)
  dummy = slc + 0.00001
  return  math.log(dummy)
end  



