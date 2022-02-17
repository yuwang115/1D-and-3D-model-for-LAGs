
-- ##########################################################
-- ## Material parameters

yearinsec = 365.25 * 24.0 * 60.0 * 60.0
Pa2MPa = 1.0E-06
MPainPa = 1.0e6

rhoi_si = 917.0  -- ## ice density in kg/m^3
rhoi    = rhoi_si / (1.0e6 * yearinsec^2)

gravity = -9.81 * yearinsec^2
gravity_abs = -gravity

-- ## ocean water density
rhoo = 1027.0 / (1.0e6 * yearinsec^2.0)

-- ## fresh water density
rhow = 1000.0 / (1.0e6 * yearinsec^2.0)

-- ## Ice fusion latent heat)
Lf = 334000.0 -- ## Joules per kg
-- ## Lf = 334.0 -- ## Joules per gram
-- ## Lf = 0.335e06 * yearinsec^2.0 -- ## some old Elmer units conversion thing?  Not used I think...

-- ## Sea level elevation
zsl = 0.0

-- ## Specific heat of ocean water (4218 J kg−1 K−1)
cw = 4218.0 * yearinsec^2.0

-- ## prescribed salinity at ice base for calculating ocean pressure melting temperature.  PSU.
Salinity_b = 35.0 

--  For Glen's flow law (Paterson 2010)
n  = 3.0
ng = n
m  = 1.0/n
A1_SI = 3.985e-13
A2_SI = 1.916e03
A1 = A1_SI*yearinsec*1.0e18
A2 = A2_SI*yearinsec*1.0e18
Q1 = 60.0e3
Q2 = 139.0e3
Tlim = -10.0

-- ## GLToleranceInit=1.0e-1
GLTolerance=1.0e-3

--  Temperature of the simulation in Celsius
--  with the formula for A works only if T > -10
Tc=-1.0


-- ##########################################################
-- ## hard coded paths to forcing data

VELOCITY_DATA = "/projappl/project_2002875/data/antarctic/antarctica_m2slim.nc"
BETA_GUESS    = "/projappl/project_2002875/data/antarctic/aa_v3_e8_l11_beta.nc"
SMB_DATA      = "/projappl/project_2002875/data/antarctic/smbref_1995_2014_mar.nc"

datadir = "/projappl/project_2002875/data/antarctic/"
outdir  = "./vtuoutputs"


-- ##########################################################
-- ## functions for material parameters and conditions

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


-- ## convert a viscosity enhancement factor to a 
-- ## flow enhancement factor, with limits.
function ConvertEF(viscEF,lowerLimit)
  if (viscEF < lowerLimit) then
    viscEF = lowerLimit
  end
  flowEF = 1.0 / viscEF
  if (flowEF < lowerLimit) then
    flowEF = lowerLimit
  end
  return flowEF
end

-- ## Scale the drag coefficient to tune thinning rates...
function ScaleDragCoef(coef)
  scaledCoef = coef*0.1
  return scaledCoef
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

-- ## function to scale a normal velocity slip coefficient 
-- ## at the ice upper surface to restrict emergence 
-- ## velocity (stronger constraint in slow flowing regions).
function ControlEmergVel(vx,vy,upplim)
  ScalingSpeed = 10.0
  speed = math.sqrt(vx*vx + vy*vy)
  SlipCoef = upplim*(1.0 - math.tanh(speed/ScalingSpeed))
end

-- ## function for setting an upper limit to mesh size based on distance
-- ## (e.g. distance from grounding line)
function refinebydist(distance)
  factor = distance/distlim
  if (factor < 0.0) then
    factor = 0.0
  end	   
  if (factor > 1.0) then
    factor = 1.0
  end	   
  Mmax = Mmaxclose*(1.0-factor) + Mmaxfar*factor
  return Mmax
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

-- ## variable timestepping
function timeStepCalc2(nt, dt_init, dt_max, dt_incr)
  dt = dt_init * dt_incr^nt 
  if ( dt > dt_max ) then
    dt = dt_max
  end
  return dt
end

-- ## thermal properties
function conductivity(T)
 k=9.828*math.exp(-5.7E-03*T)
 return k
end

-- ## heat capacity
function capacity(T)
  c=146.3+(7.253*T)
  return c
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
  Trel = TinC + 9.8e-08*1.0e06*pe
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

function pressuremelting_salinity(p)
  pe = p
  if (pe < 0.0) then
    pe = 0.0
  end  

  Tpmp = 273.15 - 5.73e-02*Salinity_b + 9.39e-02 - 7.53e-08*1.0e06*pe

  return Tpmp
end

function initMu(TempRel)
  if (TempRel < Tlim) then
    AF = A1_SI * math.exp( -Q1/(8.314 * (273.15 + TempRel)))
  elseif (TempRel > 0) then
    AF = A2_SI * math.exp( -Q2/(8.314 * (273.15)))
  else
    AF = A2_SI * math.exp( -Q2/(8.314 * (273.15 + TempRel)))
  end
  glen = (2.0 * AF)^(-1.0/n)
  viscosity = glen * yearinsec^(-1.0/n) * Pa2MPa
  return viscosity
end
--  mu = math.sqrt(viscosity)
--  return mu

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
