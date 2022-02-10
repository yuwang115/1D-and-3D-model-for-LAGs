-- ## input the upper and lower surface normal volocity to get vertical profile. Simulate smb and bmb. need to set Upper Surf level

function SmbBmbProfile(Z, uppVel, lowVel)
  vel = uppVel * Z + lowVel*(1-Z)
  return vel
end


-- This function records the uniform change of surface mass balance from Jetty Peninsula to AM04.

-- ## thermal properties
function conductivity(T)
 k=9.828*math.exp(-5.7E-03*T)
 return k
end

-- ## test conductivity for permeable marine ice layer
function conductivity_above195(T)
 k=6.727*math.exp(-4.1E-03*T)
 return k
end

-- ## 
function bmb_upstream_AM04(time)
 if (time > 775) then
    bmb = -1.3
 else 
    bmb = 0.0
 end
 return bmb
end

-- ## 
function smb_upstream_AM04(time)
 if (time > 775) then
    smb = ((0.24/224)*time)-0.7214
 else 
    smb = 0.1
 end
 return smb
end

-- ## 
function smb_AM04_AM01(time)
    smb = ((0.19/100)*time)+0.35
 return smb
end

-- ## 
function basal_temp(time) --764 years
    temp = ((1.3/764)*time)+269.45
 return temp
end

-- ## 
function bmb_AM03(time)
if (time > 999) then
    bmb = ((-9)/764)*time+21.78
 else  
    bmb = 0
 end
 return bmb
end

-- ## 
function T_GZ_AM03(timestep)--759years 
   T = (1.4/760)*timestep+269.35
 return T
end


-- ## 
function T_MG_GZ(timestep)--1339
   T = (-1.1/1340)*timestep+270.45
 return T
end


function Strain_MG_GZ(Z, time) --1000years
  uppVel = (-1.197e-08*time^2+0.0001213*time-0.3167)/2
  lowVel = -(-1.197e-08*time^2+0.0001213*time-0.3167)/2
  vel = uppVel * Z + lowVel*(1-Z)
  return vel
end

function Strain_GZ_AM03(Z, time)  --759years
  uppVel = (-1.191e-08*time^2+7.523e-05*time-0.18)/2
  lowVel = -(-1.191e-08*time^2+7.523e-05*time-0.18)/2
  vel = uppVel * Z + lowVel*(1-Z)
  return vel
end
