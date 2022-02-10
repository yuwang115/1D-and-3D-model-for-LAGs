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

function conductivity_WGJP(T)
 if (T<270.65) then
  k=9.828*math.exp(-5.7E-03*T)
 else
  k=2.457*math.exp(-5.7E-03*T)
 end
 return k
end


function conductivity_JPAM05(T)
 if (T<270) then
  k=9.828*math.exp(-5.7E-03*T)
 else
  k=2.457*math.exp(-5.7E-03*T)
 end
 return k
end

-- ## test conductivity for permeable marine ice layer
function conductivity_above195(T)
 k=6.727*math.exp(-4.1E-03*T)
 return k
end

-- ## heat capacity
function capacity(T)
  cpct = 146.3 + ( 7.253*T )
  return cpct
end

function capacity_WGJP(T)
 if (T<270.65) then
  cpct = 146.3 + ( 7.253*T )
 else
  cpct = (146.3 + ( 7.253*T ))*2
 end
  return cpct
end

function capacity_JPAM05(T)
 if (T<270) then
  cpct = 146.3 + ( 7.253*T )
 else
  cpct = (146.3 + ( 7.253*T ))*2
 end
  return cpct
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
function basal_temp(time)
 if (time > 999) then
    temp = ((1.3/764)*time)+267.75
 else  
    temp = 269.45
 end
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

function bmb_origin_AM04(time)
 if (time > 1176) then
    bmb = -1.3
 else 
    bmb = 0.0
 end
 return bmb
end

-- ## 
function smb_JP_AM05(time)
 smb = 0.07+0.08*(time/130)
 return smb
end

function smb_AM05_AM04(time)
 smb = 0.15+0.06*(time/93)
 return smb
end

function smb_AM04_AM01(time)
 smb = 0.21+0.11*(time/103)
 return smb
end

--------------

function Temp_UP_JP(time)
 T = ((2.92-2.2)/863*time)+270.23
 return T
end

function Temp_WG_JP(time)
 T = ((-1.1)/2429*time)+271.85  -- -1.3 to -2.4
 return T
end

function Temp_JP_AM05(time)
 T = (0.1/130*time)+270.75  -- -2.4 to -2.3
 return T
end

function Temp_AM05_AM04(time)
 T = (0.1/93*time)+270.85  -- -2.3 to -2.2
 return T
end

 
 
