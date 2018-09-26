#======================================================================================#
#
# This code has been translated from CYGNUS/code/LabFuncs.f95
#
# Currently it contains: 
# Lab velocity:
# LabVelocity: Full lab velocity in (N,W,Z) with Earth rotation
# JulianDay: JulianDay at dd-mm-yyyy hh:hh
#
# Solar direction:
# EarthSunDistance: Distance between Earth and Sun as a function of time
# SolarDirection: Direction of the sun at a given time
#
# Co-ordinate transformations
# eqt2lab: Equatorial system to laboratory system
# gal2eqt: Galactic system to equatorial system
# gal2lab: Galactic system to lab system
#=======================================================================================!


import numpy as np
from numpy import cos, sin, pi, floor



#---------------------------------------------------------------------------------------#
# Relevant velocities
vv_earthrev = 29.8
v_LSR = 220.0
v_pec = np.array([11.1,12.2,7.3])

# Earth orbital params
eccentricity = 0.016722
eccentricity_deg = 0.9574
orb_long_ecliptic = 13.0+1.0
lat_ecl_gal = np.array([-5.5303,59.575,29.812])
long_ecl_gal = np.array([266.141,-13.3485,179.3212])

# Other constants
AstronomicalUnit = 1.49597892e11 # Astronomical Unit
EarthRadius = 6371.01*1000.0 # Earth Radius
Integral_inv_EarthSun_sq = 4.468864372000642e-23 # integral(1/R^2) over 1 year
Jan1 = 2457755.0 # Julian Day of January 1 2018


#===================================Lab Velocity========================================#
#---------------------------------------------------------------------------------------#
# FULL LAB VELOCITY AS A FUNCTION OF DAY IN THE YEAR 
# (i.e day = 0.5 is midday on Jan 1 2018)
def LabVelocity(JD, Location):

    lat = Location[0]
    lon = Location[1]

    # Convert day into phase of Earth rotation t_lab
    #JD = day+Jan1
    UT = 24*(JD+0.5-floor(JD+0.5)) #Universal time
    MJD = JD - 2400000.5 #Modified Julian Day
    T_0 = (floor(MJD)-55197.5)/36525.0
    t_GAST = (101.0308 + 36000.770*T_0 + 15.04107*UT)/15.0
    t_lab = t_GAST + lon/15
    t_lab = 15*t_lab #Lab time in degrees

    
    # Galactic (LSR) Rotation
    vtemp = np.array([0.0,v_LSR,0.0])
    v_galrot = gal2lab(vtemp,t_lab, lat) #transform to lab co-ords
    
    # Peculiar solar Motion
    vtemp1 = v_pec
    v_solar = gal2lab(vtemp1,t_lab, lat) # transform to lab co-ords
    
    #Earth's revolution (first calculate in galactic frame then transform)
    e = eccentricity
    lambda_0 = orb_long_ecliptic
    L = 281.0298 + 36000.77*T_0 + 0.04107*UT
    g = 357.9258 + 35999.05*T_0 + 0.04107*UT
    lambda_sun = L + (1.915 - 0.0048*T_0)*sin(g*pi/180.0) \
         + 0.020*sin(2*g*pi/180.0)
    beta = lat_ecl_gal
    lambda_i = long_ecl_gal
    v_earthrev1 = vv_earthrev*(1-e*sin(pi/180.0*(lambda_sun-lambda_0)))* \
         (cos(beta*pi/180.0)*sin(pi/180.0*(lambda_sun-lambda_i)))
    v_earthrev = gal2lab(v_earthrev1,t_lab, lat) #transform to lab co-ords
    
    # Earth's rotation
    v_earthrot = 0.465102*cos(lat*pi/180)*np.array([0.0,-1.0,0.0]) #already in lab co-ords
    
    # Add them all together (delete as needed)
    v_lab = np.array([0.,0.,0.])
    v_lab += v_earthrot
    v_lab += v_earthrev
    v_lab += v_solar
    v_lab += v_galrot
    
    return v_lab
    
    
def JulianDay(month, day, year, hour):  #Calculates time in JD for a given date
    year_r = year+4800-floor((14-month)/12.0)
    month_r = month+12*floor((14-month)/12.0)-3
    JulianDay = day + floor((153*month_r+2)/5.0) + 365*year_r + floor(year_r/4.0) - floor(year_r/100.0) + floor(year_r/400.0) - 32045 + (hour-12.0)/24.0
    return JulianDay

    
        
#===================================Solar direction=====================================#
#---------------------------------------------------------------------------------------#
def EarthSunDistance(day): # Earth-sun distance at time = day
    JD = day+Jan1
    D = JD-2451545.0
    g = 357.529 + 0.98560028*D
    g = g*pi/180.0
    r_es = 1.00014 - 0.01671*cos(g) - 0.00014*cos(2*g)
    r_es = r_es*AstronomicalUnit
    return r_es
    
#---------------------------------------------------------------------------------------#
def SolarDirection(JD,Location): # Solar direction in lab coords at time = day

    lat = Location[0]
    lon = Location[1]

    # Compute RA and dec of Sun
    #JD = day+Jan1
    n = JD - 2451545.0
    Omega = 2.1429-0.0010394594*n
    L = 4.8950630 + 0.017202791698*n
    g = 6.2400600 + 0.0172019699*n
    ll = L+0.03341607*sin(g) + 0.00034894*sin(2*g) - 0.0001134 - 0.0000203*sin(Omega)
    ep = 0.4090928 - 6.214e-9*n + 0.0000396*cos(Omega)
    ra = np.arctan2((cos(ep)*sin(ll)),cos(ll)) # Right ascension of Sun
    dec = np.arcsin(sin(ep)*sin(ll)) # Declination of sun

    # Solar vector
    x_sun1 = np.array([0.,0.,0.])
    x_sun1[0] = cos(dec)*cos(ra)
    x_sun1[1] = cos(dec)*sin(ra)
    x_sun1[2] = sin(dec)

    # Lab time conversion
    UT = 24*(JD+0.5-floor(JD+0.5))
    MJD = JD - 2400000.5
    T_0 = (floor(MJD)-55197.5)/36525.0
    t_GAST = (101.0308 + 36000.770*T_0 + 15.04107*UT)/15.0
    t_lab = t_GAST + lon/15.0
    t_lab = 15*t_lab # DEGREES

    # Convert vector from equatorial system into lab system
    x_sun = eqt2lab(x_sun1,t_lab,lat)
    return x_sun
#---------------------------------------------------------------------------------------#






#-------------------------------Coordinate trans.---------------------------------------#   
def eqt2lab(vp,t_lab, lat): # Equatorial (x_e,y_e,z_e) to Laboratory (N,W,Z)
    t = t_lab*pi/180.0
    latr = lat*pi/180.0
    v = vp*0.0
    v[0] = -cos(t)*sin(latr)*vp[0] - sin(t)*sin(latr)*vp[1] + cos(latr)*vp[2]
    v[1] = sin(t)*vp[0] - cos(t)*vp[1]
    v[2] = cos(t)*cos(latr)*vp[0] + cos(latr)*sin(t)*vp[1] + sin(latr)*vp[2]
    return v

def gal2eqt(vp): # Galactic (x_g,y_g,z_g) to Equatorial (x_e,y_e,z_e) 
    v = 0.0*vp
    v[0] = -0.06699*vp[0] + 0.4927*vp[1] - 0.8676*vp[2]
    v[1] = -0.8728*vp[0] -0.4503*vp[1] -0.1884*vp[2]
    v[2] = -0.4835*vp[0] + 0.7446*vp[1] + 0.4602*vp[2]
    return v

def gal2lab(v,t_lab, lat): # Galactic (x_g,y_g,z_g) to Laboratory (N,W,Z)
    vp = gal2eqt(v)
    return eqt2lab(vp, t_lab, lat)
#---------------------------------------------------------------------------------------#
#=======================================================================================#