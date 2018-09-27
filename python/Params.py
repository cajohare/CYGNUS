from numpy import array, sqrt, pi, exp
from scipy.special import erf


#==============================================================================#
class setWIMP:
    def __init__(self,m_chi,sigma_p):
        self.Mass = m_chi
        self.SICrossSection = sigma_p

#------------------------------------------------------------------------------#
class Nucleus:
    def __init__(self,xi,N,Z,J,Sp,Sn):
        self.IsotopicFraction = xi
        self.NumberOfNeutrons = N
        self.NumberOfProtons = Z
        self.MassNumber = N+Z
        self.NuclearSpin = J
        self.ExpProtonSpin = Sp
        self.ExpNeutronSpin = Sp
        self.SDEnhancement = (4.0/3.0)*((J+1.0)/J)*(Sp-Sn)**2.0

#              (xi,      N,   Z,    J,     Sp,      Sn)
F19 =   Nucleus(1.0,    10,   9,  0.5,  0.421,   0.045)
Xe129 = Nucleus(0.265,  75,  54,  0.5,  0.046,   0.293)
Xe131 = Nucleus(0.212,  75,  54,  1.5, -0.038,   0.242)

#==============================================================================#


#==============================================================================#
class Halo:
    def __init__(self,rho_0,v_LSR,sig_v,v_esc,v_pec):
        self.Density = rho_0
        self.RotationSpeed = v_LSR
        self.Dispersion = sig_v
        self.EscapeSpeed =  v_esc
        self.PeculiarVelocity = v_pec
        self.Normalisation = erf(v_esc/(sqrt(2.0)*sig_v))-\
                            sqrt(2.0/pi)*(v_esc/sig_v)*\
                            exp(-v_esc**2.0/(2.0*sig_v**2.0))

SHM = Halo(0.3,
        220.0,
        156.0,
        544.0,
        array([11.1,12.2,7.3]))

SHMpp = Halo(0.55,
        233.0,
        164.8,
        580.0,
        array([11.1,12.2,7.3]))
#------------------------------------------------------------------------------#


#==============================================================================#
class Location:
    def __init__(self,lat,lon):
        self.Latitude = lat
        self.Longitude = lon

Boulby = Location(54.5591,0.8310)
GranSasso = Location(42.4691, 13.5654)
Kamioka = Location(36.2381, 137.1863)
SNOlab = Location(46.4719, -81.1868)
Stawell = Location(-37.0576, 142.7754)
Oahu = Location(21.4389, -158.0001)
GuantanamoBay = Location(20.0117, -75.1216)
Pyongyang = Location(39.0392, 125.7625)
#------------------------------------------------------------------------------#
