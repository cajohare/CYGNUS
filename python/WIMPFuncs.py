import numpy as np
from numpy import pi, sqrt, exp, zeros, size, shape
from numpy.linalg import norm
from scipy.special import erf
import LabFuncs
import Params

#================================WIMPFuncs.f95=================================#
# Just contains a bunch of boring but useful/essential stuff

# Contents:
# 1. Halo models
#    SHM: Loads in old SHM parameters
#    SHMpp: Loads in updated SHM paramers

# 2. Loading recoil distribution for WIMPs (RD_wimp)
#    WIMPRecoilDistribution: Loads RD_wimp to be used by likelihood
#    WIMPRD_Energy: Energy dependent recoil distribution
#    WIMPRD_3D: Direction dependent recoil distribution


# 3. Energy dependnt and direction dependent rates
#    WIMPRate_Direction: Directional recoil rate d^2R(t)/dEdO
#    WIMPRate_Energy: Non-directional recoil rate dR(t)/dE
#	 MaxWIMPEnergy: Maximum recoil energy for a given date
#==============================================================================#

m_p = 0.9315*1e6
c_cm = 3.0e8*100.0 # speed of light in cm/s
GeV_2_kg = 1.0e6*1.783e-33 # convert GeV to kg
Jan1 = 2457755.0 # Julian Day of January 1 2018


#==============================================================================#
#-------------------- Energy-Time dependent recoil rate------------------------#
def MinimumWIMPSpeed(E_r,A,m_chi):
    m_N = m_p*A # mass of nucleus
    mu_p = 1.0e6*m_chi*m_p/(1.0e6*m_chi + m_p) # reduced proton mass
    m_N_keV = A*0.9315*1.0e6 # nucleus mass in keV
    mu_N_keV = 1.0e6*m_chi*m_N_keV/(1.0e6*m_chi + m_N_keV) # reduced nucleus mass
    v_min = (sqrt(2.0*m_N_keV*E_r)/(2.0*mu_N_keV))*3.0e8/1000.0 # vmin in km/s
    return v_min

def MaxWIMPEnergy(A,v_e,m_chi):
    m_N = m_p*A
    mu_N = 1.0e6*m_N*m_chi/(1.0e6*m_chi+m_N)
    E_max_lim = 2.0*mu_N*mu_N*2.0*((v_esc+sqrt(sum(v_lab**2.0)))*1000.0/3.0e8)**2.0/m_N
    return E_max_lim


def WIMPRate_Energy(E_r,t,WIMP,HaloModel,Nuc,Loc):
    # relevant constants
    A = Nuc.MassNumber # mass number of nucleus
    m_chi = WIMP.Mass
    mu_p = 1.0e6*m_chi*m_p/(1.0e6*m_chi + m_p)
    sigma_p = WIMP.SICrossSection
    v_0 = sqrt(2.0)*HaloModel.Dispersion
    v_esc = HaloModel.EscapeSpeed
    rho_0 = HaloModel.Density
    N_esc = HaloModel.Normalisation
    FF = LabFuncs.FormFactorHelm(E_r,A)**2.0
    v_min = MinimumWIMPSpeed(E_r,A,m_chi)
    R0 = (c_cm*c_cm)*((rho_0*1.0e6*A*A*sigma_p)/(2*m_chi*GeV_2_kg*mu_p*mu_p))

    # init
    ne = size(E_r)
    nt = size(t)
    dR = zeros(shape=(size(t),size(E_r)))
    for i in range(0,nt):
        v_e = norm(LabFuncs.LabVelocity(t[i]+Jan1, Loc, HaloModel))

        # Mean inverse speed
        x = v_min/v_0
        y = v_e/v_0
        z = v_esc/v_0
        gvmin = zeros(ne)
        gvmin[(x<abs(y-z))&(z<y)] = (1.0/(v_0*y))
        gvmin[(x<abs(y-z))&(z>y)] = (1.0/(2.0*N_esc*v_0*y))\
                                *(erf(x[(x<abs(y-z))&(z>y)]+y)\
                                -erf(x[(x<abs(y-z))&(z>y)]-y)\
                                -(4.0/sqrt(pi))*y*exp(-z**2))
        gvmin[(abs(y-z)<x)&(x<(y+z))] = (1.0/(2.0*N_esc*v_0*y))\
                                *(erf(z)-erf(x[(abs(y-z)<x)&(x<(y+z))]-y)\
                                -(2/sqrt(pi))*(y+z-x[(abs(y-z)<x)&(x<(y+z))])\
                                *exp(-z**2))
        gvmin[(y+z)<x] = 0.0
        gvmin = gvmin/(1000.0*100.0) # convert to cm^-1 s

        # Compute rate
        dRdE = gvmin*FF # correct for form factor
        dR[i,:] =  dRdE

    dR = dR*3600*24*365*1000.0 # convert to per ton-year
    return dR


#-------------------- Direction dependent recoil rate--------------------------------------#
def WIMPRate_Direction(E,t,WIMP,HaloModel,Nuc,Loc):
    E_r = sqrt(E[:,0]**2 + E[:,1]**2 + E[:,2]**2) # Recoil energy
    x = zeros(shape=shape(E))
    x[:,0] = E[:,0]/E_r # Recoil direction
    x[:,1] = E[:,1]/E_r
    x[:,2] = E[:,2]/E_r

    # relevant constants
    A = Nuc.MassNumber # mass number of nucleus
    m_chi = WIMP.Mass
    mu_p = 1.0e6*m_chi*m_p/(1.0e6*m_chi + m_p)
    sigma_p = WIMP.SICrossSection
    sig_v = HaloModel.Dispersion
    v_esc = HaloModel.EscapeSpeed
    rho_0 = HaloModel.Density
    N_esc = HaloModel.Normalisation
    FF = LabFuncs.FormFactorHelm(E_r,A)**2.0
    v_min = MinimumWIMPSpeed(E_r,A,m_chi)
    R0 = (c_cm*c_cm)*((rho_0*1.0e6*A*A*sigma_p)/(4*pi*m_chi*GeV_2_kg*mu_p*mu_p))

    # init
    ne = size(E_r)
    nt = size(t)
    dR = zeros(shape=(size(t),size(E_r)))
    for i in range(0,nt):
        v_lab = LabFuncs.LabVelocity(t[i], Loc, HaloModel)
        vlabdotq = (x[:,0]*v_lab[0]+x[:,1]*v_lab[1]+x[:,2]*v_lab[2]) # recoil projection

        # Radon transform
        fhat = zeros(shape=shape(E_r))
        fhat[((v_min+vlabdotq)<(v_esc))] = (1/(N_esc*sqrt(2*pi*sig_v**2.0)))\
                                            *(exp(-(v_min[((v_min+vlabdotq)<(v_esc))]\
                                            +vlabdotq[((v_min+vlabdotq)<(v_esc))])\
                                            **2.0/(2*sig_v**2.0))\
                                            -exp(-v_esc**2.0/(2*sig_v**2.0)))
        fhat = fhat/(1000.0*100.0) # convert to cm^-1 s

        # Compute rate
        dRdEdO = fhat*FF # correct for form factor
        dR[i,:] =  dRdEdO

    dR = dR*3600*24*365*1000.0 # convert to per ton-year
    return dR
