import numpy as np
from numpy import erf, cos, sin, pi, sqrt
import LabFuncs
#================================WIMPFuncs.f95==========================================!
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
#=====================================================================================#
#=====================================================================================#


#=====================================================================================!
#-------------------- Energy-Time dependent recoil rate-------------------------------!
def WIMPRate_Energy(E_r,nucleussig_v,rho_0,JD):
    # relevant constants
    A = sum(nucleus) # mass number of nucleus
    m_p = 0.9315*1e6 # mass of proton
    m_N = m_p*A # mass of nucleus
    m_chi_kg = m_chi*1e6*1.783e-33 # mass of dm in kg
    mu_p = 1.0e6*m_chi*m_p/(1.0d6*m_chi + m_p) # reduced proton mass
    c_cm = 3.0e8*100.0 # speed of light in cm/s
    m_N_keV = A*0.9315*1.0e6 # nucleus mass in keV
    mu_N = 1.0d6*m_chi*m_N_keV/(1.0e6*m_chi + m_N_keV) # reduced nucleus mass
    v_min = (sqrt(2.0*m_N_keV*E_r)/(2.0*mu_N))*3.0e8/1000.0 # vmin in km/s
      
       
    # Compute g(vmin) mean inverse speed halo integral
    N_esc = erf(v_esc/(sqrt(2.0d0)*sig_v))-sqrt(2.0/pi)*(v_esc/sig_v)*exp(-v_esc**2.0d0/(2.0d0*sig_v**2.0d0))
    v_0 = sqrt(2.0)*sig_v
    x = v_min/v_0
    v_e = sqrt(sum(v_lab*v_lab))
    y = v_e/v_0
    z = v_esc/v_0
    if (x<abs(y-z)):
        if (z>y):
            gvmin = (1.0d0/(2*N_esc*v_0*y))*(erf(x+y)-erf(x-y)\
                -(4.0d0/sqrt(pi))*y*exp(-z**2.0d0))
        else:
            gvmin = (1.0d0/(v_0*y))
            
    elif: ((abs(y-z).lt.x).and.(x.lt.(y+z))) then
        gvmin = (1.0d0/(2.0d0*N_esc*v_0*y))*(erf(z)-erf(x-y)\
            -(2.0d0/sqrt(pi))*(y+z-x)*exp(-z**2.0d0))
    else:
        gvmin = 0.0
    
    gvmin = gvmin/(1000.0*100.0) # convert to cm^-1 s

    # compute rate
    dRdE = (c_cm*c_cm)*((rho_0*1.0e6*A*A*sigma_p)/(2*m_chi_kg*mu_p*mu_p))*gvmin
    dRdE = dRdE*3600*24*365*1000.0 # convert to per ton-year
    dRdE = dRdE*FormFactorHelm(E_r,A)**2.0 # apply form factor
    return dRdE


    !-------------------- Direction dependent recoil rate--------------------------------------!
def WIMPRate_Direction(E,tbin) result(dRdEdO)
    # Calculate relavant directions
  v_lab = v_lab_all(tbin,:)	! lab velocity   
    E_r = sqrt(sum(E**2.0)) ! Recoil energy
    x = (/E(1)/E_r,E(2)/E_r,E(3)/E_r/) ! Recoil direction
    vlabdotq = sum(x*v_lab) ! Lab-recoil projection


    ! Relevant constants
    A = sum(nucleus) 
    m_p = 0.9315*1e6
    m_chi_kg = m_chi*1.0d6*1.783d-33
    mu_p = 1.0d6*m_chi*m_p/(1.0d6*m_chi + m_p)
    c_cm = 3.0d8*100.0
    m_N = m_p*A
    m_N_keV = A*0.9315*1.0d6
    mu_N = 1.0d6*m_chi*m_N_keV/(1.0d6*m_chi + m_N_keV)
    v_min = (sqrt(2.0d0*m_N_keV*E_r)/(2.0d0*mu_N))*3.0d8/1000.0d0

    ! Compute Radon transform fhat
    N_esc = erf(v_esc/(sqrt(2.0)*sig_v))&
    -sqrt(2.0/pi)*(v_esc/sig_v)*exp(-v_esc**2.0/(2*sig_v**2.0))
    if ((v_min+vlabdotq).le.(v_esc)) then
    fhat = (1/(N_esc*sqrt(2*pi*sig_v**2.0)))*&
    (exp(-(v_min+vlabdotq)**2.0/(2*sig_v**2.0))-&
     exp(-v_esc**2.0/(2*sig_v**2.0)))
    else 
    fhat = 0.0d0
    end if
    fhat = fhat/(1000.0d0*100.0d0) ! convert to cm^-1 s

    # Compute rate
    dRdEdO = (c_cm**2.0)*((rho_0*1.0d6*A**2.0*sigma_p)/(4*pi*m_chi_kg*mu_p**2.0))*fhat
    dRdEdO = dRdEdO*3600*24*365*1000.0d0 # convert to per ton-year
    dRdEdO = dRdEdO*FormFactorHelm(E_r,A)**2.0d0 # correct for form factor

    return dRdEdO
		
			
!=================================Max Energy==================================!
function MaxWIMPEnergy(A,v_lab) result(E_max_lim)
	double precision :: E_max_lim,mu_N,m_N,v_lab(3)
	integer :: A
	m_N = 0.9315*1.0d6*A
	mu_N = 1.0d6*m_N*m_chi/(1.0d6*m_chi+m_N)
	E_max_lim = 2.0d0*mu_N**2.0d0*((v_esc+sqrt(sum(v_lab**2.0)))*1000.0/3.0d8)**2.0d0/m_N
end function MaxWIMPEnergy
	
end module WIMPFuncs
