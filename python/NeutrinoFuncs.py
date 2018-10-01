from numpy import pi, sqrt, exp, zeros, size, shape, array
from numpy import trapz, interp, loadtxt
from numpy.linalg import norm
import LabFuncs
import Params
from Params import nufile_root, nufile_dir, nuname, n_Enu_vals
from Params import mono, NuMaxEnergy, NuFlux, NuUnc, whichsolar, n_nu_tot


#================================NeutrinoFuncs=================================#
# Contents:

N_A = 6.02214e23 # Avocado's constant
sinTheta_Wsq = 0.2387e0 # sin(Theta_W) weak mixing angle
G_F_GeV = 1.16637e-5 # GeV**-2 ! Fermi constan in GeV

#========================================Neutrino data=========================#
def GetNuFluxes(E_th,Nuc):
    # Reads each neutrino flux data file
    # the energies are stored in E_nu_all, fluxes in Flux_all

    # Figure out which backgrounds give recoils above E_th
    E_r_max = MaxNuRecoilEnergies(Nuc) # Max recoil energy for neutrino
    n_nu = size(E_r_max>E_th)
    sel = range(1,n_nu_tot+1)*(E_r_max>E_th)
    sel = sel[sel!=0]-1
    E_nu_all = zeros(shape=(n_Enu_vals,n_nu))
    Flux_all = zeros(shape=(n_Enu_vals,n_nu))
    Flux_err = zeros(shape=(n_nu))
    Flux_norm = zeros(shape=(n_nu))
    Solar = zeros(n_nu,dtype=bool)

    ii = 0
    for s in sel:
        if mono[s]:
            E_nu_all[0,ii] = NuMaxEnergy[s]
            Flux_all[0,ii] = NuFlux[s]
        else:
            data = loadtxt(nufile_dir+nuname[s]+nufile_root,delimiter=',')
            E_nu_all[:,ii],Flux_all[:,ii] = data[:,0],data[:,1]
            Flux_all[:,ii] = Flux_all[:,ii]*NuFlux[s]

        Flux_norm[ii] = NuFlux[s]
        Flux_err[ii] = NuUnc[s] # Select rate normalisation uncertainties
        Solar[ii] = whichsolar[s]
        ii = ii+1
    NuBG = Params.Neutrinos(n_nu,Solar,E_nu_all,Flux_all,Flux_norm,Flux_err)
    return NuBG

 #-----------------------------------------------------------------------------#
def MaxNuRecoilEnergies(Nuc): # Max recoil energies
    m_N = 0.93141941*(Nuc.MassNumber)*1.0e6
    E_r_max = 2*m_N*(1000.0*NuMaxEnergy)**2.0/(m_N+1000*NuMaxEnergy)**2.0
    return E_r_max


#===================================nu spectra=================================#
def NuRate(E_r,t,Expt,NuBG,nu_i): # Time-Energy
    E_nu_all = NuBG.Energy
    Flux_all = NuBG.Flux
    n_nu = NuBG.NumberOfNeutrinos
    Solar = NuBG.SolarLabel

    Nuc = Expt.Nucleus
    Loc = Expt.Location

    if Expt.Directional:
        ne = size(E_r)/3
        dR = zeros(shape=(ne,n_nu))
        if Solar[nu_i]:
            dR[:,nu] = dRdEdO_solarnu(E_r,t,E_nu_all[:,nu_i],Flux_all[:,nu_i],Nuc,Loc)
        else:
            dR[:,nu] = dRdEdO_isonu(E_r,E_nu_all[:,nu_i],Flux_all[:,nu_i],Nuc)
    else:
        ne = size(E_r)
        dR = zeros(shape=(ne,n_nu))
        dR[:,nu] = dRdE_nu(E_r,t,Solar[nu_i],E_nu_all[:,nu_i],Flux_all[:,nu],Nuc)
    return dR


def dRdE_nu(E_r,t,sol,E_nu,Flux,Nuc):
    N = Nuc.NumberOfNeutrons
    Z = Nuc.NumberOfProtons
    Q_W = 1.0*N-(1-4.0*sinTheta_Wsq)*Z # weak nuclear hypercharge
    m_N_GeV = 0.93141941*(N+Z) # nucleus mass in GeV
    m_N_keV = m_N_GeV*1.0e6 # nucleus mass in keV

    dRdE = zeros(shape=shape(E_r))
    FF = LabFuncs.FormFactorHelm(E_r,N+Z)**2.0
    ne = size(E_r)

    if Flux[1]>0.0:
        for i in range(0,ne):
            diff_sigma = (G_F_GeV**2.0 /(4.0*pi))*(Q_W**2.0)*m_N_GeV*(1.0 \
                        -(m_N_keV*E_r[i])/(2.0*(E_nu*1000.0)**2.0))*\
                        (0.197e-13)**2.0*(1.0e-6)*1000.0/(1.0*N+1.0*Z)*(N_A)
            diff_sigma[diff_sigma<0.0] = 0.0
            dRdE[i] = trapz(diff_sigma*Flux*FF[i],x=E_nu)
    else:
        for i in range(0,ne):
            diff_sigma = (G_F_GeV**2.0 /(4.0*pi))*(Q_W**2.0)*m_N_GeV*(1.0 \
                        -(m_N_keV*E_r[i])/(2.0*(E_nu[0]*1000.0)**2.0))*\
                        (0.197e-13)**2.0*(1.0e-6)*1000.0/(1.0*N+1.0*Z)*(N_A)
            if diff_sigma>0:
                dRdE[i] = diff_sigma*Flux[0]*E_nu[0] # for monochromatic nu's

    if sol:
        fMod = LabFuncs.EarthSunDistanceMod(t)
    else:
        fMod = 1.0

    # Convert into /ton/year/keV
    dRdE = fMod*dRdE*(365.0*3600.0*24.0*1000.0)
    return dRdE


def dRdEdO_isonu(E,E_nu,Flux,Nuc):
    E_r = sqrt(E[:,0]**2 + E[:,1]**2 + E[:,2]**2) # Recoil energy
    dR = dRdE_nu(E_r,0.0,False,E_nu,Flux,Nuc)/(4*pi)
    return dR

def dRdEdO_solarnu(E,t,E_nu,Flux,Nuc,Loc): # Directional CEnuNS for Solar
    N = Nuc.NumberOfNeutrons
    Z = Nuc.NumberOfProtons
    Q_W = N-(1-4.0*sinTheta_Wsq)*Z # weak nuclear hypercharge
    m_N_GeV = 0.93141941*(N+Z) # nucleus mass in GeV
    m_N_keV = m_N_GeV*1.0e6 # nucleus mass in keV
    m_N_kg = m_N_keV*1.783e-33
    E_nu_keV = E_nu*1e3
    x_sun = LabFuncs.SolarDirection(t,Loc)

    E_r = sqrt(E[:,0]**2 + E[:,1]**2 + E[:,2]**2) # Recoil energy
    x = zeros(shape=shape(E))
    x[:,0] = E[:,0]/E_r # Recoil direction
    x[:,1] = E[:,1]/E_r
    x[:,2] = E[:,2]/E_r
    ne =size(E_r)
    dRdEdO = zeros(shape=ne)
    cos_th_sun = -(x_sun[0]*x[:,0]+x_sun[1]*x[:,1]+x_sun[2]*x[:,2])
    FF = LabFuncs.FormFactorHelm(E_r,N+Z)**2.0


    # CHROMATIC NEUTRINOS
    if Flux[1]>0.0:
        E_max = 2*m_N_keV*E_nu_keV[-1]**2.0/(m_N_keV+E_nu_keV[-1])**2
        i_range = range(0,ne)*(E_r<=E_max)
        i_sel = i_range[i_range!=0]
        for i in i_sel:
            costh = cos_th_sun[i]
            E_nu_min = sqrt(m_N_keV*E_r[i]/2.0)
            if costh>(E_nu_min/m_N_keV):
                Eps = 1.0/(costh/E_nu_min - 1.0/m_N_keV)
                diff_sigma = (G_F_GeV**2/(4*pi))*Q_W**2*m_N_GeV*\
                            (1-(m_N_keV*E_r[i])/(2*Eps**2))*(0.197e-13)**2.0\
                            *1e-6*1000/(N+1.0*Z)*(N_A)
                Eps = Eps*(Eps>E_nu_min)
                Eps = Eps*(Eps<E_nu_keV[-1])
                F_value = interp(Eps,E_nu_keV,Flux)
                dRdEdO[i] = diff_sigma*F_value*Eps**2.0/(1000*E_nu_min)*FF[i] # /kg/keV

    # MONOCHROMATIC NEUTRINOS
    else:
        E_max = 2*m_N_keV*E_nu_keV[0]**2.0/(m_N_keV+E_nu_keV[0])**2
        i_range = range(0,ne)*(E_r<=E_max)
        i_sel = i_range[i_range!=0]
        for i in i_sel:
            costh = cos_th_sun[i]
            E_nu_min = sqrt(m_N_keV*E_r[i]/2.0)
            costh_r = ((E_nu_keV[0]+m_N_keV)/E_nu_keV[0])*sqrt(E_r[i]/(2*m_N_keV))

            # just need to accept angles close enough to costh_r to be accurate
            # around 0.01 is enough to be disguised by most angular resolutions
            if abs((costh)-(costh_r))<0.01:
                Eps = E_nu_keV[0]
                diff_sigma = (G_F_GeV**2/(4*pi))*Q_W**2*m_N_GeV*\
                            (1-(m_N_keV*E_r[i])/(2*Eps**2))*(0.197e-13)**2.0\
                            *1e-6*1000/(N+1.0*Z)*(N_A)*FF[i]
                dRdEdO[i] = diff_sigma*(Flux[0]/1000)*E_nu_keV # /kg/keV


    fMod = LabFuncs.EarthSunDistanceMod(t)
    dRdEdO = fMod*dRdEdO*3600*24*365*1000/(2*pi) # /ton/year
    return dRdEdO
