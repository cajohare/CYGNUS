from numpy import pi, sqrt, exp, zeros, size, shape, array
from numpy import trapz, interp, loadtxt
from numpy.linalg import norm
import LabFuncs
import Params

#================================NeutrinoFuncs=================================#
# Contents:
#
#==============================================================================#
# ORDER OF NEUTRINOS
# Neutrno types and file locations
nufile_root = "-1000.txt"
nufile_dir = "../neutrinos/"
nuname = ["" for x in range(0,11)]
nuname[0] = "pp"
nuname[1] = "pep"
nuname[2] = "hep"
nuname[3] = "7Be1"
nuname[4] = "7Be2"
nuname[5] = "8B"
nuname[6] = "13N"
nuname[7] = "15O"
nuname[8] = "17F"
nuname[9] = "DSNB"
nuname[10] = "Atm"
mono = zeros(11,dtype=bool)
mono[[1,3,4]] = True

#==============================================================================#
# Neutrino max energies (MeV):
NuMaxEnergy = array([0.42341,1.44,18.765,0.3843,0.8613,16.34,1.193,\
                    1.7285,1.7365,91.201,981.75])


# Neutrino fluxes (cm-2 s-1 MeV-1) and uncertainties (%):
# (from Vinyoles et al (2017) Barcelona GS98 SSM)
NuFlux_B17GS98 = array([5.98e10,1.44e8,7.98e3,4.93e8,4.50e9,5.16e6,\
                        2.78e8,2.05e8,5.29e6,85.7,10.54])
NuUnc_B17GS98 = array([0.006, 0.01, 0.3,0.06, 0.06, 0.02, 0.15 ,\
                        0.17 ,0.2 ,0.5, 0.25])

N_A = 6.02214e23 # Avocado's constant
sinTheta_Wsq = 0.2387e0 # sin(Theta_W) weak mixing angle
G_F_GeV = 1.16637e-5 # GeV**-2 ! Fermi constan in GeV

#========================================Neutrino data========================================#
def GetNuFluxes(E_th,Nuc):
    # Reads each neutrino flux data file
    # Each flux file has 1000 rows apart from monochromatic ones
    # the energies are stored in E_nu_all, fluxes in Flux_all

    # CURRENTLY USING BARCELONA FLUXES
    NuFlux = NuFlux_B17GS98
    NuUnc = NuUnc_B17GS98
    nvals = 1000 # All neutrino backgrounds saved with 1000 entries

    # Monochromatic neutrinos (2, 4, 5) have a negative value for E_nu which is
    # Used to tell the rate formula to use monochromatic result

    # Figure out which backgrounds give recoils above E_th
    E_r_max = MaxNuRecoilEnergies(Nuc) # Max recoil energy for neutrino
    n_nu = size(E_r_max>E_th)
    sel = range(1,12)*(E_r_max>E_th)
    sel = sel[sel!=0]-1
    E_nu_all = zeros(shape=(1000,n_nu))
    Flux_all = zeros(shape=(1000,n_nu))
    Flux_err = zeros(shape=(1000,n_nu))
    # Load in all backgrounds with sel #= 0
    ii = 0
    for s in sel:
        if mono[s]:
            E_nu_all[0,ii] = NuMaxEnergy[s]
            Flux_all[0,ii] = NuFlux[s]
        else:
            data = loadtxt(nufile_dir+nuname[s]+nufile_root,delimiter=',')
            E_nu_all[:,ii],Flux_all[:,ii] = data[:,0],data[:,1]
            Flux_all[:,ii] = Flux_all[:,ii]*NuFlux[s]
        Flux_err[ii] = NuUnc[s] # Select rate normalisation uncertainties
        ii = ii+1
    NuBG = Params.Neutrinos(n_nu,E_nu_all,Flux_all,Flux_err)
    return NuBG

 #-----------------------------------------------------------------------------#
def MaxNuRecoilEnergies(Nuc): # Max recoil energies
    m_N = 0.93141941*(Nuc.MassNumber)*1.0e6
    E_r_max = 2*m_N*(1000.0*NuMaxEnergy)**2.0/(m_N+1000*NuMaxEnergy)**2.0
    return E_r_max


#===================================nu spectra=================================#
def NuRate_Energy(E_r,t,NuBG,Nuc):
    E_nu_all = NuBG.Energy
    Flux_all = NuBG.Flux
    n_nu = NuBG.NumberOfNeutrinos

    ne = size(E_r)
    nt = size(t)
    dR = zeros(shape=(nt,ne,n_nu))
    fMod = LabFuncs.EarthSunDistanceMod(t)
    for nu in range(0,n_nu):
        if nu>=(n_nu-1):
            fMod[:] = 1.0

        # Compute rate
        dRdE = dRdE_nu(E_r,E_nu_all[:,nu],Flux_all[:,nu],Nuc) # correct for form factor
        for i in range(0,nt):
            dR[i,:,nu] =  fMod[i]*dRdE

    return dR

def dRdE_nu(E_r,E_nu,Flux,Nuc):
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

    # Convert into /ton/year/keV
    dRdE = dRdE*(365.0*3600.0*24.0*1000.0)
    return dRdE

def NuRate_Direction(E,t,NuBG,Nuc,Loc):
    E_nu_all = NuBG.Energy
    Flux_all = NuBG.Flux
    n_nu = NuBG.NumberOfNeutrinos

    ne = size(E_r)/3
    nt = size(t)
    dR = zeros(shape=(nt,ne,n_nu))
    fMod = LabFuncs.EarthSunDistanceMod(t)
    for nu in range(0,n_nu):
        if nu>=(n_nu-1):
            fMod[:] = 1.0

        # Compute rate
        dRdEdO = dRdEdO_nu(E,t,E_nu,Flux,Nuc,Loc) # correct for form factor
        for i in range(0,nt):
            dR[i,:,nu] =  fMod[i]*dRdEdO

    return dR

def dRdEdO_nu(E,t,E_nu,Flux,Nuc,Loc): # Directional CEnuNS for Solar
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


    dRdEdO = dRdEdO*3600*24*365*1000/(2*pi) # /ton/year
    #dRdEdO[dRdEdO==NaN]=0.0

    return dRdEdO
