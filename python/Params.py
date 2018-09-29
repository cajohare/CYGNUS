from numpy import array, sqrt, pi, exp, interp, loadtxt, zeros
from scipy.special import erf


#==============================================================================#
class WIMP:
    def __init__(self,m_chi,sigma_p):
        self.Mass = m_chi
        self.SICrossSection = sigma_p

    def RecoilDistribution(self,RD):
        self.RD = RD
#==============================================================================#





#==============================================================================#
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
He4 =   Nucleus(1.0,     2,    2,0.01,  0.000,   0.000)
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
#==============================================================================#










#==============================================================================#
nufile_root = "-1000.txt"
nufile_dir = "../neutrinos/"
n_nu_tot = 11
nuname = ["" for x in range(0,n_nu_tot)]
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
mono = zeros(n_nu_tot,dtype=bool)
mono[[1,3,4]] = True
n_Enu_vals = 1000 # All neutrino backgrounds saved with 1000 entries
# Monochromatic neutrinos (2, 4, 5) have a negative value for E_nu which is
# Used to tell the rate formula to use monochromatic result
# Neutrino max energies (MeV):
NuMaxEnergy = array([0.42341,1.44,18.765,0.3843,0.8613,16.34,1.193,\
                    1.7285,1.7365,91.201,981.75])
# Neutrino fluxes (cm-2 s-1 MeV-1) and uncertainties (%):
# (from Vinyoles et al (2017) Barcelona GS98 SSM)
NuFlux = array([5.98e10,1.44e8,7.98e3,4.93e8,4.50e9,5.16e6,\
                        2.78e8,2.05e8,5.29e6,85.7,10.54])
NuUnc = array([0.006, 0.01, 0.3,0.06, 0.06, 0.02, 0.15 ,\
                        0.17 ,0.2 ,0.5, 0.25])

class Neutrinos:
    def __init__(self,n_nu,energies,fluxes,normlisations,uncertainties):
        self.Flux = fluxes
        self.Energy = energies
        self.Uncertainties = uncertainties
        self.Normalisations = normlisations
        self.NumberOfNeutrinos = n_nu

    def RecoilDistribution(self,RD):
        self.RD = RD
#==============================================================================#








#==============================================================================#
froot_er = "-EnergyRes.txt"
fdir_er = "../readouts/energyres/"
froot_eff = "-Efficiency.txt"
fdir_eff = "../readouts/efficiency/"
froot_ar = "-AngRes.txt"
fdir_ar = "../readouts/angres/"
froot_ht = "-HeadTail.txt"
fdir_ht = "../readouts/headtail/"
class Detector:
    def __init__(self,Nuc,ne=50,nt=10,nside=4,Eoff=False,ReadOut_Name="Ideal"):
        self.Gas = Nuc

        Ebins = logspace(log10(E_th),log10(E_max),ne+1)
        self.Energies = Ebins
        self.Times = linspace(Jan1,Jan1+365.0,nt)
        self.TotalNumberOfBins = ne*nt
        if nside>0:
            if nside==2:
                npix = 12*nside**2
                x_pix[:,0],x_pix[:,1],x_pix[:,2] = loadtxt("../pixels/x_pix2.txt")
            elif nside==4:
                x_pix[:,0],x_pix[:,1],x_pix[:,2] = loadtxt("../pixels/x_pix4.txt")
            elif nside==4:
                x_pix[:,0],x_pix[:,1],x_pix[:,2] = loadtxt("../pixels/x_pix8.txt")
            self.Directions = x_pix
            self.TotalNumberOfBins = ne*nt*npix
            if Eoff:
                self.EnergyOff = True
                self.TotalNumberOfBins = ne*nt


        data_er = loadtxt(fdir_er+ReadOut_Name+froot_er)
        data_eff = loadtxt(fdir_eff+ReadOut_Name+froot_eff)
        data_ar = loadtxt(fdir_ar+ReadOut_Name+froot_ar)
        data_ht = loadtxt(fdir_ht+ReadOut_Name+froot_ht)

        if Nuc.NumberOfProtons==19:
            icol = 1
        elif Nuc.NumberOfProtons==2:
            icol = 2

        self.EnergyRes = interp(Ebins,data_er[:,0],data[:,icol])
        self.Efficiency = interp(Ebins,data_eff[:,0],data[:,icol])
        self.AngRes = interp(Ebins,data_ar[:,0],data[:,icol])
        self.HeadTail = interp(Ebins,data_ht[:,0],data[:,icol])

#==============================================================================#




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
