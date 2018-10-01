from numpy import pi, sqrt, exp, zeros, size, shape, array
from numpy import trapz, interp, loadtxt
from iminuit import minimize
from numpy.linalg import norm
from scipy.special import gammaln
import LabFuncs
import Params
import NeutrinoFuncs
import WIMPFuncs

#=================================Like.py======================================#
# Modeule for doing all of the Likelihood analysis
#
# Contents:
#
#===============================================================================#


#=========================Generate limits======================================#
#------------------------------CYGNUS------------------------------------------#
def CYGNUSLimit(m_vals,sigma_vals,Halo,Detector_F,Detector_He,Volume,TotTime):
    # Code for running both fluorine and helium limits

    # Put CYGNUS at Boubly
    # (this has no impact on DLs but we need to set it somewhere)
    Loc = Boulby

    # Calculate exposure for specified TPC Volume x Time
    # 1000 m^3 of SF6 at 20 torr or He at 740 torr is 0.16 tons
    Detector_F.Exposure = Volume*TotTime*(0.16/1000.0) # Convert m^3-years into ton-years
    Detector_He.Exposure = Volume*TotTime*(0.16/1000.0) # Convert m^3-years into ton-years

    # Calculate Helium limits
    DLHe = GetLimits(m_vals,sigma_vals,HaloModel,Detector_He)

    # Calculate Fluorine limits
    DLF = GetLimits(m_vals,sigma_vals,HaloModel,Detector_F)

    # Save Data
    file = open(filename,"w")
    file.write(m_vals)
    file.write(DLF)
    file.write(DLHe)
    file.close()
    print "written to:",filename
    print "----------------------------------------------------"
    return DLHe, DLF


#===================================USEFUL SUMS================================#
def lnPF(Nob,Nex): # SUM OF LOG(POISSON PDF)
    L = sum(Nob*log(Nex) - Nex - lgamma(Nob+1.0))
    return L

def lnGF(x,mu,sig): # SUM OF LOG(GAUSSIAN PDF)
    L = sum(-1.0*log(sig)-0.5*log(2.0*pi)-(x-mu)**2.0/(2.0*sig**2.0))
    return L





#===================================LIKELIHOODS================================#
#---------------------------------SIGNAL+BACKGROUND----------------------------#
def llhood1(X,N_obs,Signal,Background):
    # Signal events
    N_exp = Signal.RD*(10.0**X[0])

    # Background events
    RD_bg = Background.RD
    for i in range(0,n_bg):
    	N_exp += X[i+1]*RD_bg[:,i]

    # Likelihood
    LL = -1.0*lnPF(N_exp1,N_exp1)
    LL += -1.0*lnGF(X[1:],Background.Normalisations,Background.Uncertainties)
    return LL

#---------------------------------BACKGROUND-----------------------------------#
def llhood0(X,N_obs,Signal,Background):
    # Background events
    RD_bg = Background.RD
    N_exp = X[0]*RD_bg[:,0] # Sum over backgrounds
    for i in range(1,n_bg):
        N_exp = N_exp + X[i]*RD_bg[:,i]

    LL = -1.0*lnPF(N_obs,N_exp0)
    LL += -1.0*lnGF(X,Background.Normalisations,Background.Uncertainties)
    return llhood0






#===============================Discovery Limit================================#
def GetLimits(m_vals,sigma_vals,HaloModel,Expt):
    # Load neutrino backgrounds
    Background = GetNuFluxes(Expt.EnergyThreshold,Experiment.Gas)
    n_bg = Background.NumberOfNeutrinos
    R_bg = Background.Normalisations
    R_bg_err = Background.Uncertainties
    RD_nu = zeros(shape=(Expt.TotalNumberOfBins,n_bg))
    for i in range(0,n_bg):
        RD_nu[:,i] = BinEvents(Expt,NeutrinoFuncs.NuRate,Args=(NuBG,i))
        RD_nu[:,i] *= (1.0/R_bg[i])
    Background.RecoilDistribution(RD_nu)

    # Fix parameters for scan
    nm = size(m_vals)
    ns = size(sigma_vals)
    DL = zeros(shape=nm)
    X_in1 = zeros(shape=(n_bg+1))
    N_bg = zeros(shape=Expt.TotalNumberOfBins)
    for i in range(0,n_bg):
        N_bg = N_bg + sum(Background.RD[:,i])


    # MASS SCAN:
    for im in range(0,nm):
        Signal = WIMP(m_vals[im],1.0e-45)
        RD_wimp = BinEvents(Expt,WIMPFuncs.WIMPRate,Args=(Signal,HaloModel))
        Signal.RecoilDistribution(RD_wimp/sigma_p)

    # CROSS SECTION SCAN
    for j in range(0,ns):
        sigma_p = sigma_p_vals[j]
        N_signal = Signal.RD*sigma_p
        if sum(N_signal)>0.5: # Generally need >0.5 events to see DM
            # ------ Asimov dat -----------#
            N_obs = N_signal + N_bg

            #------ Signal + Background ------#
            X_in1[0] = log10(sigma_p)
            X_in1[1:n_bg] = R_bg
            L1 = llhood1(X_in1,Signal,Background)

            #------ Background only ------#
            X_in0 = R_bg
            step = R_bg*R_bg_err
            res = minimize(llhood0, X_in0, args=(N_obs,Signal,Background)\
                            ,options={'xtol': 1e-4, 'eps': step})

            # Test statistic
            D01 = -2.0*(L1-L0)
            if D01>9.0: # Median 3sigma detection -> D = 9
            	# Do interpolation to find discovery limit cross section
            	DL[im] = 10.0**(9.0,interp1D(array([D_prev,D01]),\
                        array([log10(s_prev),log10(sigma_p)])))
            	break
        s_prev = sigma_p # Reset for interpolation
        D_prev = D01
    return DL
