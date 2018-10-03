from __future__ import print_function
from numpy import pi, sqrt, exp, zeros, size, shape, array, append
from numpy import trapz, interp, loadtxt, log10, log, savetxt, vstack
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
def CYGNUSLimit(filename,m_vals,sigma_vals,Volume,TotTime,HaloModel=Params.SHM,ne=20,nt=1,nside=4,Eoff=False,ReadOut_Name="Ideal"):
    # Code for running both fluorine and helium limits
    print("="*50)
    print("CYGNUS limits")
    E_th_F = 3.0
    E_th_He = 1.8
    E_max = 200.0

    Loc = Params.Boulby

    # Calculate exposure for specified TPC Volume x Time
    # 1000 m^3 of SF6 at 20 torr or He at 740 torr is 0.16 tons
    Exposure = Volume*TotTime*(0.16/1000.0)

    # Set detector
    Detector_F = Params.Detector(E_th_F,E_max,Params.F19,Loc,Exposure,ne,nt,nside,Eoff,ReadOut_Name)
    Detector_He = Params.Detector(E_th_He,E_max,Params.He4,Loc,Exposure,ne,nt,nside,Eoff,ReadOut_Name)

    # Calculate Fluorine limits

    DLF = GetLimits(m_vals,sigma_vals,HaloModel,Detector_F)
    print("-"*50)

    # Calculate Helium limits
    print("-"*50)
    DLHe = GetLimits(m_vals,sigma_vals,HaloModel,Detector_He)
    print("-"*50)
    # Save Data
    savetxt(filename,vstack((m_vals,DLF,DLHe)), delimiter='\t')

    print("written to:",filename)
    print("="*50)
    return DLHe, DLF


#===================================USEFUL SUMS================================#
def lnPF(Nob,Nex): # SUM OF LOG(POISSON PDF)
    L = sum(Nob*log(Nex) - Nex - gammaln(Nob+1.0))
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
    for i in range(0,size(RD_bg,1)):
    	N_exp += X[i+1]*RD_bg[:,i]

    # Likelihood
    LL = -1.0*lnPF(N_obs,N_exp)
    LL += -1.0*lnGF(X[1:],Background.Normalisations,Background.Uncertainties)
    return LL

#---------------------------------BACKGROUND-----------------------------------#
def llhood0(X,N_obs,Signal,Background):
    # Background events
    RD_bg = Background.RD
    N_exp = X[0]*RD_bg[:,0] # Sum over backgrounds
    for i in range(1,size(RD_bg,1)):
        N_exp = N_exp + X[i]*RD_bg[:,i]

    LL = -1.0*lnPF(N_obs,N_exp)
    LL += -1.0*lnGF(X,Background.Normalisations,Background.Uncertainties)
    return LL






#===============================Discovery Limit================================#
def GetLimits(m_vals,sigma_vals,HaloModel,Expt):
    # Load neutrino backgrounds
    Background = NeutrinoFuncs.GetNuFluxes(Expt.EnergyThreshold,Expt.Nucleus)
    n_bg = Background.NumberOfNeutrinos
    R_bg = Background.Normalisations
    R_bg_err = Background.Uncertainties
    RD_nu = zeros(shape=(Expt.TotalNumberOfBins,n_bg))
    for i in range(0,n_bg):
        RD_nu[:,i] = LabFuncs.BinEvents(Expt,NeutrinoFuncs.NuRate,Background,i)
        RD_nu[:,i] *= (1.0/R_bg[i])
    Background.RecoilDistribution(RD_nu)

    # Fix parameters for scan
    nm = size(m_vals)
    ns = size(sigma_vals)
    DL = zeros(shape=nm)
    X_in1 = zeros(shape=(n_bg+1))
    N_bg = zeros(shape=Expt.TotalNumberOfBins)
    for i in range(0,n_bg):
        N_bg = N_bg + R_bg[i]*Background.RD[:,i]


    # MASS SCAN:
    for im in range(0,nm):
        Signal = Params.WIMP(m_vals[im],1.0e-45)
        RD_wimp = LabFuncs.BinEvents(Expt,WIMPFuncs.WIMPRate,Signal,HaloModel)
        Signal.RecoilDistribution(RD_wimp/1.0e-45)

        # CROSS SECTION SCAN
        for j in range(0,ns):
            sigma_p = sigma_vals[j]
            N_signal = Signal.RD*sigma_p
            D_prev = 0.0
            s_prev = sigma_p
            if sum(N_signal)>0.5: # Generally need >0.5 events to see DM
                # ------ Asimov dat -----------#
                N_obs = N_signal + N_bg

                #------ Signal + Background ------#
                X_in1 = append(log10(sigma_p),R_bg)
                L1 = llhood1(X_in1,N_obs,Signal,Background)

                #------ Background only ------#
                X_in0 = R_bg
                step = R_bg_err
                res = minimize(llhood0, X_in0, args=(N_obs,Signal,Background)\
                                ,options={'xtol': R_bg_err, 'eps': 2*step})
                L0 = res.fun
                #L0 = llhood0(X_in0,N_obs,Signal,Background)

                # Test statistic
                D01 = -2.0*(L1-L0)
                #print(j,sigma_p,D01,sum(N_signal),L1,L0)
                if D01>9.0: # Median 3sigma detection -> D = 9
                	# Do interpolation to find discovery limit cross section
                	DL[im] = 10.0**(interp(9.0,array([D_prev,D01]),\
                            array([log10(s_prev),log10(sigma_p)])))
                	break
                s_prev = sigma_p # Reset for interpolation
                D_prev = D01
        #Params.printProgressBar(im, nm)
        print("m_chi = ",m_vals[im],"sigma_p = ",DL[im],sum(N_signal),sum(N_bg))
    return DL
