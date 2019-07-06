from numpy import *
from numpy.random import *
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
#from mpl_toolkits.basemap import Basemap
import matplotlib.gridspec as gridspec
from matplotlib.colors import ListedColormap
from matplotlib import colors
import pandas
from scipy.spatial import Delaunay
from scipy.spatial import ConvexHull
from scipy.stats import zscore,chi2,multivariate_normal
from sklearn import mixture
from scipy.special import erfinv
from scipy.stats import gaussian_kde
import cartopy.crs as ccrs
import matplotlib.ticker as mticker

# Galpy
from galpy.orbit import Orbit
from galpy.potential import MWPotential2014
from mpl_toolkits.mplot3d import Axes3D
from astropy import units
from skimage import measure

Sun = array([8.122,0.0,0.005])
Sun_cyl = array([Sun[0],0.0,Sun[2]])

from astropy import units as u
from astropy.coordinates import SkyCoord, get_constellation

# cygnus_stars = array(['β','η','γ','α','γ','δ','ι','κ','ι','δ','γ','ε','ζ'])
# #cygnus_stars = ['Deneb','gamma cyg']
# nst = size(cygnus_stars)
# cyg = zeros(shape=(nst,2))
# for i in range(0,nst):
#     c = SkyCoord.from_name(cygnus_stars[i]+' Cyg').galactic
#     cyg[i,:] = array([c.l.degree,c.b.degree])


def reverse_colourmap(cmap, name = 'my_cmap_r'):
    """
    In: 
    cmap, name 
    Out:
    my_cmap_r

    Explanation:
    t[0] goes from 0 to 1
    row i:   x  y0  y1 -> t[0] t[1] t[2]
                   /
                  /
    row i+1: x  y0  y1 -> t[n] t[1] t[2]

    so the inverse should do the same:
    row i+1: x  y1  y0 -> 1-t[0] t[2] t[1]
                   /
                  /
    row i:   x  y1  y0 -> 1-t[n] t[2] t[1]
    """        
    reverse = []
    k = []   

    for key in cmap._segmentdata:    
        k.append(key)
        channel = cmap._segmentdata[key]
        data = []

        for t in channel:                    
            data.append((1-t[0],t[2],t[1]))            
        reverse.append(sorted(data))    

    LinearL = dict(zip(k,reverse))
    my_cmap_r = mpl.colors.LinearSegmentedColormap(name, LinearL) 
    return my_cmap_r

        
limSI_face = [0.0, 0.62, 0.38]
limSI_edge = "darkgreen"
limSD_face = [0.66, 0.11, 0.03]
limSD_edge = [0.5, 0.0, 0.13]

froot = "../data/CYGNUS"


def getlimSI(Vstr,fname,RS,maxy):
    dat = loadtxt(froot+Vstr+'-'+fname+'.txt')
    m = dat[0,:]
    F = dat[1,:]
    He = dat[2,:]
    He[He==0] = 1.0e-34
    F[F==0] = 1.0e-34
    sig1 = minimum(He,F)
    sig1[sig1>=maxy] = nan
    return m,sig1/RS

def getlimSD(Vstr,fname,RS,maxy):
    dat = loadtxt(froot+Vstr+'-'+fname+'.txt')
    m = dat[0,:]
    F = dat[1,:]
    sig1 = F
    sig1[sig1==0] = 1.0
    sig1[sig1>=maxy] = nan
    return m,sig1/RS


def get1wimpSI(Vstr,fname,RS,maxy):
    dat = loadtxt(froot+Vstr+'-'+fname+'.txt')
    m = dat[0,:]
    F = dat[3,:]
    He = dat[4,:]
    He[He==0] = 1.0e-34
    F[F==0] = 1.0e-34
    sig1 = minimum(He,F)
    sig1[sig1>=maxy] = nan
    return m,sig1/RS

def get1wimpSD(Vstr,fname,RS,maxy):
    dat = loadtxt(froot+Vstr+'-'+fname+'.txt')
    m = dat[0,:]
    F = dat[3,:]
    sig1 = F
    sig1[sig1==0] = 1.0
    sig1[sig1>=maxy] = nan
    return m,sig1/RS
    
    
def MakeLimitPlot_SI(NuFloor=True,NuLabels=False,Annotations=True,Collected=False,\
                     xmin=0.1,xmax=1.0e4,ymax=1.0e-36,ymin=1.0e-51):
    plt.rcParams['axes.linewidth'] = 2.5
    fig = plt.figure(figsize=(13,12))
    ax = fig.add_subplot(111)
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')

    if NuFloor==True:
        # Xe and F neutrino floors
        nu_Xe = loadtxt("../data/WIMPLimits/SI/nufloor-Xe.txt")
        nu_F = loadtxt("../data/WIMPLimits/SI/nufloor-F.txt")
        plt.fill_between(nu_F[:,0], nu_F[:,1], y2=ymin,facecolor='dimgrey',alpha=0.4)
        plt.fill_between(nu_Xe[:,0], nu_Xe[:,1], y2=ymin,facecolor='dimgrey',alpha=0.4)
        #plt.text(0.14,8e-46,r"$\nu$-{\bf floor}: F",fontsize=19,rotation=22,color="k")
        #plt.text(0.11,2e-45,r"$\nu$-{\bf floor}: He",fontsize=19,rotation=22,color="k")

        if NuLabels:
            plt.text(0.5,5e-48,r"{\bf Solar}-$\nu$",fontsize=25)
            plt.text(18.0,0.3e-50,r"{\bf SN}-$\nu$",fontsize=25)
            plt.text(700.0,1e-49,r"{\bf Atm}-$\nu$",fontsize=25)

        nu_F = loadtxt("../data/WIMPLimits/SI/nuFloorSimple-F.txt")
        #plt.plot(nu_F[0,:],nu_F[1,:],'b-')

        nu_F = loadtxt("../data/WIMPLimits/SI/nuFloorSimple-Xe.txt")
        #plt.plot(nu_F[0,:],nu_F[1,:],'b-')

        nu_He = loadtxt("../data/WIMPLimits/SI/nuFloorSimple-He.txt")
        plt.fill_between(nu_He[0,:], nu_He[1,:], y2=ymin,facecolor='dimgrey',alpha=0.4)

    if Collected==False:
        if Annotations:
            plt.text(0.14,2.0e-38,r"{\bf CRESST}",color=[0.8, 0.25, 0.33],fontsize=20,rotation=0)
            plt.text(0.35,1.0e-39,r"{\bf CDMSlite}",color="blue",fontsize=20,rotation=0)
            plt.text(1.5,3.0e-42,r"{\bf DarkSide}",color="green",fontsize=20,rotation=-35)
            plt.text(2300.0,5.0e-44,r"{\bf PandaX}",color="navy",fontsize=20,rotation=16)
            plt.text(2300.0,1.5e-44,r"{\bf Xenon1T}",color="darkgreen",fontsize=20,rotation=16)
            plt.text(10.0,1.0e-43,r"{\bf EDELWEISS}",color=[0.67, 0.31, 0.32],fontsize=13,rotation=-20)
            plt.text(2000.0,2.0e-42,r"{\bf PICO60}",color=[0.5, 0.0, 0.13],fontsize=20,rotation=18)
            plt.text(2000.0,2.3e-40,r"{\bf PICO2L}",color=[0.5, 0.0, 0.13],fontsize=20,rotation=18)
            plt.text(22.0,1e-39,r"{\bf DAMA}",color='darkgreen',fontsize=20)
            plt.text(1200.0,5.0e-41,r"{\bf COSINE-100}",color="gold",fontsize=20,rotation=18)
            plt.arrow(0.48, 0.79, -0.04, -0.04, transform=ax.transAxes,
                  length_includes_head=True,
                  head_width=0.01, head_length=0.01, overhang=0.4,
                  edgecolor='k',facecolor='darkgreen')
            plt.arrow(0.54, 0.79, 0.04, -0.08, transform=ax.transAxes,
                  length_includes_head=True,
                  head_width=0.01, head_length=0.01, overhang=0.4,
                  edgecolor='k',facecolor='darkgreen')

        # Expt limits
        CRESST = loadtxt("../data/WIMPLimits/SI/CRESST.txt")
        plt.fill_between(CRESST[:,0], CRESST[:,1],edgecolor=None,y2=ymax,facecolor=[0.74, 0.56, 0.56])
        plt.plot(CRESST[:,0], CRESST[:,1],color=[0.8, 0.25, 0.33],linewidth=3)

        CDMSLite = loadtxt("../data/WIMPLimits/SI/CDMSLite.txt")
        plt.fill_between(CDMSLite[:,0], CDMSLite[:,1],edgecolor=None,y2=ymax,facecolor=[0.27, 0.51, 0.71])
        plt.plot(CDMSLite[:,0], CDMSLite[:,1],color="blue",linewidth=3)

        DarkSide = loadtxt("../data/WIMPLimits/SI/DarkSide.txt")
        plt.fill_between(DarkSide[:,0], DarkSide[:,1],edgecolor=None,y2=ymax,facecolor="forestgreen",alpha=0.5)
        plt.plot(DarkSide[:,0], DarkSide[:,1],color="green",linewidth=3)

        PandaX = loadtxt("../data/WIMPLimits/SI/PandaX.txt")
        plt.fill_between(PandaX[:,0], PandaX[:,1],edgecolor=None,y2=ymax,facecolor="teal")
        plt.plot(PandaX[:,0], PandaX[:,1],color="navy",linewidth=3)

        XENON1T = loadtxt("../data/WIMPLimits/SI/XENON1T.txt")
        plt.fill_between(XENON1T[:,0], XENON1T[:,1],edgecolor=None,y2=ymax,facecolor=limSI_face)
        plt.plot(XENON1T[:,0], XENON1T[:,1],color=limSI_edge,linewidth=3)
        

        EDELWEISS = loadtxt("../data/WIMPLimits/SI/EDELWEISS.txt")
        plt.plot(EDELWEISS[:,0], EDELWEISS[:,1],color=[0.67, 0.31, 0.32],linewidth=3)

        PICO60 = loadtxt("../data/WIMPLimits/SI/PICO60.txt")
        plt.plot(PICO60[:,0], PICO60[:,1],color=[0.5, 0.0, 0.13],linewidth=3)

        PICO2L = loadtxt("../data/WIMPLimits/SI/PICO2L.txt")
        plt.plot(PICO2L[:,0], PICO2L[:,1],color=[0.5, 0.0, 0.13],linewidth=3)


        DAMA1 = loadtxt("../data/WIMPLimits/SI/DAMA1.txt")
        DAMA2 = loadtxt("../data/WIMPLimits/SI/DAMA2.txt")
        plt.fill_between(DAMA1[:,0], DAMA1[:,1],edgecolor=None,y2=1.0e-45,facecolor='forestgreen')
        plt.fill_between(DAMA2[:,0], DAMA2[:,1],edgecolor=None,y2=1.0e-45,facecolor='forestgreen')
        plt.plot(DAMA1[:,0], DAMA1[:,1],color='darkgreen',linewidth=3)
        plt.plot(DAMA2[:,0], DAMA2[:,1],color='darkgreen',linewidth=3)


        COSINE = loadtxt("../data/WIMPLimits/SI/COSINE-100.txt")
        plt.plot(COSINE[:,0], COSINE[:,1],color="gold",linewidth=3)
   
    else:
        # Expt limits
        AllLimits = loadtxt("../data/WIMPLimits/SI/AllLimits-2018.txt")
        plt.fill_between(AllLimits[:,0], AllLimits[:,1],edgecolor=None,y2=ymax,facecolor=limSI_face)
        plt.plot(AllLimits[:,0], AllLimits[:,1],color=limSI_edge,linewidth=3)

    
    
    
    # Labels
    plt.xlabel(r"WIMP mass [GeV$/c^2$]",fontsize=35)
    plt.ylabel(r"SI WIMP-proton cross section [cm$^2$]",fontsize=35)
    plt.yscale('log')
    plt.xscale('log')
    ax.set_xlim(left=xmin, right=xmax)
    ax.set_ylim(bottom=ymin, top=ymax)
    plt.xticks(fontsize=25)
    plt.yticks(10.0**arange(log10(ymin),log10(ymax),1),fontsize=25)
    ax.tick_params(which='major',direction='in',width=2,length=13,right=True,top=True,pad=7)
    ax.tick_params(which='minor',direction='in',width=1,length=10,right=True,top=True)

    return fig,ax


def MakeLimitPlot_SD(NuFloor=True,NuLabels=False,Annotations=True,Collected=False,\
                     xmin = 0.1,xmax = 1.0e4,ymax = 1.0e-34,ymin = 1.0e-48):
    plt.rcParams['axes.linewidth'] = 2.5
    fig = plt.figure(figsize=(13,12))
    ax = fig.add_subplot(111)
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')

    if NuFloor==True: 
        # Xe and F neutrino floors
        nu_Xe = loadtxt("../data/WIMPLimits/SD/nufloor-Xe.txt")
        nu_F = loadtxt("../data/WIMPLimits/SD/nufloor-F.txt")
        plt.fill_between(nu_F[:,0], nu_F[:,1], y2=ymin,facecolor='dimgrey',alpha=0.4)
        plt.fill_between(nu_Xe[:,0], nu_Xe[:,1], y2=ymin,facecolor='dimgrey',alpha=0.4)
        #plt.text(0.12,3e-38,r"$\nu$-{\bf floor}: Xe",fontsize=20,rotation=22,color="k")
        #plt.text(0.12,1.5e-42,r"$\nu$-{\bf floor}: F",fontsize=20,rotation=25,color="k")

        if NuLabels:
            plt.text(0.5,5e-45,r"{\bf Solar}-$\nu$",fontsize=25)
            plt.text(18.0,0.2e-46,r"{\bf SN}-$\nu$",fontsize=25)
            plt.text(700.0,1e-46,r"{\bf Atm}-$\nu$",fontsize=25)      

    if Collected==False:
        if Annotations:
            plt.text(200.0,3.0e-37,r"{\bf PICASSO}",color='orange',fontsize=20,rotation=20)
            plt.text(2500.0,5e-39,r"{\bf PICO60}",color=[ 0.775 ,  0.55  ,  0.6085],fontsize=20,rotation=20)
            plt.text(2500.0,7.5e-38,r"{\bf PICO2L}",color=[ 0.775 ,  0.55  ,  0.6085],fontsize=20,rotation=20)
            plt.text(200.0,1e-38,r"{\bf COUPP}",color=[1.0, 0.55, 0.41],fontsize=20,rotation=18)
            plt.text(3000.0,7e-37,r"{\bf KIMS}",color=[0.81, 0.44, 0.69],fontsize=20,rotation=18)
            plt.text(2000.0,5e-40,r"{\bf IceCube} $\tau\bar{\tau}$",color=[0.55, 0.71, 0.0],fontsize=20,rotation=35)
            plt.text(15.0,5e-41,r"{\bf SK} $\tau\bar{\tau}$",color=[0.24, 0.71, 0.54],fontsize=20,rotation=0)


        # Expt limits
        PICASSO = loadtxt("../data/WIMPLimits/SD/PICASSO.txt")
        plt.fill_between(PICASSO[:,0], PICASSO[:,1],edgecolor=None,y2=ymax,facecolor='orange')
        plt.plot(PICASSO[:,0], PICASSO[:,1],color='orangered',linewidth=3)

        PICO60 = loadtxt("../data/WIMPLimits/SD/PICO60.txt")
        plt.fill_between(PICO60[:,0], PICO60[:,1],edgecolor=None,y2=ymax,facecolor=[0.66, 0.11, 0.03])
        plt.plot(PICO60[:,0], PICO60[:,1],color=[0.5, 0.0, 0.13],linewidth=3)

        PICO2L = loadtxt("../data/WIMPLimits/SD/PICO2L.txt")
        plt.plot(PICO2L[:,0], PICO2L[:,1],color=[0.5, 0.0, 0.13],linewidth=3)

        COUPP = loadtxt("../data/WIMPLimits/SD/COUPP.txt")
        plt.plot(COUPP[:,0], COUPP[:,1],color=[1.0, 0.55, 0.41],linewidth=3)

        KIMS = loadtxt("../data/WIMPLimits/SD/KIMS.txt")
        plt.plot(KIMS[:,0], KIMS[:,1],color=[0.81, 0.44, 0.69],linewidth=3)

        IceCube = loadtxt("../data/WIMPLimits/SD/IceCube-tt.txt")
        plt.plot(IceCube[:,0], IceCube[:,1],'--',color=[0.55, 0.71, 0.0],linewidth=3)

        SK = loadtxt("../data/WIMPLimits/SD/SuperK-tt.txt")
        plt.plot(SK[:,0], SK[:,1],'--',color=[0.24, 0.71, 0.54],linewidth=3)
   
    else:
        # Expt limits
        AllLimits = loadtxt("../data/WIMPLimits/SD/AllLimits-2018.txt")
        plt.fill_between(AllLimits[:,0], AllLimits[:,1],edgecolor=None,y2=ymax,facecolor=limSD_face)
        plt.plot(AllLimits[:,0], AllLimits[:,1],color=limSD_edge,linewidth=3)

    
    
    
    # Labels
    plt.xlabel(r"WIMP mass [GeV$/c^2$]",fontsize=35)
    plt.ylabel(r"SD WIMP-proton cross section [cm$^2$]",fontsize=35)
    plt.yscale('log')
    plt.xscale('log')
    ax.set_xlim(left=xmin, right=xmax)
    ax.set_ylim(bottom=ymin, top=ymax)
    plt.xticks(fontsize=25)
    plt.yticks(10.0**arange(log10(ymin),log10(ymax),1),fontsize=25)
    ax.tick_params(which='major',direction='in',width=2,length=13,right=True,top=True,pad=7)
    ax.tick_params(which='minor',direction='in',width=1,length=10,right=True,top=True)

    return fig,ax



def MySquarePlot(xlab='',ylab='',\
                 lw=2.5,lfs=45,tfs=25,size_x=13,size_y=12,Grid=False):
    plt.rcParams['axes.linewidth'] = lw
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif',size=tfs)
    
    fig = plt.figure(figsize=(size_x,size_y))
    ax = fig.add_subplot(111)
    
    ax.set_xlabel(xlab,fontsize=lfs)
    ax.set_ylabel(ylab,fontsize=lfs) 
    
    ax.tick_params(which='major',direction='in',width=2,length=13,right=True,top=True,pad=7)
    ax.tick_params(which='minor',direction='in',width=1,length=10,right=True,top=True)
    if Grid:
        ax.grid()
    return fig,ax

def MyDoublePlot(xlab1='',ylab1='',xlab2='',ylab2='',\
                 wspace=0.25,lw=2.5,lfs=45,tfs=25,size_x=20,size_y=11,Grid=False):
    plt.rcParams['axes.linewidth'] = lw
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif',size=tfs)
    
    fig, axarr = plt.subplots(1, 2,figsize=(size_x,size_y))
    gs = gridspec.GridSpec(1, 2)
    gs.update(wspace=wspace)
    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])
    ax1.tick_params(which='major',direction='in',width=2,length=13,right=True,top=True,pad=7)
    ax1.tick_params(which='minor',direction='in',width=1,length=10,right=True,top=True)
    ax2.tick_params(which='major',direction='in',width=2,length=13,right=True,top=True,pad=7)
    ax2.tick_params(which='minor',direction='in',width=1,length=10,right=True,top=True)
    
    ax1.set_xlabel(xlab1,fontsize=lfs)
    ax1.set_ylabel(ylab1,fontsize=lfs) 
    
    ax2.set_xlabel(xlab2,fontsize=lfs)
    ax2.set_ylabel(ylab2,fontsize=lfs) 
    
    if Grid:
        ax1.grid()
        ax2.grid()
    return fig,ax1,ax2


def MyTriplePlot(xlab1='',ylab1='',xlab2='',ylab2='',xlab3='',ylab3='',\
                 wspace=0.25,lw=2.5,lfs=45,tfs=25,size_x=20,size_y=7,Grid=False):
    plt.rcParams['axes.linewidth'] = lw
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif',size=tfs)
    
    fig, axarr = plt.subplots(1, 3,figsize=(size_x,size_y))
    gs = gridspec.GridSpec(1, 3)
    gs.update(wspace=wspace)
    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])
    ax3 = plt.subplot(gs[2])
    
    ax1.tick_params(which='major',direction='in',width=2,length=13,right=True,top=True,pad=7)
    ax1.tick_params(which='minor',direction='in',width=1,length=10,right=True,top=True)
    
    ax2.tick_params(which='major',direction='in',width=2,length=13,right=True,top=True,pad=7)
    ax2.tick_params(which='minor',direction='in',width=1,length=10,right=True,top=True)
    
    ax3.tick_params(which='major',direction='in',width=2,length=13,right=True,top=True,pad=7)
    ax3.tick_params(which='minor',direction='in',width=1,length=10,right=True,top=True)
    
    ax1.set_xlabel(xlab1,fontsize=lfs)
    ax1.set_ylabel(ylab1,fontsize=lfs) 
    
    ax2.set_xlabel(xlab2,fontsize=lfs)
    ax2.set_ylabel(ylab2,fontsize=lfs) 
    
    ax3.set_xlabel(xlab3,fontsize=lfs)
    ax3.set_ylabel(ylab3,fontsize=lfs) 
    
    if Grid:
        ax1.grid()
        ax2.grid()
        ax3.grid()
    return fig,ax1,ax2,ax3



def MollweideMap(ax,TH,PH,fv0,cmin,cmax,nlevels,cmap,tfs,PlotCygnus=False,gridlinecolor='k',GalacticPlane=False):
    plt.rcParams['axes.linewidth'] = 3
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif',size=15)

    
    ax.contourf(rad2deg(PH), rad2deg(TH), fv0,nlevels,transform=ccrs.PlateCarree(),cmap=cmap,vmin=cmin,vmax=cmax)
    gl = ax.gridlines(color=gridlinecolor,linewidth=1.5, linestyle='--',alpha=0.5)
    gl.ylocator = mticker.FixedLocator([-90,-60, -30, 0, 30, 60,90])
    ax.outline_patch.set_linewidth(3)
   

    tx = array([r'$-60^\circ$',r'$-30^\circ$',r'$0^\circ$',r'$+30^\circ$',r'$+60^\circ$']) 
    xtx = array([0.17,0.05,-0.01,0.05,0.18])
    ytx = array([0.08,0.26,0.49,0.72,0.9])
    
    for i in range(0,size(xtx)):
        plt.text(xtx[i],ytx[i],tx[i],transform=ax.transAxes,horizontalalignment='right',verticalalignment='center',fontsize=tfs)


    if PlotCygnus==True:
        ax.plot(-cyg[0:4,0],cyg[0:4,1],'-',color='crimson',transform=ccrs.PlateCarree())
        ax.plot(-cyg[4:,0],cyg[4:,1],'-',color='crimson',transform=ccrs.PlateCarree())
        ax.plot(-cyg[:,0],cyg[:,1],'.',color='k',ms=5,transform=ccrs.PlateCarree())

    if GalacticPlane==True:
        ax.plot([-181,181],[0,0],'-',color=gridlinecolor,lw=1.5,transform=ccrs.PlateCarree())
        ax.text(125,4,'Galactic',color=gridlinecolor,transform=ccrs.PlateCarree(),fontsize=int(tfs*0.8))
        ax.text(135,-10,'plane',color=gridlinecolor,transform=ccrs.PlateCarree(),fontsize=int(tfs*0.8))
    return





def PointScatter(xin,yin):
    dens = gaussian_kde(vstack([xin,yin]))(vstack([xin,yin]))
    idx = dens.argsort()
    x, y, dens = xin[idx], yin[idx], dens[idx]
    return x,y,dens


def in_hull(p, hull):
    if not isinstance(hull,Delaunay):
        hull = Delaunay(hull)
    return hull.find_simplex(p)>=0

def chaikins_corner_cutting(x_edge,y_edge, refinements=3):
    edge=zeros(shape=(size(x_edge)+1,2))
    edge[:-1,0] = x_edge
    edge[:-1,1] = y_edge
    edge[-1,0] = x_edge[0]
    edge[-1,1] = y_edge[0]
    coords = array(edge)
    for _ in range(refinements):
        L = coords.repeat(2, axis=0)
        R = empty_like(L)
        R[0] = L[0]
        R[2::2] = L[1:-1:2]
        R[1:-1:2] = L[2::2]
        R[-1] = L[-1]
        coords = L * 0.75 + R * 0.25

    return coords

def col_alpha(col,alpha=0.1):
    rgb = colors.colorConverter.to_rgb(col)
    bg_rgb = [1,1,1]
    return [alpha * c1 + (1 - alpha) * c2
            for (c1, c2) in zip(rgb, bg_rgb)]



def MollweideMap1(ax,TH,PH,fv0,cmin,cmax,nlevels,cmap,tfs,PlotCygnus=False,gridlinecolor='k',GalacticPlane=False):
    plt.rcParams['axes.linewidth'] = 3
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif',size=15)

    
    ax.contourf(rad2deg(PH), rad2deg(TH),fv0,nlevels,transform=ccrs.PlateCarree(),\
                cmap=cmap,vmin=cmin,vmax=cmax,linestyles='none')
    ax.contourf(rad2deg(PH), rad2deg(TH),fv0,nlevels,transform=ccrs.PlateCarree(),\
                cmap=cmap,vmin=cmin,vmax=cmax,linestyles='none')
    ax.contourf(rad2deg(PH), rad2deg(TH),fv0,nlevels,transform=ccrs.PlateCarree(),\
                cmap=cmap,vmin=cmin,vmax=cmax,linestyles='none')

    gl = ax.gridlines(color=gridlinecolor,linewidth=1.5, linestyle='--',alpha=0.5)
    gl.ylocator = mticker.FixedLocator([-90,-60, -30, 0, 30, 60,90])
    ax.outline_patch.set_linewidth(3)
   

    tx = array([r'$-60^\circ$',r'$-30^\circ$',r'$0^\circ$',r'$+30^\circ$',r'$+60^\circ$']) 
    xtx = array([0.17,0.05,-0.01,0.05,0.18])
    ytx = array([0.08,0.26,0.49,0.72,0.9])
    
    for i in range(0,size(xtx)):
        plt.text(xtx[i],ytx[i],tx[i],transform=ax.transAxes,horizontalalignment='right',verticalalignment='center',fontsize=tfs)


    if PlotCygnus==True:
        ax.plot(-cyg[0:4,0],cyg[0:4,1],'-',color='crimson',transform=ccrs.PlateCarree())
        ax.plot(-cyg[4:,0],cyg[4:,1],'-',color='crimson',transform=ccrs.PlateCarree())
        ax.plot(-cyg[:,0],cyg[:,1],'.',color='k',ms=5,transform=ccrs.PlateCarree())

    if GalacticPlane==True:
        ax.plot([-181,181],[0,0],'-',color=gridlinecolor,lw=1.5,transform=ccrs.PlateCarree())
        ax.text(125,4,'Galactic',color=gridlinecolor,transform=ccrs.PlateCarree(),fontsize=int(tfs*0.8))
        ax.text(135,-10,'plane',color=gridlinecolor,transform=ccrs.PlateCarree(),fontsize=int(tfs*0.8))
    return

