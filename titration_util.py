#!/usr/bin/python
import sys
import numpy as np
import os

import matplotlib
#matplotlib.use("TkAgg")
from matplotlib import pyplot as plt

import scipy.stats as stats
import scipy.optimize as op
import glob

from emcee_tools import helper_emcee as he
from emcee_tools import chain_plots as cp

from co2sys import co2sys
#from co2sys.cc_py_mod import py_interface_mod as cc


rho_a = 1.000
S = 0.05

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def get_data_from_file( filename):
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ''' Reads data from PC_LIMS reports
    returns
        mass_s  mass of sample
        conc_a  concentration of acid
        data_array
            vol, pH, temp, time
    '''
    fp = open(filename, errors='ignore')
    lines = fp.readlines()

    # find the mass of sample and the concentration of the acid
    for i,line in enumerate(lines):
        ll = line.rstrip().split('\t')

        if(ll[0]=='$S Sample data V1'):
            next_line=lines[i+1].split('\t')
            mass_s=float(next_line[2])

        if(ll[0]=='$S Titrant1 V1'):
            next_line=lines[i+1].split('\t')
            conc_a=float(next_line[1])*float(next_line[3])

    time=[]
    vol =[]
    pH  =[]
    temp=[]
    titr_output=False
    for i,line in enumerate(lines):
        ll = line.rstrip().split('\t')
        if( ll[0]=='$E'):
           titr_output=False

        if( titr_output ):
            vol.append (ll[1])
            pH.append  (ll[2])
            temp.append(ll[5])
            time.append(ll[4])
            

        if(ll[0]=='$S Mode 2' and ll[2]=='DET pH' and ll[3]=='V1.0'):
            titr_output=True

    fp.close()
    return mass_s, conc_a,np.array([vol,pH,temp,time],dtype=float)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def get_multi_data( filenames, S=33, const=10):
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ''' sets up a array containg the data for a set files
    last column of vol_data contians the index i fomr the filenamed
    '''
    nf = len(filenames)

    aux_data = np.zeros([nf,4])        # contains the mass_s, conc_a, S
    vol_data = None
    for i,filename in enumerate(filenames):
        mass_s, conc_a, data = get_data_from_file(filename)

        nc,no = np.shape(data) #
        datat = np.append(data,[np.zeros(no)+i],axis=0)
        if( vol_data is None ):
            vol_data = datat
        else:
            vol_data = np.append(vol_data,datat,axis=1)
        aux_data[i,:] = [mass_s, conc_a, S, const ]
    return vol_data, aux_data


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def calc_curve( theta, data, aux ):
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ''' simulate the pH for a titration 
        theta :
            pH slope correcrtion 
            pH  offset
            K1_f, K2_f multiplicative correction for equlibirum coefficents
            DIC, TA    in mM  ( just to check scales ~1)
        data  :
            [:,0]     vol of acid in mL
            [:,2]     Temperatiuree in deg C
        args
            optional arguements
    
    '''
    # uppack the data

    a_pH,b_pH, K1_f,K2_f = theta[0:4]
    DIC, TA              = theta[4:6]*1e-3      # from mM to M

    mass_s, conc_a, S, const = aux[0:4]

    nd,no = np.shape(data)      

    pHs = np.zeros(no)

    pH0=3.8 # initial guess at pH, low as beyond endpoint at 4.5
    for i in reversed(range(no)):
        vol= data[0,i]
        TK = data[2,i]+273.15
        mass_a      =   vol*rho_a
        DIC_i       =    mass_s*DIC                 /(mass_s + mass_a )
        TA_i        =  ( mass_s*TA - mass_a*conc_a )/(mass_s + mass_a )
        S_i         =  ( mass_s*S  + mass_a*0.     )/(mass_s + mass_a )
        #pHs[i],outp = co2sys.CC_solve( DIC_i, TA_i, TK, S_i, const=10, K1_f=K1_f, K2_f=K2_f, pHi=pH0)
        pHs[i] = co2sys.CC_solve_pH( DIC_i, TA_i, TK, S_i, const=const, TP=36.e-6,TSi=0.,TB=None,TS=None,TF=None,K1_f=K1_f, K2_f=K2_f, pHi=pH0)
        pH0         = pHs[i]
    pHs  = a_pH*(pHs - 7.)+7. + b_pH
    return pHs

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def calc_curve_arr( theta, data, aux ):
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ''' simulate the pH for a titration 
        theta :
            pH slope correcrtion 
            pH  offset
            K1_f, K2_f multiplicative correction for equlibirum coefficents
            DIC, TA    in mM  ( just to check scales ~1)
        data  :
            [:,0]     vol of acid in mL
            [:,2]     Temperatiuree in deg C
        args
            optional arguements
    
    '''
    # uppack the data
    TP,TSi,TNH3,THS2 = 0
    a_pH,b_pH, K1_f,K2_f = theta[0:4]
    DIC, TA              = theta[4:6]*1e-3      # from mM to M

    mass_s, conc_a, S, const = aux[0:4]

    nd,no = np.shape(data)

    pHs = np.zeros(no)

    pH0=3.8 # initial guess at pH, low as beyond endpoint at 4.5
    for i in reversed(range(no)):
        vol= data[0,i]
        TK = data[2,i]+273.15
        mass_a      =   vol*rho_a
        DIC_i       =    mass_s*DIC                 /(mass_s + mass_a )
        TA_i        =  ( mass_s*TA - mass_a*conc_a )/(mass_s + mass_a )
        S_i         =  ( mass_s*S  + mass_a*0.     )/(mass_s + mass_a )
        #pHs[i],outp = co2sys.CC_solve( DIC_i, TA_i, TK, S_i, const=10, K1_f=K1_f, K2_f=K2_f, pHi=pH0)
        #pHs[i] = co2sys.CC_solve_pH( DIC_i, TA_i, TK, S_i, const=const, TP=36.e-6,TSi=0.,TB=None,TS=None,TF=None,K1_f=K1_f, K2_f=K2_f, pHi=pH0)
        pHs[i] = coc2sys.CC_solve_pH_arr( S,TK,P, TP,TSi,TNH3,TH2S, TA,DIC, const,K1_f,K2_f, pH0 )
        pH0         = pHs[i]
    pHs  = a_pH*(pHs - 7.)+7. + b_pH
    return pHs

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def log_like( theta, data, *args ):
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    sg2  = theta[-1]
    pHo  = data[1]

    pHc  = calc_curve(theta,data, *args)
    #pHc   = cc.calc_curve( theta,data,args)
    return sum(stats.norm.logpdf( pHo-pHc, scale=sg2 ))

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def unpack_theta( theta, i ):
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    theta_i  = np.zeros(6)

    ndim = np.ndim(theta)
    if( ndim==1 ):
        theta_i[0:4] = theta[0:4]
        theta_i[4:6] = theta[4+2*i:4+2*(i+1)]
    else:
        theta_i = theta[i]
    return theta_i

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def log_like_m( theta, vol_data, aux_data ):
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    sg2  = theta[-1]

    ll = 0.0
    for i,adata in enumerate(aux_data):
        theta_i = unpack_theta( theta, i )
        vdata   = vol_data[:,vol_data[4]==i]
        pHo     = vdata[1]
        pHc     = calc_curve( theta_i, vdata, adata )
        #pHc  = cc.calc_curve( theta_i, vdata, adata )
        ll += sum(stats.norm.logpdf( pHo-pHc, scale=sg2 ))
    return ll

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def prior_list( th, i):
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if( i==0 ): return stats.lognorm.logpdf(th,s=0.1)
    if( i==1 ): return stats.norm.   logpdf(th,scale=0.1)
    if( i==2 ): return stats.lognorm.logpdf(th,s=0.4)
    if( i==3 ): return stats.lognorm.logpdf(th,s=0.4)
    return 0
    

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def prior( theta ):
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    return np.sum([ prior_list(theta[i],i) for i in range(len(theta))])

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def log_prob(theta, *args):
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    lprior = prior(theta)
    if( np.isfinite(lprior) ):
        llike = log_like_m( theta, *args )
        if( np.isfinite(llike) ):
            return lprior + llike
        else:
            return -np.inf
    else:
        return -np.inf

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def nll( theta, *args ):
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    return -log_like_m( theta, *args)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def plot_fit( vol_data, aux_data, theta=None, xaxis="vol" ):
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    fig,ax = plt.subplots(1,1)
    ncurve = np.shape(aux_data)[0]
    colormap = plt.cm.gist_ncar

    if( xaxis == "vol"):
        xfact = 1.0
        ax.set_xlabel("Vol acid (mL)")
    elif( xaxis=="alk"):
        ax.set_xlabel(r'Alkalinity ($\mu$M)')
    ax.set_ylabel("pH")


    color=iter(plt.cm.rainbow(np.linspace(0,1,ncurve)))
    for i,aux in enumerate(aux_data):
        c=next(color)
        # select this run's data
        vol = vol_data[:,vol_data[4]==i]
        mass_s, conc_a = aux[0:2]
        if( xaxis == "alk"):
            xfact=conc_a/mass_s*1e6

        plt.plot( vol[0]*xfact, vol[1],'o',color=c)

        if( theta is not None ):
            theta_i = unpack_theta( theta, i )
            #pH_m = cc.calc_curve( theta_i, vol, [mass_s, conc_a, S] )
            pH_m  = calc_curve(theta_i,vol, aux )
            plt.plot( vol[0]*xfact, pH_m,'-',color=c )
    plt.show()
    return

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def get_theta_from_file( fname, th_file):
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    theta=[]
    try:
        names= np.genfromtxt(th_file, usecols=[0],dtype="U")
    except:
        return theta
    im = np.where( names==fname )
    if( len(im) != 0 and len(im[0])!=0 ):
        irow = im[0][-1]
        thetas = np.genfromtxt( th_file, usecols=range(1,8))
        theta = thetas[irow]
    return  theta

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def gen_theta_guess(  vol_data, aux_data ):
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    nf,na = np.shape(aux_data)
    theta=np.zeros( 4+nf*2 + 1)
    theta[0:4] = [ 0.9,0.01, 1.09,1.09 ]
    for i in range(nf):
        vdata = vol_data[:,vol_data[4]==i]
        mass_s, conc_a, S, const = aux_data[i,0:4]

        ind = np.abs(vdata[1]-4.35).argmin()
        TA = vdata[0,ind]/mass_s*conc_a

        ind = np.abs(vdata[1]-8.05).argmin()
        DIC = TA-vdata[0,ind]/mass_s*conc_a
        theta[4+2*i:4+2*(i+1)] = [DIC*1e3,TA*1e3]
    theta[-1] = 0.05
    return theta


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def plot_DIC_Alk( chain, filenames=None ):
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    import datetime as dt
    nw,ns,nt = np.shape(chain)
    nr = int((nt-5)/2)
    if( filenames is not None ):
        rel_times = np.zeros(nr)
        t0 = dt.datetime(2017,9,12,00,00)
        for i,filename in enumerate(filenames):
            t1=dt.datetime.strptime( filename[-19:-4], "%Y%m%d-%H%M%S")
            rel_times[i] = (t1-t0)/dt.timedelta(hours=1)
    else:
        rel_times  = np.arange(int((nt-5)/2)+1)
    for i in range(nr):
        t = rel_times[i]
        ii=4+2*i
        percs=np.percentile( chain[:,:,ii].flatten(), [2.5,25,50,75,97.5])*1000
        plt.errorbar( t, percs[2], yerr=[[percs[2]-percs[1]],[percs[3]-percs[2]]], capsize=0,lw=4, color="red" )
        plt.errorbar( t, percs[2], yerr=[[percs[2]-percs[0]],[percs[4]-percs[2]]], capsize=0,lw=2, color="red" )
       
        percs=np.percentile( chain[:,:,ii+1].flatten(), [2.5,25,50,75,97.5])*1000
        plt.errorbar( t, percs[2], yerr=[[percs[2]-percs[1]],[percs[3]-percs[2]]], capsize=0,lw=4, color="blue" )
        plt.errorbar( t, percs[2], yerr=[[percs[2]-percs[0]],[percs[4]-percs[2]]], capsize=0,lw=2, color="blue" )
    plt.axhline(1780.,ls='--',color='red')
    plt.xlabel("Time (hours)")
    plt.ylabel(r'Concentraion ($\mu  M$)')
    plt.grid()
    plt.show()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -                                        
def save_file( filenames, theta, f ):
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    for i,filename in enumerate(filenames):
        print(("{} "+(" {:12.4f}"*7)).format(filename,*theta[0:4],              \
                                                      *theta[4+2*i:4+2*(i+1)],  \
                                                       theta[-1] ), file=f)
    return

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
def report_chain( chain, filenames, f=sys.stdout):
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    fstr = "{:40s} "+(" {:8.4f}"*5)
    print("#variable                  2.5        25         50       75    97.5",file=f)
    for i,var in enumerate(["pH_slope","pH_off", "K1_f","K2_f"]):
        perc=np.percentile( chain[:,:,i].flatten(),[2.5,25,50,75,97.5])
        print(fstr.format(var,*perc),file=f)
    print("#DIC",file=f)
    fstr = "{:40s} "+(" {:8.1f}"*5)    
    for i,fname in enumerate(filenames):
        perc = np.percentile( chain[:,:,4+2*i].flatten(),[2.5,25,50,75,97.5])*1000
        print(fstr.format(fname[15:-4],*perc),file=f)
    print("#Alk",file=f)
    for i,fname in enumerate(filenames):
        perc = np.percentile( chain[:,:,4+2*i+1].flatten(),[2.5,25,50,75,97.5])*1000
        print(fstr.format(fname[15:-4],*perc),file=f)
    return

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def get_directory():
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -    
    from tkinter import Tk
    from tkinter import filedialog
    import os

    root = Tk()
    root.update()
    root.withdraw()

    current_directory = filedialog.askdirectory()
    root.update()
    return current_directory
