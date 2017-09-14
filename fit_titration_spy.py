#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import glob

from co2sys import titration_util as tu
from emcee_tools import helper_emcee as he

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# get some files to process
filenames = glob.glob("PC_LIMS*.txt")

vol_data, aux_data = tu.get_multi_data( filenames)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# look for previous parameters
#
theta = tu.gen_theta_guess( vol_data, aux_data)
for i,filename in enumerate(filenames):
    # generate the guess
    theta_i = tu.get_theta_from_file( filename,"theta_save.txt")
    if( len(theta_i)!=0 ):
        theta[0:4] = theta_i[0:4]
        theta[4+2*i:4+2*(i+1)] = theta_i[4:6]
        theta[-1] = theta_i[-1]

# plot the fi
tu.plot_fit( vol_data, aux_data, theta )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# optimise
#
ch= input('Optimise? n')
while( ch not in ['n','N'] ):
    theta = he.optimise_ll( tu.log_prob, theta, vol_data, aux_data, live_plot=True, nloop=30 )
    plt.close('all')
    tu.plot_fit( vol_data, aux_data, theta )
    ch= input('Optimise? n')
plt.ioff()
plt.close()

# save the results
f = open('theta_new.txt', 'w')
tu.save_file( filenames,theta, f)
f.close()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# MCMC to do the stats
treat="pooled"
ch = input('sample chain? n')
while( not(ch in ["n","N"]) ):
    ch = input("Restart from chain? Y")
    if( ch in ["y","Y"] ):
        savefile="{}_chain.p".format(treat)
        if( os.path.isfile(savefile) ):
            probs,chain = he.load_chain(savefile)
            print("Chain loaded")
            p0 = chain[:,-1,:]
            samp = he.MCMC_restart( tu.log_prob, p0, vol_data,aux_data, nsteps=100, live_plot=True )
            fname = "{}_chain.p".format(treat)
            he.save_chain( fname, samp )
        else:
            print("chain file not found")
    else:
        samp = he.MCMC_all( tu.log_prob, theta, vol_data,aux_data , nsteps=30, nwalker=200 )
        fname = "{}_chain.p".format(treat)
        he.save_chain( fname, samp )
    ch = input("Keep sampling? Y")
plt.ioff()
plt.close()

ch = input('Load the chain? Y')
if( ch in ["y","Y"] ):
    fname = "{}_chain.p".format(treat)
    [prob,chain] = he.load_chain( fname )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# plot results and report
#
tu.report_chain( chain, filenames)
tu.plot_DIC_Alk( chain, filenames)




