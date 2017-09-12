import numpy as np
import matplotlib.pyplot as plt
import matplotlib.transforms as transforms

from co2sys import co2sys
from co2sys import cc_py_mod as cc

consts = [10,8]
files  = ["co2sys_results_Luecker.txt", "co2sys_results_Millero_fresh.txt"]

for const,filen in zip( consts, files) :
    print(const,filen)

    co2_xls = np.genfromtxt( filen, comments="#" )
    # calculate with the new library
    nres,nsamp = np.shape(co2_xls)


    def format_vec( arr, fmt=""):
        return ((" {:"+fmt+"}")*len(arr)).format(*arr)


    co2_lib = np.zeros([nres,nsamp])
    co2_fort= np.zeros([nres,nsamp])
    for j,res in enumerate(co2_xls):
        T,S,DIC,TA, TP,TSi= list(map(float,res[[1,0,6,5, 3,4]]))
        #print("T=",T,"S=",S)
        pH,outp = co2sys.CC_solve( DIC*1e-6,TA*1e-6, T+273.15,S, TP=TP*1e-6, TSi=TSi*1e-6, const=const )
        co2_lib[j,[0,1, 7 ]]=[S,T, pH ]
        co2_lib[j,[6,5]] = [ DIC,TA ] 
        co2_lib[j,[  10,11,12,14,15,16]] = [              outp['HCO3'], outp['CO3'],outp['CO2'],outp['OH'],outp['AP'],outp['ASi']]
        co2_lib[j,5:7] = co2_lib[j,5:7]
        co2_lib[j,10:] = co2_lib[j,10:]*1e6

        outp = cc.py_interface_mod.cc_solve( [ DIC*1e-6,TA*1e-6, T+273.15,S,TP*1e-6,TSi*1e-6 ],const )
        co2_fort[j,:2] = [S,T]
        co2_fort[j,5:17] = outp[0:12]

        co2_fort[j,5:7] = co2_fort[j,5:7]*1e6
        co2_fort[j,9:] = co2_fort[j,9:]*1e6

        #print((" {:12.4g}"*14).format(*co2_fort[j,:14]))

    print( np.average(co2_lib-co2_xls,axis=0))

    # the x coords of this transformation are data, and the
    # y coord are axes

    plc = np.linspace(-0.025,+0.025,nres)
    val_index=[  7,   9,     10,    11,    12,   14,   15,  16]
    var_index=["pH","pCO2","HOC3","CO3","CO2","OH","Alk_P","Alk_Si"]
    cols=['red','darkgreen','blue']
    fig,axs = plt.subplots(3,4)

    for i,ax in enumerate( fig.axes):
        if( i < len(var_index) ):
            trans = transforms.blended_transform_factory(\
                        ax.transData, ax.transAxes)
            ax.set_title( var_index[i] )
            ax.text( 0.9, 0.90, "py-xls",   transform=trans, color=cols[0])
            ax.text( 1.0, 0.90, "fort-xls", transform=trans, color=cols[1])
            ax.text( 1.1, 0.90, "fort-py",  transform=trans, color=cols[2])
            k = val_index[i]
            ax.plot( 0.9+plc, co2_lib [:,k]-co2_xls[:,k], 'o',color=cols[0])
            ax.plot( 1.0+plc, co2_fort[:,k]-co2_xls[:,k], 'o',color=cols[1])
            ax.plot( 1.1+plc, co2_fort[:,k]-co2_lib[:,k], 'o',color=cols[2])
    plt.show()
    sdfjk
