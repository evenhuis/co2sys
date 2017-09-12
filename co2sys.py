import numpy as np
import numba

@numba.jit("f8(f8)",nopython=True)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def T_B( S ):
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Uppstrom, L., Deep-Sea Research 21:161-162, 1974:
    #return 4.157E-4 * S / 35.
    # Lee, Kim, Byrne, Millero, Feely, Yong-Ming Liu. 2010.  Geochimica Et Cosmochimica Acta 74 (6)
    return 0.0004326 * S / 35 

@numba.jit("f8(f8)",nopython=True)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def T_S( S ):
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Morris, A. W., and Riley, J. P., Deep-Sea Research 13:699-705, 1966:
    return (0.14 / 96.062) * (S / 1.80655 )

@numba.jit("f8(f8)",nopython=True)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def T_F( S) :
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    return 0.000067 * S / 18.9984 / 1.80655

@numba.jit("f8(f8,f8,i4)",nopython=True)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def K1_H2CO3( T, S, const=10 ):
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if( const== 8 ):
        K1_H2CO3 = np.exp(290.9097 - 14554.21 / T - 45.0575 * np.log(T))
    if( const==10 ):
        K1_H2CO3 = 3633.86/T - 61.2172 + 9.6777*np.log(T) - 0.011555*S + 0.0001152*S**2
        K1_H2CO3 = 10**(-K1_H2CO3)
    return K1_H2CO3

@numba.jit("f8(f8,f8,i4)",nopython=True)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def K2_H2CO3( T, S, const=10 ):
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if( const== 8 ): 
        K2_H2CO3 = np.exp(207.6548 - 11843.79 / T - 33.6485 * np.log(T))
    if( const==10 ):
        K2_H2CO3 = 471.8/T + 25.9290 - 3.16967*np.log(T) - 0.01781*S + 0.0001122*S**2
        K2_H2CO3 = 10**(-K2_H2CO3)
    return K2_H2CO3

@numba.jit("f8(f8,f8)",nopython=True)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def K_S( T, S ):
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ''' The equlibrium constant for Bisulphate ion
            K_SO4 [HSO4-] = [H+]_F [SO4--]
     
      Zeebe p259
      pH scale is free scale , pH_F
      This is used to convert 
          pH_F  <->   pH_T  <->     pH_SWS
         'free'     Hansson     seawater scale '''
    # Dickson, A. G., J. Chemical Thermodynamics, 22:113-127, 1990"
    I =  19.924*S          \
            / (1000. - 1.005*S )   

    logT = np.log(T)
    K_S =  -4276.1/T + 141.328 -   23.093*logT                \
          + (-13856./T + 324.57  -  47.986*logT ) * np.sqrt(I) \
          + ( 35474./T - 771.54  + 114.723*logT ) * I          \
                -2698./T * I**1.5 + 1776./T * I**2

    K_S = np.exp(K_S) * (1 - 0.001005 * S)

    return K_S

@numba.jit("f8(f8,f8)",nopython=True)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def K_F( T, S ):
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ''' equilibrium constant for 
        [HF] K_F = [H+] [F-] 
     
      from Zeebe p 260
      pH scale is Hansson, pH_T
    ''' 
    
    I  =     19.924*S      \
         /(1000. - 1.005*S)
    
    TS = T_S(S)
    KS = K_S(T,S)
   
    ## Dickson, A. G. and Riley, J. P., Marine Chemistry 7:89-99, 1979:
    #K_F = 1590.2/T - 12.641 + 1.525*np.sqrt(I)
    #                                   # this term converts 
    #                                   # 'free' -> Hansson
    #K_F = np.exp(K_F) * (1.0 - 0.001005*S)* (1.+ TS/KS)

    #  Perez,F.F. and Fraga, F., Marine Chemistry 21(2):161-168, 1987
    #S =10-40   T=9-33 oC
    lnKF = -874. / T - 0.111 * np.sqrt(S)+ 9.68 
    lnKF = -lnKF
    K_F = np.exp(lnKF)
    
    return K_F

@numba.jit("f8(f8,f8)",nopython=True)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def K1_P( T, S ):
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    K1_P = -4576.752/T +115.525 -18.453*np.log(T) \
          +(-106.736/T+0.69171)*np.sqrt(S) + (-0.65643/T-0.01844)*S
    K1_P = np.exp(K1_P)
    return K1_P

@numba.jit("f8(f8,f8)",nopython=True)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def K2_P( T, S ):
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    K2_P = -8814.715/T + 172.1033 -27.927*np.log(T) \
          +(-160.34/T+1.35661)*np.sqrt(S) + (0.37355/T-0.05778)*S
    K2_P = np.exp(K2_P)
    return K2_P

@numba.jit("f8(f8,f8)",nopython=True)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def K3_P( T, S ):
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    K3_P = -3070.75/T -18.126 + ( 17.27039/T + 2.81197 )*np.sqrt(S) \
            + (-44.994846/T - 0.09984)*S
    K3_P = np.exp(K3_P)
    return K3_P

@numba.jit("f8(f8,f8)",nopython=True)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def K_W( T, S ):
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    logT = np.log(T)
    K_W =  148.9802 - 13847.26/T  - 23.6521 * logT     \
         +(118.67/T  - 5.977 + 1.0495*logT)*np.sqrt(S) - 0.01615*S
    K_W = np.exp(K_W)
    return K_W

@numba.jit("f8(f8,f8)",nopython=True)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def K_Si( T, S ):
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    I  =     19.924*S      \
        /(1000. - 1.005*S)

#   K_Si = -8904.2/T + 117.385 -19.334*np.log(T)                        \
    K_Si = -8904.2/T + 117.4   -19.334*np.log(T)                        \
        + ( 3.5913  - 458.79 /T )*np.sqrt(I) + ( 188.74/T - 1.5998)*I   \
        + ( 0.07871 - 12.1652/T )*I**2

    K_Si = np.exp(K_Si)*(1.0-0.001005*S)
    return K_Si

@numba.jit("f8(f8,f8)",nopython=True)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def K_B( T, S ) :
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Dickson, A. G., Deep-Sea Research 37:755-766, 1990:
    sqrtS = np.sqrt(S)
    K_B = ( -8966.90 - 2890.53*sqrtS - 77.942*S + 1.728*S*sqrtS  \
           -0.0996*S**2 ) / T                                    \
           +148.0248 + 137.1942*sqrtS + 1.62142*S                \
           -(24.4344 +   25.085*sqrtS + 0.2474*S)*np.log(T)      \
          +0.053105*sqrtS*T
    #! Needed for converting from Hansson pH_T -> seawater pH_SWS
    TS = T_S( S ) ; KS = K_S(T,S)         # pH_F 
    TF = T_F( S ) ; KF = K_F(T,S)         # pH_T    
    K_B = np.exp( K_B ) * (1. + TS/KS)/( 1. + TS/KS + TF/KF )
    return K_B


@numba.jit("f8(  f8,  f8, f8,f8, i4,    optional(f8), optional(f8), optional(f8),optional(f8),optional(f8),f8,f8,optional(f8))",nopython=True)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def CC_solve_pH( DIC, TA, T,S, const=10, TP=0., TSi=0, TB=None, TS=None, TF=None, \
              K1_f=1.0, K2_f=1.0, pHi=None ):
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #! Needed for converting from Hansson pH_T -> seawater pH_SWS
    TS = T_S( S ) ; KS = K_S(T,S)         # pH_F 
    TF = T_F( S ) ; KF = K_F(T,S)         # pH_T

    if( TS is None ): TS = T_S(S)
    if( TF is None ): TF = T_F(S)

    SWS_2_T  = (1. + TS/KS)/( 1. + TS/KS + TF/KF )
    Free_2_T =  1. + TS/KS

    K1 = K1_H2CO3( T,S, const=const) *K1_f#/ SWS_2_T
    K2 = K2_H2CO3( T,S, const=const) *K2_f#/ SWS_2_T

    KW = K_W( T,S)                # pH_T
    KB = K_B( T,S )/SWS_2_T
    if( TB is None ): TB = T_B(S)

    K1P = K1_P(T,S)#/SWS_2_T
    K2P = K2_P(T,S)#/SWS_2_T
    K3P = K3_P(T,S)#/SWS_2_T
    KSi = K_Si(T,S)#/SWS_2_T

    #@numba.jit("f8(f8)",nopython=True)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def obj( pH ):
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        h = 10**(-pH)
        h_free = h/Free_2_T
        y =  DIC*(K1*h+2.*K1*K2)/(h*h+K1*h+K1*K2)   \
            - h_free + KW/h                         \
            - TA                                    \
           +TB  /(1.+h/KB)                          \
           +TP*(K1P*K2P*h+2*K1P*K2P*K3P-h**3)/(h**3+K1P*h**2+K1P*K2P*h+K1P*K2P*K3P) \
           -TF /(1.+KF/h_free)                       \
           +TSi/(1.+h/KSi)                          \
           -TS /(1.+KS/h_free)
        y = y*1e6
        return y

    #@numba.jit("f8(f8)",nopython=True)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def Dobj( pH ):
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        h = 10**(-pH)
        dy =   DIC*(K1 +2*K1*K2)/(h**2+K1*h+K1*K2) \
           -DIC*(K1*h+2*K1*K2)/(h**2+K1*h+K1*K2)**2*(2*h+K1) \
           -TB *1./(1+h/KB)**2 / KB                         \
           -KW/h**2 - 1./Free_2_T                         \
           +TP*(K1P*K2P             -3.*h**2)/(h**3+K1P*h**2+K1P*K2P*h+K1P*K2P*K3P) \
           -TP*(K1P*K2P*h+2*K1P*K2P*K3P-h**3)/(h**3+K1P*h**2+K1P*K2P*h+K1P*K2P*K3P)**2\
              *(3.*h**2+2.*K1P*h+K1P*K2P) \
           -TF  /(h+KF)**2 * KF    \
           +TSi /(KSi+h/KSi)**2      \
           -TS  /(h+KS*Free_2_T)**2*KS*Free_2_T    
        dy =dy*1e6  * (-np.log(10.)*10**(-pH))
        return dy

    if( pHi is None ):
       # find a rough range by bisection
       h0 = 12.0
       h1 = 7.0
       h2 = 3.0
       for i in range(8):
          h1 = (h0+h2)/2.0
          f0 = obj(h0)
          f1 = obj(h1)
          f2 = obj(h2)
          if(      ( f0<0. and f1>0.  and f2>0. ) or \
                   ( f0>0. and f1<0.  and f2<0. ) ):
              h2 = h1
          elif( ( f0<0. and f1<0.  and f2>0. ) or \
                ( f0>0. and f1>0.  and f2<0. ) ):
              h0 = h1
          else:
              break
       pH0=h1
    else:
       pH0=pHi

#   from scipy.optimize import newton
#   pH_opt=newton( obj,pH0,              tol=1e-6 )
#   pH_opt=newton( obj,pH0,fprime=Dobj,  tol=1e-6 )
#   pH0 = pH_opt
    for i in range(100):
       f0  =  obj( pH0 )
       df0 = Dobj( pH0 )
 
       pH1  = pH0 - f0/df0
 
       if( abs(f0)<1e-6 ): break
       pH0 = pH1

    return pH0

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def CC_solve( DIC, TA, T,S, const=10, TP=0., TSi=0, TB=None, TS=None, TF=None,  \
              K1_f=1.0, K2_f=1.0, pHi=None ):
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #! Needed for converting from Hansson pH_T -> seawater pH_SWS
    TS = T_S( S ) ; KS = K_S(T,S)         # pH_F 
    TF = T_F( S ) ; KF = K_F(T,S)         # pH_T

    if( TS==None ): TS = T_S(S)
    if( TF==None ): TF = T_F(S)

    SWS_2_T  = (1. + TS/KS)/( 1. + TS/KS + TF/KF )
    Free_2_T =  1. + TS/KS

    K1 = K1_H2CO3( T,S, const=const) *K1_f#/ SWS_2_T
    K2 = K2_H2CO3( T,S, const=const) *K2_f#/ SWS_2_T

    KW = K_W( T,S)                # pH_T
    KB = K_B( T,S )/SWS_2_T
    if( TB==None ): TB = T_B(S)

    K1P = K1_P(T,S)#/SWS_2_T
    K2P = K2_P(T,S)#/SWS_2_T
    K3P = K3_P(T,S)#/SWS_2_T
    KSi = K_Si(T,S)#/SWS_2_T

    #@numba.jit("f8(f8)",nopython=True)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def obj( pH ):
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        h = 10**(-pH)
        h_free = h/Free_2_T
        y =  DIC*(K1*h+2.*K1*K2)/(h*h+K1*h+K1*K2)   \
            - h_free + KW/h                         \
            - TA                                    \
           +TB  /(1.+h/KB)                          \
           +TP*(K1P*K2P*h+2*K1P*K2P*K3P-h**3)/(h**3+K1P*h**2+K1P*K2P*h+K1P*K2P*K3P) \
           -TF /(1.+KF/h_free)                       \
           +TSi/(1.+h/KSi)                          \
           -TS /(1.+KS/h_free)
        y = y*1e6
        return y

    #@numba.jit("f8(f8)",nopython=True)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    def Dobj( pH ):
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        h = 10**(-pH)
        dy =   DIC*(K1 +2*K1*K2)/(h**2+K1*h+K1*K2) \
           -DIC*(K1*h+2*K1*K2)/(h**2+K1*h+K1*K2)**2*(2*h+K1) \
           -TB *1./(1+h/KB)**2 / KB                         \
           -KW/h**2 - 1./Free_2_T                         \
           +TP*(K1P*K2P             -3.*h**2)/(h**3+K1P*h**2+K1P*K2P*h+K1P*K2P*K3P) \
           -TP*(K1P*K2P*h+2*K1P*K2P*K3P-h**3)/(h**3+K1P*h**2+K1P*K2P*h+K1P*K2P*K3P)**2\
              *(3.*h**2+2.*K1P*h+K1P*K2P) \
           -TF  /(h+KF)**2 * KF    \
           +TSi /(KSi+h/KSi)**2      \
           -TS  /(h+KS*Free_2_T)**2*KS*Free_2_T    
        dy =dy*1e6  * (-np.log(10.)*10**(-pH))
        return dy

    if( pHi is None ):
       # find a rough range by bisection
       h0 = 12.0
       h1 = 7.0
       h2 = 3.0
       for i in range(8):
          h1 = (h0+h2)/2.0
          f0 = obj(h0)
          f1 = obj(h1)
          f2 = obj(h2)
          if(      ( f0<0. and f1>0.  and f2>0. ) or \
                   ( f0>0. and f1<0.  and f2<0. ) ):
              h2 = h1
          elif( ( f0<0. and f1<0.  and f2>0. ) or \
                ( f0>0. and f1>0.  and f2<0. ) ):
              h0 = h1
          else:
              break
       pH0=h1
    else:
       pH0=pHi

#   from scipy.optimize import newton
#   pH_opt=newton( obj,pH0,              tol=1e-6 )
#   pH_opt=newton( obj,pH0,fprime=Dobj,  tol=1e-6 )
#   pH0 = pH_opt
    for i in range(100):
       f0  =  obj( pH0 )
       df0 = Dobj( pH0 )
 
       pH1  = pH0 - f0/df0
 
       if( abs(f0)<1e-6 ): break
       pH0 = pH1


    H = 10**(-pH0)
    H2 = H*H
    H3 = H*H*H

    outp={}
    denom = (H2+K1*H+K1*K2)
    outp[  'CO2'] = DIC*H2      /denom
    outp[ 'HCO3'] = DIC*H *K1   /denom
    outp[  'CO3'] = DIC   *K1*K2/denom
    outp[   'OH'] = KW/H

    denom = ( H3 + K1P*H2+ H*K1P*K2P+K1P*K2P*K3P )
    outp['H3PO4'] = TP*H3            /denom
    outp['H2PO4'] = TP*H2*K1P        /denom
    outp[ 'HPO4'] = TP*H *K1P*K2P    /denom
    outp[  'PO4'] = TP   *K1P*K2P*K3P/denom
    outp[   'AP'] = outp['HPO4'] +2*outp['PO4'] -outp['H3PO4']

    outp[  'ASi'] = TSi*KSi/(KSi+H)

    return pH0,outp
