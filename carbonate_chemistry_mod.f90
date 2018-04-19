!--------------------------------------------------------------------------------
module carbonate_chemistry_mod
!--------------------------------------------------------------------------------
implicit none

double precision, private, parameter :: T0K=-273.15d0

contains






! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
elemental subroutine CC_solve_DIC_Talk_old( DIC, Talk, T, S,              & ! INPUT
                                        CO2,  H,                      & ! OUTPUT
                                        TP, TSi,                      & ! OPTIONAL
                                        HCO3, CO3, OH,                & ! OPTIONAL
                                        omega_cal, omega_arg, pH,     & ! OPTIONAL
                                        AB, AP, ASi,                  & ! OPTIONAL
                                        const,                        & ! OPTIONAL
                                        mask )                          ! OPTIONAL
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! This solver takes more iterations than the one for CO2 and Talk
! Other choices for variables to solve for do not converge of physical solutions
!
! TODO
! Perhaps implementing fifth order equation for TA and DIC in Zeebe  pg 277
! and using Newton-Raphson to solve for the root would be faster.  Apparently
! this equation has only 1 real, positive root so it should be very stable.
!
double precision, intent(in)  :: DIC,         & ! Concentration CO2 aq        mol / kg-soln
                       Talk,        & ! Total Alkinity              mol / kg-soln
                       T,           & ! Temperature                 Kelvin
                       S              ! Salinity                    psu

double precision, intent(in),optional  ::     &
                       TP,          & ! Total Phosphate             mol / kg-soln
                       TSi            ! Total Silicate              mol / kg-soln
integer,intent(in),optional :: const  ! which set of constants to use
logical,intent(in),optional :: mask    ! if false, exit
                                          

double precision,intent(out)  :: H,           & ! Hydrogen Ion conc. Hansson  mol / kg-soln
                       CO2            ! Concentration CO2 aq        mol / kg-soln
double precision, intent(out),optional ::     &
                       HCO3,        & ! Bicarbonate ion conc        mol / kg-soln
                        CO3,        & ! Carbonate   ion conc        mol / kg-soln
                         OH,        & ! OH conc                     mol / kg-soln 
                       omega_cal,   & ! Calcite   saturation state  1
                       omega_arg,   & ! Aragonite saturation state  1
                       pH,          & ! pH                          Hansson scale
                       AB,          & ! Borate   Allalinity         mol / kg-soln
                       AP,          & ! Phospate Alkalinity         mol / kg-soln
                       ASi            ! Silicate Alkalinity         mol / kg-soln

double precision :: a, b, c,   & ! used in solution of quadratic equation
          CO3i,      & ! Carbonate concentration
          CO3i_prev, HCO3i

double precision :: K1,        & ! first  dissociation constant CO2 
          K2,        & ! second dissociation constant CO2
          Kw,        & ! dissociaton of H2O
          KB, TB,    & ! Boron    constant, total Boron
          KS, TS,    & ! Sulphate constant, total sulphate
          KF, TF,    & ! Fluroide constant, total fluoride
          K1P, K2P, K3P, APi, & ! Phosphate 
          KSi, ASii, & ! Silicate
          ABi,       & ! Boron
          H_F,       & ! The free H concentration
          H_T,       & ! The Hanson H concentration
          SWS_2_T      ! Convert _SWS to _T

integer :: n
integer,parameter :: nmax=200

if( present(mask) )then
   if( mask ) return
endif

! Needed for converting from Hansson pH_T -> seawater pH_SWS
TS = T_S( S ) ; KS = K_S(T,S)         ! pH_F 
TF = T_F( S ) ; KF = K_F(T,S)         ! pH_T
SWS_2_T = (1. + TS/KS)/( 1. + TS/KS + TF/KF )

! Get the equilibrium constants
K1 = K1_H2CO3(T,S)!*SWS_2_T          ! pH_SWS -> pH_T
K2 = K2_H2CO3(T,S)!*SWS_2_T          ! pH_SWS -> pH_T
     
KW = K_W( T,S )               ! pH_T
KB = K_B( T,S )                     ! pH_T

TB  = T_B(S)          

if( present(TP) )then   ! the CO2sys worksheet state these are in SWS unit
                        ! Zeebe says Total scale
                        ! the original source says nothing. Doesn't seems to make a difference
   K1P = K1_P(T,S)!*SWS_2_T          ! pH_T
   K2P = K2_P(T,S)!*SWS_2_T          ! pH_T
   K3P = K3_P(T,S)!*SWS_2_T          ! pH_T
else
   APi = 0.
endif

if( present(TSi))then
   KSi  = K_Si(T,S)                 ! pH_T
else
   ASii = 0.
endif


! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Iterate solution
! 0. Starting point
CO3i_prev = 0.
CO3i = Talk - DIC
n=0
do while ( abs((CO3i-CO3i_prev)/CO3i) > 0.00001 .and. n<nmax)
   ! 1. Solve for H
   a =  1.d0/(K1*K2)
   b =  1.d0/K2
   c = 1.d0-DIC/CO3i
   H = ( -b + sqrt(b**2 - 4*a*c ) )/(2.d0*a )

   H_F = H/( 1+TS/KS )  ! free H concentration

   ! 2. Solve for CO3

   ! See Zeebe p 39 for a disscussion of this odd definition of alkalinity
   if( present(TP )) APi  = TP*(-H/K1p      +   K2P/H + 2.*K2P*K3P/H**2)/ & 
                               ( H/K1p + 1. +   K2P/H +    K2P*K3P/H**2)
   if( present(TSi)) ASii = TSi/(1.+KSi/H)


   ! See Zeebe p 39 for a disscussion of this odd definition of alkalinity
   ABi = TB*KB/(H+KB)
   if( present(TSi)) ASii = TSi/(1.+H/KSi)
   if( present(TP )) APi  = TP*(-H**3            + K1P*K2P*H + 2.*K1P*K2P*K3P)/ &
                               ( H**3 + K1P*H**2 + K1P*K2P*H +    K1P*K2P*K3P)

   CO3i_prev = CO3i
   !               B(OH4)-   OH-  H_F^+    HSO4-             HF
   CO3i = ( Talk  -ABi    -KW/H  +H_F  +TS/(1.+KS/H_F)  +TF/(1.+KF/H)  &
                  -APi    -ASii              )*K2/(H+2.d0*K2)


  n=n+1
enddo

CO2  = CO3i*H**2/(K1*K2)

! optional outputs
if( present(pH)       ) pH        = -log10(H)
if( present(HCO3)     ) HCO3      = CO3i * H / K2
if( present( CO3)     ) CO3       = CO3i
if( present(omega_cal)) omega_cal = T_Ca(S) * CO3i / Ksp_cal(T,S)
if( present(omega_arg)) omega_arg = T_Ca(S) * CO3i / Ksp_arg(T,S)
if( present(AB       )) AB  = ABi
if( present(AP       )) AP  = APi
if( present(ASi      )) ASi = ASii
if( present(OH       )) OH  = KW/H 

!checks for other equations
HCO3i = CO3i * H / K2
   
return
end subroutine CC_solve_DIC_Talk_old

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
subroutine CC_solve_DIC_Talk( DIC, Talk, T, S,                            & ! INPUT
                                            CO2,  H,                      & ! OUTPUT
                                            TB, TP, TSi, TNH3, TS, TF,    & ! OPTIONAL INPUT
                                            HCO3, CO3, OH,                & ! OPTIONAL OUTPUT
                                            omega_cal, omega_arg, pH,     & ! OPTIONAL
                                            AB, AP, ASi,                  & ! OPTIONAL
                                            const,pCO2,                   & ! OPTIONAL
                                            mask, K1_f, K2_f, pHi )          ! OPTIONAL INPUT
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Attempt to solve the fifth order equation
!
! Perhaps implementing fifth order equation for TA and DIC in Zeebe  pg 277
! and using Newton-Raphson to solve for the root would be faster.  Apparently
! this equation has only 1 real, positive root so it should be very stable.
!
double precision, intent(in)  :: DIC,         & ! Concentration CO2 aq        mol / kg-soln
                       Talk,        & ! Total Alkinity              mol / kg-soln
                       T,           & ! Temperature                 Kelvin
                       S              ! Salinity                    psu

double precision, intent(in),optional  ::     &
                       TB,          & ! Total Borate                mol / kg-soln
                       TP,          & ! Total Phosphate             mol / kg-soln
                       TS,          & ! Total Sulphate              mol / kg-soln
                       TF,          & ! Total Fluoride              mol / kg-soln
                       TSi,         & ! Total Silicate              mol / kg-soln
                       TNH3,        & ! Total NH3                   mol / kg-sol
                       K1_f,        & ! factor for equlibrium constant 1, 
                       K2_f,        &
                       pHi             ! intial guess at pH
logical,intent(in),optional ::mask    ! if false, exit
integer,intent(in),optional :: const   ! which set of contsants to use
                                          

double precision,intent(out)  :: H,           & ! Hydrogen Ion conc. Hansson  mol / kg-soln
                       CO2            ! Concentration CO2 aq        mol / kg-soln
double precision, intent(out),optional ::     &
                       HCO3,        & ! Bicarbonate ion conc        mol / kg-soln
                        CO3,        & ! Carbonate   ion conc        mol / kg-soln
                         OH,        & ! OH conc                     mol / kg-soln 
                       omega_cal,   & ! Calcite   saturation state  1
                       omega_arg,   & ! Aragonite saturation state  1
                       pCO2,        & 
                       pH             ! pH                          Hansson scale

double precision, intent(out), optional :: &
                       AB,          & ! Borate   Allalinity         mol / kg-soln
                       AP,          & ! Phospate Alkalinity         mol / kg-soln
                       ASi            ! Silicate Alkalinity         mol / kg-soln

double precision :: a, b, c,   & ! used in solution of quadratic equation
          CO3i,      & ! Carbonate concentration
          CO3i_prev, HCO3i

double precision ::  &
          K1,              & ! first  dissociation constant CO2 
          K2,              & ! second dissociation constant CO2
          Kw,              & ! dissociaton of H2O
          KB, TB_i,  ABi,  & ! Boron    constant, total Boron
          KS, TS_i,        & ! Sulphate constant, total sulphate
          KF, TF_i,        & ! Fluroide constant, total fluoride
          K1P, K2P, K3P, TP_i,  APi,     & ! Phosphate 
          KSi,           TSi_i ,ASii,    & ! Silicate
          KNH3,          TNH3_i,ANH3i,  & ! Ammonia
          H_F,             & ! The free H concentration
          H_T,             & ! The Hanson H concentration
          SWS_2_T,         & ! Convert  _SWS to _T
          Free_2_T,        & ! Convert _Free to _T
          K1_fin,K2_fin 

integer :: n, i, &
           const_i
integer,parameter :: nmax=200

! variables for the secant method root finder
double precision :: h0, h1, h2, &
                    f0, f1, f2, df0
                     
const_i = 10      ! default to Lueker
if( present(const) ) const_i=const

if( present(mask) )then
   if( mask ) return
endif


! Needed for converting from Hanssion pH_T -> seawater pH_SWS
TS_i = T_S( S ) ; KS = K_S(T,S)         ! pH_F 
TF_i = T_F( S ) ; KF = K_F(T,S)         ! pH_T
if( present(TS) ) TS_i = TS
if( present(TF) ) TF_i = TF

 SWS_2_T = (1. + TS_i/KS)/( 1. + TS_i/KS + TF_i/KF )
Free_2_T =  1. + TS_i/KS

! Get the equilibrium constants
K1_fin =1.d0
if( present(K1_f) ) K1_fin = K1_f
K2_fin =1.d0
if( present(K2_f) ) K2_fin = K2_f

K1 = K1_H2CO3( T,S, const=const_i)*K1_fin!/SWS_2_T             ! pH_SWS -> pH_T
K2 = K2_H2CO3( T,S, const=const_i)*K2_fin!/SWS_2_T             ! pH_SWS -> pH_T
     
KW = K_W( T,S, const=const_i)                  ! pH_T
KB = K_B( T,S )/SWS_2_T   
TB_i = T_B(S)
select case(const_i)
case( 8 ) ; TB_i  = 0
end select
if( present(TB) ) TB_i = TB

APi  = 0.
TP_i =0.d0
if( present(TP) ) TP_i = TP
                        ! Zeebe says Total scale
                        ! the original source says nothing. Doesn't seems to make a difference
K1P = K1_P(T,S)!*SWS_2_T          ! pH_T
K2P = K2_P(T,S)!*SWS_2_T          ! pH_T
K3P = K3_P(T,S)!*SWS_2_T          ! pH_T

KSi  = K_Si(T,S)!/SWS_2_T         ! pH_T
TSi_i = 0.
if( present(TSi)) TSi_i = TSi

!write(19,*)"K1,K2",K1,K2
!write(19,*)"KW",KW
!write(19,*)"TB,KB",TB_i,KB
!write(19,*)"TS,KS",TS_i,KS
!write(19,*)"TF,KF",TF_i,KF
!stop 1

if(  present(pHi) )then
   h1 = pHi
else
   ! find a rough range by bisection
   h0 = 14.0  
   h1 = 9.0   
   h2 = 1.0   
   do i = 1, 8
      h1 = (h0+h2)/2.d0
      f0 = obj( 10**(-h0) )
      f1 = obj( 10**(-h1) )
      f2 = obj( 10**(-h2) ) 
      !write(19,*) " - - -"
      !write(19,*)  h0, f0
      !write(19,*)  h1, f1
      !write(19,*)  h2, f2

      if(      ( f0<0. .and. f1>0.  .and. f2>0. ) .or. &
               ( f0>0. .and. f1<0.  .and. f2<0. ) )then
          h2 = h1 
       elseif( ( f0<0. .and. f1<0.  .and. f2>0. ) .or. &
               ( f0>0. .and. f1>0.  .and. f2<0. ) )then
          h0 = h1
       else 
          exit
       endif
   enddo
endif

!write(19,*),"#Newton"
!write(19,*) "TS",T,S
!write(19,*) K1,K2
!write(19,*) free_2_T
!write(19,*) KW
!write(19,*) "H",10**(-h1)
!write(19,*) "obj",h1,obj(10**(-h1)),Dobj(10**(-h1))
h0 = h1
do i = 1, 100
    f0 =  obj( 10**(-h0) ) 
   df0 = Dobj( 10**(-h0) ) * (-log(10.)*10**(-h0))

   h1  = h0 - f0/df0

   !write(19,*) ,i,h0,f0
   if( abs(f0)<1e-4 ) exit
   h0 = h1 
enddo


! fill in the solution now
H   = 10**(-h1)
CO2 = DIC/( 1. + K1/H + K1*K2/H**2 )
CO3i = DIC/(1.+H/K2 + H**2/K1/K2 )
APi  = +TP_i*(K1P*K2P*h+2*K1P*K2P*K3P-h**3)/(h**3+K1P*h**2+K1P*K2P*h+K1P*K2P*K3P)
ASii = TSi_i/(1+h/KSi)

! optional outputs
if( present(pH)       ) pH        = h1
if( present(HCO3)     ) HCO3      = DIC/(1.+H/K1 + K2/H)
if( present( CO3)     ) CO3       = CO3i
if( present(omega_cal)) omega_cal = T_Ca(S) * CO3i / Ksp_cal(T,S)
if( present(omega_arg)) omega_arg = T_Ca(S) * CO3i / Ksp_arg(T,S)
if( present(OH       )) OH  = KW/H 
if( present(AB       )) AB  = TB_i/(1.+H/KB)
if( present(AP       )) AP  = APi
if( present(ASi      )) ASi = ASii
if( present(pCO2     ))then
   pCO2 = CO2/K0_CO2(T,S) / CO2_fugacity_const(T)
endif


!checks for other equations
HCO3i = CO3i * H / K2
   
return

contains

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function obj(h) result( y )
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! this is the objective function that we need to find the roots of
double precision, intent(in) :: h
double precision             :: y

double precision :: h_free


! This is from page 277 of Zeebe
!y =   DIC*( K1*h**2 +2.*K1*K2*h )*(KB+h)  &
!    -( (KB+h)*(Talk*h-Kw+h**2) - KB*TBi*h )*(h**2+K1*h+K1*K2)  

! old
!h_free=h/Free_2_T
!
! y =  DIC*(K1*H+2.*K1*K2)/(H*H+K1*h+K1*K2)   &
!     - h_free + Kw/h                         &
!     - Talk                                  &
!    +TB_i  /(1.+h/KB)                        &
!    +TP_i*(K1P*K2P*h+2*K1P*K2P*K3P-h**3)/(h**3+K1P*h**2+K1P*K2P*h+K1P*K2P*K3P)  &
!    -TF_i /(1.d0+KF/h_free)                        &
!    +TSi_i/(1.d0+h/Ksi)                       &
!    -TS_i /(1.d0+KS*Free_2_T/h_free)              
!
!y = y*1d6

h_free = h/Free_2_T
y =  DIC*(K1*h+2.*K1*K2)/(h*h+K1*h+K1*K2)   &
    - h_free + KW/h                         &
    - Talk                                  &
   +TB_i  /(1.+h/KB)                          &
   +TP_i*(K1P*K2P*h+2*K1P*K2P*K3P-h**3)/(h**3+K1P*h**2+K1P*K2P*h+K1P*K2P*K3P) &
   -TF_i /(1.+KF/h_free)                      &
   +TSi_i/(1.+h/KSi)                          &
   -TS_i /(1.+KS/h_free)

y = y*1d6

return
end function obj

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
pure function Dobj(h) result( dy )
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! this is the objective function that we need to find the roots of
double precision, intent(in) ::  h
double precision             :: dy

double precision :: h_free

! This is from page 277 of Zeebe
h_free = h/Free_2_T
!dy =    DIC*( 2.*K1*h    +2.*K1*K2   )*(KB+h) &
!       +DIC*(    K1*h**2 +2.*K1*K2*h )        &
!       -( (KB+h)*(Talk  +2.*h   ) &
!               + (Talk*h-Kw+h**2) - KB*TB_i   )*(   h**2+K1*h+K1*K2)  &
!       -( (KB+h)*(Talk*h-Kw+h**2) - KB*TB_i*h )*(2.*h   +K1        )
!

 dy =   DIC*(K1 +2*K1*K2)/(H**2+K1*h+K1*K2) &
       -DIC*(K1*h+2*K1*K2)/(H**2+K1*h+K1*K2)**2*(2*h+K1) &
       -TB_i *1./(1+h/KB)**2 / KB                         &
       -Kw/h**2 - 1./Free_2_T                         &
       +TP_i*(K1P*K2P             -3.*h**2)/(h**3+K1P*h**2+K1P*K2P*h+K1P*K2P*K3P) &
       -TP_i*(K1P*K2P*h+2*K1P*K2P*K3P-h**3)/(h**3+K1P*h**2+K1P*K2P*h+K1P*K2P*K3P)**2 &
        *(3.*h**2+2.*K1P*h+K1P*K2P)    &
       -TF_i  /(1.+KF/h_free)**2 * KF*Free_2_T    &
       +TSi_i /(1.+h /KSi   )**2 / KSi   &
       -TS_i  /(1.+KS/h_free)**2*KS*Free_2_T
         

dy = dy*1d6

return
end function Dobj

end subroutine CC_solve_DIC_Talk

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
elemental subroutine CC_solve_DIC_Talk_CO2sys( DIC, Talk, T, S,              & ! INPUT
                                        CO2,  H,                      & ! OUTPUT
                                        TP, TSi,                      & ! OPTIONAL
                                        HCO3, CO3, OH,                & ! OPTIONAL
                                        omega_cal, omega_arg, pH,     & ! OPTIONAL
                                        AB, AP, ASi,                  & ! OPTIONAL
                                        const,K1_f,                   & ! OPTIONAL
                                        mask )                          ! OPTIONAL
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! a reimplementation of the CO2sys version of this subroutine
double precision, intent(in)  :: DIC,         & ! Concentration CO2 aq        mol / kg-soln
                       Talk,        & ! Total Alkinity              mol / kg-soln
                       T,           & ! Temperature                 Kelvin
                       S              ! Salinity                    psu

double precision, intent(in),optional  ::     &
                       TP,          & ! Total Phosphate             mol / kg-soln
                       TSi,         & ! Total Silicate              mol / kg-soln
                       K1_f           ! factor for equlibrium constant 1
logical,intent(in),optional ::mask    ! if false, exit
integer,intent(in),optional :: const   ! which set of contsants to use


double precision,intent(out)  :: H,           & ! Hydrogen Ion conc. Hansson  mol / kg-soln
                       CO2            ! Concentration CO2 aq        mol / kg-soln
double precision, intent(out),optional ::     &
                       HCO3,        & ! Bicarbonate ion conc        mol / kg-soln
                        CO3,        & ! Carbonate   ion conc        mol / kg-soln
                         OH,        & ! OH conc                     mol / kg-soln 
                       omega_cal,   & ! Calcite   saturation state  1
                       omega_arg,   & ! Aragonite saturation state  1
                       pH             ! pH                          Hansson scale

double precision, intent(out), optional :: &
                       AB,          & ! Borate   Allalinity         mol / kg-soln
                       AP,          & ! Phospate Alkalinity         mol / kg-soln
                       ASi            ! Silicate Alkalinity         mol / kg-soln

double precision :: a, b, c,   & ! used in solution of quadratic equation
          CO3i,      & ! Carbonate concentration
          CO3i_prev, HCO3i, OHi

double precision ::  &
          K1,        & ! first  dissociation constant CO2 
          K2,        & ! second dissociation constant CO2
          ACi,       & ! Carbon alkalitnithy
          Kw,            & ! dissociaton of H2O
          KB, TB, ABi,   & ! Boron    constant, total Boron
          K1P, K2P, K3P, TPi, APi,  & ! Phosphate 
          KSi,           TSii,ASii, & ! Silicate
          KS, TS,    & ! Sulphate constant, total sulphate
          KF, TF,    & ! Fluroide constant, total fluoride
          H_F,       & ! The free H concentration
          H_T,       & ! The Hanson H concentration
          SWS_2_T,   & ! Convert _SWS to _T
          K1_fin

double precision :: deltapH, pHTol, Residual, Denom, freetotot, ln10, pHGuess, HF, HFree, &
                     HSO4, PhosBot, PhosTop, Slope

integer :: n, i
integer,parameter :: nmax=200

! variables for the secant method root finder
freetotot=1.d0

! Needed for converting from Hansson pH_T -> seawater pH_SWS
TS = T_S( S ) ; KS = K_S(T,S)         ! pH_F 
TF = T_F( S ) ; KF = K_F(T,S)         ! pH_T
SWS_2_T = (1. + TS/KS)/( 1. + TS/KS + TF/KF )
freetotot=SWS_2_T

! Get the equilibrium constants
K1 = K1_H2CO3(T,S)!*SWS_2_T          ! pH_SWS -> pH_T
K2 = K2_H2CO3(T,S)!*SWS_2_T          ! pH_SWS -> pH_T

KW = K_W( T,S )               ! pH_T

TB  = T_B(S)
KB = K_B( T,S )               ! pH_T

K1P = K1_P(T,S)!*SWS_2_T          ! pH_T
K2P = K2_P(T,S)!*SWS_2_T          ! pH_T
K3P = K3_P(T,S)!*SWS_2_T          ! pH_T
TPi=0.d0
if( present(TP) ) TPi = TP

TSii=0.d0



pHGuess = 8. !   : ' this is the first guess
pHTol = 0.0001 !: ' this is .0001 pH units
ln10 = Log(10.d0)
pH = pHGuess
do
   H = 10.d0**(-pH)
   Denom = (H**2 + K1*H + K1*K2)
   ACi = DIC * K1 * (H + 2.*K2) / Denom
   ABi = TB * KB / (KB + H)
   OHi = KW / H

   PhosTop = K1P*K2P*H + 2*K1P*K2P*K3P - H**3
   PhosBot = H**3 + K1P*H**2 + K1P*K2P*H + K1P*K2P*K3P
   APi  = TPi * PhosTop / PhosBot

   ASii = TSii * KSi / (KSi + H)
   FREEtoTOT = (1 + TS / KS)        !: ' pH scale conversion factor
   Hfree = H / FREEtoTOT            ! ' for H on the total scale
   HSO4 = TS / (1 + KS / Hfree)     ! ' since KS is on the free scale
   HF   = TF / (1 + KF / Hfree)     ! ' since KF is on the free scale
   Residual = Talk - ACi - ABi  - OHi - APi  - ASii  + Hfree + HSO4 + HF

   !  find Slope dTA/dpH:
   !  (this is not exact, but keeps all important terms):
   Slope = ln10 * (DIC*K1*H * (H**2 + K1 * K2 + 4*H*K2) / Denom**2 + ABi * H / (KB + H) + OHi + H)
   deltapH = Residual / Slope      ! ' this is Newton's method
   if( abs(deltapH)>1 ) deltapH=sign(deltapH,0.5d0)
   ! to keep the jump from being too big:

    pH = pH + deltapH
   if( deltapH<pHtol ) exit
enddo

HCO3i = DIC / ( H/K1+1.d0+K2/H ) 
 CO3i = K2*HCO3/H
CO2   = CO3i*H**2/(K1*K2)
! optional outputs
if( present(pH)       ) pH        = -log10(H)
if( present(HCO3)     ) HCO3      = HCO3i
if( present( CO3)     ) CO3       =  CO3i
if( present(omega_cal)) omega_cal = T_Ca(S) * CO3i / Ksp_cal(T,S)
if( present(omega_arg)) omega_arg = T_Ca(S) * CO3i / Ksp_arg(T,S)
if( present(AB       )) AB  = ABi
if( present(AP       )) AP  = APi
if( present(ASi      )) ASi = ASii
if( present(OH       )) OH  = KW/H


return
end subroutine CC_solve_DIC_Talk_CO2sys



! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
subroutine CC_solve_pH_Alk( pH, Talk, T, S, &      ! INPUT
                            DIC, CO2 )
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
double precision, intent(in)  :: pH,  Talk, T, S
double precision, intent(out) :: DIC, CO2

double precision :: H, KW, KB, TB, K1,K2

! Get the equilibrium constants
K1 = K1_H2CO3(T,S)!*SWS_2_T          ! pH_SWS -> pH_T
K2 = K2_H2CO3(T,S)!*SWS_2_T          ! pH_SWS -> pH_T

KW = K_W     (T,S)                  ! pH_T
KB = K_B( T,S )                     ! pH_T

TB  = T_B(S)

H    = 10.**(-pH)
CO2  = ( Talk - KB*TB/(KB+h) + Kw/h + h )/(K1/h+2*K1*K2/h**2 )
DIC =  CO2*(1+K1/H+K1*K2/H**2)

return
end subroutine CC_solve_pH_Alk


! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
subroutine CC_solve_CO2_pH( CO2, pH,   T, S,              & ! INPUT
                                        DIC, Talk,                    & ! OUTPUT
                                        TP, TSi,                      & ! OPTIONAL
                                        HCO3, CO3,                    & ! OPTIONAL
                                        omega_cal, omega_arg,         & ! OPTIONAL
                                        AB, AP, ASi,                  & ! OPTIONAL
                                        mask )                          ! OPTIONAL
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                                      ! Variable                    Units 
double precision, intent(in) :: CO2,& ! Concentration CO2 aq        mol / kg-soln
                       pH,          & ! pH                          Hansson scale
                       T,           & ! Temperature                 Kelvin
                       S              ! Salinity                    psu

double precision, intent(in),optional  ::     &
                       TP,          & ! Total Phosphate             mol / kg-soln
                       TSi            ! Total Silicate              mol / kg-soln

double precision, intent(out) :: DIC,         & ! Carbon Alkalinity mol / kg-soln
                                 Talk           ! Total Alkinity    mol / kg-soln
double precision, intent(out),optional ::     &
                       HCO3,        & ! Bicarbonate ion conc        mol / kg-soln
                        CO3,        & ! Carbonate   ion conc        mol / kg-soln
                       omega_cal,   & ! Calcite   saturation state  1
                       omega_arg,   & ! Aragonite saturation state  1
                       AB,          & ! Borate   Alkalinity         mol / kg-soln 
                       AP,          & ! Phospate Alkalinity         mol / kg-soln
                       ASi            ! Silicate Alkalinity         mol / kg-soln

logical,intent(in),optional ::mask    ! if false, exit

double precision :: a, b, c,   & ! used in solution of quadratic equation
          Calk,      & ! Carbon alkalinity
          Calk_prev, H_in    !     "                     "

double precision :: K1,           & ! first  dissociation constant CO2 
          K2,           & ! second dissociation constant CO2
          Kw,           & ! dissociaton of H2O
          KB, TB, ABi,  & ! Boron    constant, total Boron
          KS, TS,       & ! Sulphate constant, total sulphate
          KF, TF,       & ! Fluroide constant, total fluoride
          CO3i,         & ! Carbonate ion conc for internal calculations
          K1P, K2P, K3P, APi, & ! Phosphate alkalinity
          KSi, ASii,          & ! Silicate  alkalinity
          H_F,                & ! The free H concentration
                       H,           & ! Hydrogen Ion conc. Hansson  mol / kg-soln
          SWS_2_T               ! Convert _SWS to _T

if( present(mask) )then
   if( mask ) return
endif

! Needed for converting from Hanssion pH_T -> seawater pH_SWS
TS = T_S( S ) ; KS = K_S(T,S)         ! pH_F
TF = T_F( S ) ; KF = K_F(T,S)         ! pH_T
SWS_2_T = (1. + TS/KS) / ( 1. + TS/KS + TF/KF )

! Get the equilibrium constants
K1 = K1_H2CO3(T,S)                  ! pH_T
K2 = K2_H2CO3(T,S)                  ! pH_T
KW = K_W     (T,S)                  ! pH_T

TB = T_B(S)
KB = K_B( T,S )                     ! pH_T


if( present(TP) )then
   K1P = K1_P(T,S)
   K2P = K2_P(T,S)
   K3P = K3_P(T,S)
else
   APi = 0.
endif

if( present(TSi))then
   KSi  = K_Si(T,S)
else
   ASii = 0.
endif


H = 10**(-pH)
DIC   = CO2*(1.d0+K1/H + K1*K2/H**2)
Talk  = CO2*(K1/H + 2.d0*K1*K2/H**2) + TB*KB/(KB+H) + KW/H - H

!if( present(HCO3) ) HCO3 = DIC/( 1.d0 + H/K1 + K2/H )
!if( present( CO3) )  CO3 = DIC/( 1.d0 + H/K2 + H**2/(K1*K2) )

return
end subroutine CC_solve_CO2_pH


! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
subroutine CC_solve_CO2_Talk(           CO2, Talk, T, S,              & ! INPUT
                                        DIC, H,                       & ! OUTPUT
                                        TP, TSi,                      & ! OPTIONAL
                                        HCO3, CO3,                    & ! OPTIONAL
                                        omega_cal, omega_arg, pH,     & ! OPTIONAL
                                        AB, AP, ASi,                  & ! OPTIONAL
                                        mask )                          ! OPTIONAL
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                                      ! Variable                    Units 
double precision, intent(in)  :: CO2,         & ! Concentration CO2 aq        mol / kg-soln
                       Talk,        & ! Total Alkinity              mol / kg-soln
                       T,           & ! Temperature                 Kelvin
                       S              ! Salinity                    psu

double precision, intent(in),optional  ::     &
                       TP,          & ! Total Phosphate             mol / kg-soln
                       TSi            ! Total Silicate              mol / kg-soln

double precision, intent(out) :: DIC,         & ! Carbon Alkalinity           mol / kg-soln
                       H              ! Hydrogen Ion conc. Hansson  mol / kg-soln
double precision, intent(out),optional ::     & 
                       HCO3,        & ! Bicarbonate ion conc        mol / kg-soln
                        CO3,        & ! Carbonate   ion conc        mol / kg-soln
                       omega_cal,   & ! Calcite   saturation state  1
                       omega_arg,   & ! Aragonite saturation state  1
                       pH,          & ! pH                          Hansson scale
                       AB,          & ! Borate   Alkalinity         mol / kg-soln 
                       AP,          & ! Phospate Alkalinity         mol / kg-soln
                       ASi            ! Silicate Alkalinity         mol / kg-soln

logical,intent(in),optional ::mask    ! if false, exit

double precision :: a, b, c,   & ! used in solution of quadratic equation
          Calk,      & ! Carbon alkalinity
          Calk_prev    !     "                     "

double precision :: K1,           & ! first  dissociation constant CO2 
          K2,           & ! second dissociation constant CO2
          Kw,           & ! dissociaton of H2O
          KB, TB, ABi,  & ! Boron    constant, total Boron
          KS, TS,       & ! Sulphate constant, total sulphate
          KF, TF,       & ! Fluroide constant, total fluoride
          CO3i,         & ! Carbonate ion conc for internal calculations
          K1P, K2P, K3P, APi, & ! Phosphate alkalinity
          KSi, ASii,          & ! Silicate  alkalinity
          H_F,                & ! The free H concentration
          SWS_2_T               ! Convert _SWS to _T

if( present(mask) )then
   if( mask ) return
endif

! Needed for converting from Hanssion pH_T -> seawater pH_SWS
TS = T_S( S ) ; KS = K_S(T,S)         ! pH_F
TF = T_F( S ) ; KF = K_F(T,S)         ! pH_T
SWS_2_T = (1. + TS/KS) / ( 1. + TS/KS + TF/KF )

! Get the equilibrium constants
K1 = K1_H2CO3(T,S)                  ! pH_T
K2 = K2_H2CO3(T,S)                  ! pH_T
KW = K_W     (T,S)                  ! pH_T

TB = T_B(S)
KB = K_B( T,S )                     ! pH_T
         
if( present(TP) )then
   K1P = K1_P(T,S)
   K2P = K2_P(T,S)
   K3P = K3_P(T,S)
else
   APi = 0.
endif

if( present(TSi))then
   KSi  = K_Si(T,S)
else
   ASii = 0.
endif

! 0. Approximate carbon alkalinity
Calk_prev = 0.d0
Calk      = Talk

! iterate solution
do while ( abs((Calk-Calk_prev)/Calk) > 1.d-7)


! Now solve the quadratic equation for H
!  Calk * H^2  -   [ CO2 K1  ]* H - [2 K1 K2 CO2 ] = 0
      a   =  Calk
      b   =  -    K1   *CO2
      c   = -2.d0*K1*K2*CO2  != b * 2*K2
      H   =  ( -b + sqrt(b**2 - 4*a*c ) )/(2.d0*a )   ! take the  strictly positve root
      H_F = H/( 1+TS/KS )                             ! free H concentration
!     H   =  ( -b/(2*a) + sqrt( [ b/2a ]**2 - c/a )    
!   let XX = k1*CO2/(2.d0*Calk) = b/2a
!     XX   = K1*CO2/(2.d0*Calk)
!     H   =        XX   + sqrt( XX**2  - c/a )
!                  XX   + sqrt( XX**2 + 2*b*K2/a )
!                  XX   + sqrt( XX**2 + 4*(b/2a)*K2
!     H   =        XX   + sqrt( XX**2 + 4* XX * K2 )

! Now update the Carbon alkalinity
!    TA   = CA    + [B(OH)4-] + [OH-] - [H+]
!=>  CA   = TA    - TB*KB/(H+KB)  -KW/H + H
     Calk_prev = Calk
    !Calk = Talk  - TB*KB/(H+KB)  -KW/H+H

   ! See Zeebe p 39 for a disscussion of this odd definition of alkalinity
   ABi = TB*KB/(H+KB)
   if( present(TSi)) ASii = TSi/(1.+H/KSi)
   if( present(TP )) APi  = TP*(-H**3            + K1P*K2P*H + 2.*K1P*K2P*K3P)/ &
                               ( H**3 + K1P*H**2 + K1P*K2P*H +    K1P*K2P*K3P)
   
   !               B(OH4)-   OH-  H_F^+    HSO4-             HF
   Calk = ( Talk  -ABi    -KW/H  +H_F  +TS/(1.+KS/H_F)  +TF/(1.+KF/H)  &
                  -APi    -ASii              )
end do

CO3i = K1 * K2 * CO2 / H**2
DIC  = CO2 + K1*CO2/H + CO3i

! optional outputs
if( present(pH)       ) pH        = -log10(H)
if( present(HCO3)     ) HCO3      = K1 * CO2 / H
if( present( CO3)     ) CO3       = CO3i
if( present(omega_cal)) omega_cal = T_Ca(S) * CO3i / Ksp_cal(T,S)
if( present(omega_arg)) omega_arg = T_Ca(S) * CO3i / Ksp_arg(T,S)
if( present(AB       )) AB  = ABi
if( present(AP       )) AP  = APi
if( present(ASi      )) ASi = ASii

return
end subroutine CC_solve_CO2_Talk

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
elemental function rho_pw( T ) result( p )
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! The density of pure water as a function of T for Salinty=0, Press = 1atm
! From Zeebe p269
                          ! variable      units
double precision,intent(in) :: T    ! Temperature   Kelvin
double precision            :: p    ! Density       kg / m**3

double precision            :: Tc   ! Temperature   Celcius

Tc = T - 273.15

p = 999.842594 +6.793952d-2*Tc     -9.095290d-3*Tc**2 &
               +1.001685d-4*Tc**3  -1.120083d-6*Tc**4 &
               +6.536332d-9*Tc**5

return
end function rho_pw

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
elemental function rho_sw( T, S ) result( p )
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! The density kg m^-3 of sea water as a function of temperature and salinity
! From Zeebe p270
! NB. the dependence on T,S is reversed from Zeebe so it is consistent with
!     the ordering of T,S in the equlibrium constants K_*( T,S )
                          ! variable      units
double precision,intent(in) :: T, & ! Temperature   Kelvin
                     S    ! Salinity      psu
double precision            :: p    ! Density       kg / m**3

double precision            :: Tc   ! Temperature   Celcius
double precision            :: A, B, C

Tc = T - 273.15

A  =  8.24493d-1 -4.0899d-3*Tc    +7.6438d-5*Tc**2 &
                 -8.2467d-7*Tc**3 +5.3875d-9*Tc**4

B  = -5.72466d-3 +1.0227d-4*Tc    -1.6546d-6*Tc**2

C  =  4.83140d-4

p = rho_pw(T) + A*S + B*S**1.5 + C*S**2

return
end function rho_sw


! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
elemental function CO2_fugacity_const( T )  result( c )
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! calculates the fugacity of CO2 from the partial pressure
!
!     fCO2 = fugacity_const * pCOS
!
! from Zeebe p 65
double precision,intent(in) :: T
double precision            :: c

double precision ::           B, delta   
double precision,parameter :: p=101325.,&  ! pressure in Pa, 1atm = 101325Pa
                    R=8.314      ! Gas constant

B = (-1636.75 + 12.0408*T - 3.27957d-2*T**2 &
     +3.16528d-5*T**3 )   * 1.d-6
delta = ( 57.7 - 0.118*T) * 1.d-6
c = exp( p*(B+2.*delta)/(R*T ) )

return
end function CO2_fugacity_const

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
elemental function Talk_from_Sal(S) result( Talk )
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Not sure of the orgin of this relationship
! Appears to be the fit to Fig 1.2.16 of Zeebe p 50
                           ! Variable          Units
double precision,intent(in) :: S     ! Salinity          psu
double precision            :: Talk  ! Total Alkalinity  mol / (kg-soln)

Talk = (2300.+(66.3*(S-35.0))) * 1.d-6
return
end function Talk_from_Sal

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function mu_sw( T, S ) result( mu )
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! The dynamical viscosity of sea water
!  From equations  22 and 23 of 
!  Sharqawy, M. (2010). o
!  Thermophysical properties of seawater: A review of existing correlations and data. 
!  Desalination and Water, 16(10), 354–380. 
!
                           ! Variable          Units
double precision,intent(in) :: T, &  ! Temperaturature   degrees K
                     S     ! Salinity          psu            
double precision            :: mu    ! dyn. viscosity    kg m^-1 s^-1

double precision :: mu_w, T_C, A, B, Sk

T_C = T - 273.15
Sk  = S/1000.

! pure water calculation
mu_w = 4.2844E-5 + 1.d0/( 0.157*(T_C+64.993)**2 - 91.296 )

A = 1.541 + 1.998E-2*T_C - 9.520E-5*T_C**2
B = 7.974 - 7.561E-2*T_C + 4.724E-4*T_C**2

mu = mu_w * ( 1. + A*Sk + B*Sk**2 )

return
end function mu_sw

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
double precision elemental function K0_O2( T, S )
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Equlibrium constant for O2 solubilty in water, Henry's Law
!
! THIS IS FOR STANDARD ATMOSPHERIC COMPOSITION - NEAC
! and 100% WATER VAPOUR CONTENT C^{sa}
!
!  sol       gas
!
! [O2] = K0 fCO2
!
! fO2 : fugacity, or approx. the partial pressure of O2A
! From Weiss 1970
! 
! the units are:
!   O2 :  mol / kg_soln
!  fO2 :  atm
!
! From: Battino, Rettich, Tominaga
!       The Solubility of oxygen and Ozone in Liquids
!       J. Phys. Chem. Ref. Data Vol 12, No. 2 1983
! Eq 19 with Table 4 column 3
 double precision, intent(in) :: T,S
 
 double precision,parameter  :: &
      A1= -1282.8704, A2= 36619.96,    &
      A3=   223.1396, A4=    -0.354707, &
      B1= 5.957E-3,   B2=-3.7353,     B3= 3.68E-6
 
                     
 K0_O2 =  A1 + A2/T +A3*log(T) + A4*T + S*(B1+B2/T) + B3*S**2 
 K0_O2 = exp( K0_O2 )
 
 ! convert from NEAC, 20.94% oxygen to 1atm and from uM to M
 K0_O2 = K0_O2 / 0.2094  * 1E-6


!! from
!!  Oxygen solubility in sea water
!!  limnol. Oceanography 37(6) 1992 1307
!!  Garcia and Gordon
!double precision,parameter :: &
!   A0=5.80818, A1=3.20684, A2=4.11890, A3=4.93845, A4=1.01567, A5=1.41575, &
!   B0=-7.01211E-3, B1=-7.25958E-3, B2=7.93334E-3, B3=-5.54491E-3, &
!   C0= -1.32412E-7
!
!double precision :: Ts
!
!Ts = log((298.15-(T-273.15))/(273.15+(T-273.15)))
!
!K0_O2 = A0 + A1*Ts + A2*Ts**2 + A3*Ts**3 + A4*Ts**4 + A5*Ts**5 &
!       +S*( B0 + B1*Ts + B2*Ts**2 + B3*Ts**3 ) &
!       +C0*S**2
!
!K0_O2 = exp( K0_O2 ) / 0.2094  * 1E-6

return
end function K0_O2

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
double precision elemental function K0_CO2( T, S )
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! T in Kelvin
! Equlibrium constant for CO2 solubilty in water, Henry's Law
!  sol       gas
!
! [CO2] = K0 fCO2
!
! fCO2 : fugacity, or approx. the partial pressure of CO2
! Zeebe eq A.3.6 pg 256
! K0 from Weiss 1974
!
! the units are
!   CO2 : mol / kg_soln
!  fCO2 : atm
double precision,intent(in) :: T,S


K0_CO2 = 9345.17/T - 60.2409 + 23.3585*log(T/100.)    &
          + S*( 0.023517 - 0.00023656*T + 0.0047036*(T/100.)**2 )
K0_CO2 = exp( K0_CO2 )

return
end function K0_CO2

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
elemental double precision function K1_H2CO3( T, S, const )
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
double precision,intent(in) :: T,S
integer,optional,intent(in):: const    ! which constant set to use
! First dissociation constant of Carbonic acis
!  [H2CO3] = K_H2CO2_1 [H+] [HCO3-]

integer :: const_i

const_i = 10
if( present(const) ) const_i = const
!
! Zeebe pg 255

! Using the Mehrbach parameterisation  : pH_SWS
! K1_H2CO3 = 3670.7/T - 62.008 + 9.7944*log(T) - 0.0118*S   &
!             + 0.000116*S**2

select case( const_i )
case( 8) ! K1, K2 for fresh water
   K1_H2CO3 = Exp(290.9097 - 14554.21 / T - 45.0575 * log(T))

case(10)! Using the Lueker parameterisation : pH_T
   K1_H2CO3 = 3633.86/T - 61.2172 + 9.6777*log(T) - 0.011555*S + 0.0001152*S**2
   K1_H2CO3 = 10**(-K1_H2CO3)
end select

return
end function K1_H2CO3

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
elemental double precision function K2_H2CO3( T, S, const )
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
double precision,intent(in) :: T,S
integer,optional,intent(in) :: const    ! which constant set to use
! First dissociation constant of Carbonic acis
!     [HCO3-] = K_H2CO2_2 [H+] [CO3--]
!
! Zeebe pg 255
integer :: mode

! Using the Mehrbach parameterisation : pH_SWS
! K2_H2CO3 = 1394.7/T + 4.777 -0.0184*S + 0.000118*S**2
mode = 10
if( present(const) ) mode = const
select case (mode)
case( 8)
!               Millero, F. J., Geochemica et Cosmochemica Acta 43:1651-1661, 1979:
!               K1 from refit data from Harned and Davis,
!                       J American Chemical Society, 65:2030-2037, 1943.
!               K2 from refit data from Harned and Scholes,
!                       J American Chemical Society, 43:1706-1709, 1941.
!       These are the thermodynamic constants:
   K2_H2CO3 = exp(207.6548 - 11843.79 / T - 33.6485 * log(T))

case(10)

! Using the Lueker parameterisation : pH_T
  K2_H2CO3 = 471.8/T + 25.9290 - 3.16967*log(T) - 0.01781*S + 0.0001122*S**2
  K2_H2CO3 = 10**(-K2_H2CO3)
end select

return
end function K2_H2CO3

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
double precision elemental function K_S( T, S )
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! The equlibrium constant for Bisulphate ion
!     K_SO4 [HSO4-] = [H+]_F [SO4--]
!
! Zeebe p259
! pH scale is free scale , pH_F
! This is used to convert 
!    pH_F  <->   pH_T  <->     pH_SWS
!   'free'     Hansson     seawater scale
double precision,intent(in) :: T,S
double precision            :: I    ! ionic strength

I =  19.924*S              &
      / (1000. - 1.005*S )    

K_S =  -4276.1/T + 141.328 -  23.093*log(T)                  &
     + (-13856./T + 324.57  -  47.986*log(T) ) * sqrt(I)     &
     + ( 35474./T - 771.54  + 114.723*log(T) ) * I           &
         -2698./T * I**1.5 + 1776./T * I**2

K_S = exp(K_S) * (1 - 0.001005 * S)


return
end function K_S

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
double precision elemental function T_S( S )
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! total sulphate concentration in seawater as a function of salinity
!
! Units mol/kg-soln
double precision,intent(in) :: S
! Zeebe pg 260
  T_S = 0.02824 * S / 35.0             ! 8.0686E-04 * S

! Morris & Riley (1966) from previous calchem version
! T_S  = 0.14 * S /96.062 / 1.80655    ! 8.0673E-04 * S

! From CO2sys
! T_S = 19.924 * S / (1000 - 1.005 * S)
T_S = (0.14 / 96.062) * (S / 1.80655)
return
end function T_S

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
double precision elemental function K_F( T, S )
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! equilibrium constant for 
!  [HF] K_F = [H+] [F-] 
!
! from Zeebe p 260
! pH scale is Hansson, pH_T
!
double precision,intent(in) :: T, S
double precision :: I, KS, TS, lnKF


!I  =     19.924*S      &
!     /(1000. - 1.005*S)
!
!TS = T_S(S)
!KS = K_S(T,S)
!
!K_F = 1590.2/T - 12.641 + 1.525*sqrt(I) 
!                                   !  this term converts 
!                                   !  'free' -> Hansson
!K_F = exp(K_F) * (1.0 - 0.001005*S)* (1.+ TS/KS)

!  Perez,F.F. and Fraga, F., Marine Chemistry 21(2):161-168, 1987
!S =10-40   T=9-33 oC
lnKF = -874. / T - 0.111 * sqrt(S)+ 9.68
lnKF = -lnKF
K_F = exp(lnKF)

return
end function K_F

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
double precision elemental function T_F( S )
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! total fluoride concentration in seawater
double precision,intent(in) :: S  

! Zeebe p 261                          
!T_F = 7.E-5 * S / 35.0                ! 2.0000E-06 * S            

!  Riley (1965) ,previous version of calchem          
 T_F = 0.000067 * S / 18.9984 / 1.80655 ! 1.9521E-06 * S

return
end function T_F

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
double precision elemental function K1_P( T, S )
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! equilibrium constant for 
!  [H3PO4] K1_P = [H+] [H2PO4-]
!
! from Zeebe p 264
! pH scale is Hansson, pH_T
!
double precision,intent(in) :: T, S

K1_P = -4576.752/T +115.54  -18.453*log(T) &
      +(-106.736/T+0.69171)*sqrt(S) + (-0.65643/T-0.01844)*S

K1_P = exp(K1_P)

return
end function K1_P

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
double precision elemental function K2_P( T, S )
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! equilibrium constant for 
!  [H2PO4-] K2_P = [H+] [HPO4--]
!
! from Zeebe p 264
! pH scale is Hansson, pH_T
!
double precision,intent(in) :: T, S

K2_P = -8814.715/T + 172.1033 -27.927*log(T) &
      +(-160.34/T+1.3566)*sqrt(S) + (0.37355/T-0.05778)*S

K2_P = exp(K2_P)

return
end function K2_P

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
double precision elemental function K3_P( T, S )
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! equilibrium constant for 
!  [HPO4--] K3_P = [H+] [PO4---]
!
! from Zeebe p 265
! pH scale is Hansson, pH_T
!
double precision,intent(in) :: T, S

K3_P = -3070.75/T -18.126 + ( 17.27039/T + 2.81197 )*sqrt(S) &
        + (-44.99486/T - 0.09984)*S

K3_P = exp(K3_P)

return
end function K3_P

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
double precision elemental function K_Si( T, S )
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! equilibrium constant for 
!  [Si(OH)4] K_Si = [H+] [H3SiO4-]
!
! from Zeebe p 265
! pH scale is Hansson, pH_T
!
double precision,intent(in) :: T, S

double precision :: I


I  =     19.924*S      &
     /(1000. - 1.005*S)

!K_Si = -8904.2/T + 117.385 -19.334*log(T)                         &
K_Si = -8904.2/T + 117.4   -19.334*log(T)                         &   
       + ( 3.5913  - 458.79 /T )*sqrt(I) + ( 188.74/T - 1.5998)*I   &
       + ( 0.07871 - 12.1652/T )*I**2  

K_Si = exp(K_Si)*(1.0-0.001005*S)

return
end function K_Si

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
double precision elemental function T_Ca( S )
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! total calcium concentration in seawater
! From Zeebep 269
double precision,intent(in) :: S
T_Ca = 0.01028 * S / 35.0
return
end function T_Ca

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
double precision elemental function K_W( T, S, const )
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Ion product of water
!     K_W = [H+][OH-]
! 
! Zeebe p259
! pH scale is seawater scale, pH_T
double precision,intent(in) :: T,S
integer,optional,intent(in) :: const


double precision :: logT
integer :: const_i

logT=log(T)

const_i = 10
if( present(const) ) const_i = const

select case( const_i )
case(10)
   K_W =  148.9802 - 13847.26/T  - 23.6521 * logT     &
         +(118.67/T  - 5.977 + 1.0495*logT)*sqrt(S) - 0.01615*S

   K_W = exp(K_W)
case( 8 )
   K_W = 148.9802 - 13847.26 / T - 23.6521 * logT
   K_W = exp(K_W)
end select
return
end function K_W

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
double precision elemental function K_B( T, S )
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! the equlibrium constant for boric acid
!
!  K_B [B(OH)_3] = [ H+ ][ B(OH)_4^-]
!
!  From Zeebe p262
!  pH scale is Hansson, pH_T
double precision,intent(in) :: S, T

double precision :: TS,TF,KS,KF

K_B = ( -8966.90 - 2890.53*sqrt(S) - 77.942*S + 1.728*S*sqrt(S)   &
       -0.0996*S**2 ) / T                                         &
       +148.0248 + 137.1942*sqrt(S) + 1.62142*S                   &
       -(24.4344 + 25.085*sqrt(S) + 0.2474*S)*log(T)              &
      +0.053105*sqrt(S)*T
! Needed for converting from Hansson pH_T -> seawater pH_SWS
TS = T_S( S ) ; KS = K_S(T,S)         ! pH_F 
TF = T_F( S ) ; KF = K_F(T,S)         ! pH_T      
K_B = exp( K_B )* (1. + TS/KS)/( 1. + TS/KS + TF/KF )
return
end function K_B

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
double precision elemental function T_B( S ) 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Total Boron mol / kgsln
! Uppstrom (1974)
! ratio of [B mg/kg ]/[Cl parts per 1000 ] =  0.232
double precision,intent(in) :: S

double precision            :: S_cl
! from footnote 3 of Zeebe page 100
!     S = 1.80655 * Cl  where Cl is Chlorinity
!scl = s/1.80655 ! = S * 0.5524412  ! and percentage of Cl by weight is 55% according to wiki
!bt  = 0.000232 * scl/10.811        ! 10.811 is the atomic weight of Boron
                                    ! factor = 1.18787E-5
      
!! From Zeebe p263, Correction from COSys.m matlab file
!T_B = 4.157E-4 * S / 35.d0   ! factor = 1.18857E-5

! Lee, Kim, Byrne, Millero, Feely, Yong-Ming Liu. 2010.  Geochimica Et Cosmochimica Acta 74 (6)
T_B = 0.0004326 * S / 35

return
end function T_B

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
double precision elemental function Ksp_cal( T, S )
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! The solubity product of calcite
!
! Zebee p266
double precision, intent(in) :: T, S

Ksp_cal = -171.9065  -0.077993*T +2839.319/T                &
           + 71.595*log10(T)                                &
           +( -0.77712   +0.0028426*T+178.34/T )*sqrt(S)    &
              -0.07711*S +0.0041249*S**(1.5)
Ksp_cal = 10**(Ksp_cal)

return
end function Ksp_cal

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
double precision elemental function Ksp_arg( T, S )
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! The solubity product of calcite
!
! Zebee p266
double precision, intent(in) :: T, S

Ksp_arg = -171.945 -0.077993*T +2903.293/T                  &
           + 71.595*log10(T)                                &
           +( -0.068393  +0.0017276*T +88.135/T )*sqrt(S)   &
              -0.10018*S +0.0059415*S**(1.5)
Ksp_arg = 10**(Ksp_arg)

return
end function Ksp_arg

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function dCO2dH_dDdA( CO2, H, T, Sal )
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! partial derivative of CO2 and H+ in terms of DIC and Alk
double precision,intent(in) :: H, CO2, T, Sal
double precision            ::  dCO2dH_dDdA(2,2)

double precision            :: K1, K2, KW, KB, TB
double precision :: s,           & ! CO2 conc
          DIC, Ds, Dh, &
          A,   As, Ah

s = CO2

! Calculate the equilibrium constants
K1 = K1_H2CO3( T,Sal )
K2 = K2_H2CO3( T,Sal )
KW = K_W     ( T,Sal )
TB = T_B(Sal)
KB = K_B     ( T,Sal )

! Take the partial derivatives
Ds  = 1. + K1/H    +    K1*K2/H**2
Dh  = -s*( K1/H**2 + 2.*K1*K2/H**3 )

As  =      K1/H    + 2.*K1*K2/H**2
Ah  = -s*( K1/H**2 + 4.*K1*K2/H**3 ) - KB*TB/( KB+H)**2 - KW/H**2 -1.


! Assemble the inverse
dCO2dH_dDdA = reshape( (/ &
            Ah, -As,      &
           -Dh,  Ds /),(/2,2/))/( Ds*Ah-Dh*As)

return
end function dCO2dH_dDdA

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!function Revelle_factor2( H, CO2, T, Sal, TB ) result( RF0 )
 function Revelle_factor2( CO2, H, T, Sal     ) result( dCO2_dH )
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! From Zeebe p71 
double precision,intent(in) :: H, CO2, T, Sal
double precision            :: K1, K2, KW, KB, TB
double precision            :: RF0, dCO2_dH(2)

double precision :: s,           & ! CO2 conc
          DIC, Ds, Dh, &
          A,   As, Ah

K1 = K1_H2CO3( T,Sal )
K2 = K2_H2CO3( T,Sal )
KW = K_W     ( T,Sal )
TB = T_B(Sal)
KB = K_B     ( T,Sal )

!write(6,*) TB_in,  TB_in*KB/(KB+H)

s   = CO2

DIC = s*(1. + K1/H +    K1*K2/H**2 )

A   = s*(   + K1/H + 2.*K1*K2/H**2 ) + TB*KB/(KB+H) + Kw/H - H


Ds  = 1. + K1/H    +    K1*K2/H**2
Dh  = -s*( K1/H**2 + 2.*K1*K2/H**3 )


As  =      K1/H    + 2.*K1*K2/H**2
Ah  = -s*( K1/H**2 + 4.*K1*K2/H**3 ) - KB*TB/( KB+H)**2 - KW/H**2 -1.

!write(6,*) A, DIC
!write(6,*) Ds , Dh*As/Ah
!RF0 =  DIC/CO2 /(Ds - Dh*As/Ah )
dCO2_dH = (/ 1.d0/(Ds - Dh*As/Ah ),   -As/Ah  /)

!RF0 = 1.d0/RF0
     
return
end function Revelle_factor2

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
subroutine CC_run_checks()
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! test values from Zeebe
double precision :: T, S
double precision :: pCO2, fCO2,& ! partial pressure and fugacity of CO2
           CO2,      & ! Dissolved Carbon dioxide   : aqueous concentration of CO2
           DIC,      & ! Dissolved Inorganic Carbon : 
            AT,      & ! Total alkalinity
          HCO3,      & ! Bicarbonate ion            : aqueous concentration of HCO_3^{-} 
           CO3,      & ! Carbonate  ion             : aqueous concentration of  CO_3^{-2} 
         pH, H, OH,  & ! Hydrogen ion               : aqueous concentration of H^{+}
         Omega_cal,  & ! Saturation state calcite
         Omega_arg,  & ! Saturation state aragonite
         TP, AP,     & ! Total Phosphate, phosphate alkinity
         AB, ASi, TSi,&! Borate and Silicate alkalinity
         SWS_2_T, KF,TF,KS,TS

T = 25. +273.15
S = 35.
TS = T_S( S ) ; KS = K_S(T,S)         ! pH_F 
TF = T_F( S ) ; KF = K_F(T,S)         ! pH_T
SWS_2_T = (1. + TS/KS)/( 1. + TS/KS + TF/KF )
write(6,*) 
write(6,*) "Check of rate constants from Zeebe"
write(6,*)

2 format( a12,2(" |",2(1x,f8.4))," |")

write(6,'(27(" -"))')
write(6,'(a)')"Variable  pg |       p Var       |      ln Var       |"
write(6,'(a)')"             | this        Zeebe |  this       Zeebe |"
write(6,'(27(" -"))')
write(6,2)"K0_CO2  257", -log10( K0_CO2  (T,S) ),  1.5468,log(K0_CO2(T,S)), -3.5617
write(6,2)"K1_CO3  255", -log10( K1_H2CO3(T,S) ),  5.8472
write(6,2)"K2_CO3  255", -log10( K2_H2CO3(T,S) ),  8.9660
write(6,2)"K_B     262", -log10( K_B     (T,S) ),  8.5975,log(K_B   (T,S)),-19.7964
write(6,2)"K_S     260", -log10( K_S     (T,S) ),  0.9987,log(K_S   (T,S)), -2.2996
write(6,2)"K_F     261", -log10( K_F     (T,S) ),  2.5183,log(K_F   (T,S)), -5.7986
write(6,2)"K1_P    264", -log10( K1_P    (T,S) ),  1.61  ,log(K1_P  (T,S)), -3.71
write(6,2)"K2_P    264", -log10( K2_P    (T,S) ),  5.96  ,log(K2_P  (T,S)),-13.727
write(6,2)"K3_P    265", -log10( K3_P    (T,S) ),  8.79  ,log(K3_P  (T,S)),-20.24  
write(6,2)"K_Si    265", -log10( K_Si    (T,S) ),  9.38  ,log(K_Si  (T,S)),-21.61
write(6,2)"Ksp_cal    ", -log10( Ksp_cal (T,S) ),  6.3693
write(6,2)"Ksp_arg    ", -log10( Ksp_arg (T,S) ),  6.1883

T = 5. + 273.15
S = 35.

write(6,*)
write(6,*) "Equation of state for seawater"
write(6,*)
write(6,*)"rho(5, 0,0) ",rho_pw(T  )
write(6,*)"rho(5,35,0) ",rho_sw(T,S)

write(6,*) " - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"
write(6,*) "Check solver against CO2sys the excel version 2.0"
write(6,*) 
write(6,*) "CO2 constants: K1, K2 from Leker et al 2000"
write(6,*) "KHSO4 Dickson (1990b)"
write(6,*) "pH scale: Hansson, Total scale"
1 format( (a14,2("|",1x,f12.3,1x),"|",a20,2x,"|") )

S   = 35
T   = 20 + 273.15
pCO2 = 360.d-6
AT  = 2400.d-6

fCO2 = pCO2 * CO2_fugacity_const(T)
CO2  = fCO2 * K0_CO2(T,S)

call CC_solve_CO2_Talk( CO2, AT, T, S, DIC, H,   &
                        HCO3=HCO3, CO3=CO3, pH=pH,&
                        AB = AB, &
                        Omega_arg=Omega_arg,        &
                        Omega_cal=Omega_cal )
write(6,*)
write(6,*) "TEST : solve using Total alkalinity and partial pressure CO2"
write(6,*)
write(6,'(34(" -"))')
write(6,1) "* Salinity    ",S         ,   35.,    "psu"
write(6,1) "* Temperature ",T-273.15  ,   20.,    "deg C"
write(6,1) "  fCO2        ",fCO2 *1.d6,  361.111,  "mu atm"
write(6,1) "* pCO2        ",pCO2 *1.d6,  359.99 ,  "mu atm"
write(6,1) "* AT          ",AT   *1.d6, 2400.00 ,  "mu mol / kg-sol"
write(6,'(34(" -"))')
write(6,1) "  pH          ",pH        ,    8.100,   "pH T"
write(6,1) "  DIC         ",DIC  *1.d6, 2103.86 ,"mu mol / kg-sol"
write(6,1) "  HCO3        ",HCO3 *1.d6, 1879.665,"mu mol / kg-sol"
write(6,1) "   CO3        ", CO3 *1.d6,  212.571,"mu mol / kg-sol"
write(6,1) "   CO2        ", CO2 *1.d6,   11.627,"mu mol / kg-sol"
write(6,1) "  Alk B       ",  AB *1.d6,   90.390,"mu mol / kg-sol"
write(6,1) "  Omega_cal   ", Omega_cal,    5.084," 1"
write(6,1) "  Omega_arg   ", Omega_arg,    3.304," 1"
write(6,'(34(" -"))')


S    = 36
T    = 28 + 273.15
fCO2 =  360.d-6
AT   = 2500.d-6

pCO2 = fCO2 / CO2_fugacity_const(T)
CO2  = fCO2 * K0_CO2(T,S)

call CC_solve_CO2_Talk( CO2, AT, T, S, DIC, H,   &
                        HCO3=HCO3, CO3=CO3, pH=pH,&
                        Omega_arg=Omega_arg,        &
                        Omega_cal=Omega_cal, AB=AB )
write(6,*)
write(6,*) "TEST : solve using Total alkalinity and fugacity CO2"
write(6,*)
write(6,'(34(" -"))')
write(6,1) "* Salinity    ",S         , 36.,     "psu"
write(6,1) "* Temperature ",T-273.15  , 28.,     "deg C"
write(6,1) "* fCO2        ",fCO2 *1.d6, 360.000, "mu atm"
write(6,1) "  pCO2        ",pCO2 *1.d6, 361.111, "mu atm"
write(6,1) "* AT          ",AT   *1.d6, 2500.,   "mu mol / kg-sol"
write(6,'(34(" -"))')
write(6,1) "  pH          ",pH        , 8.102,   "pH T"
write(6,1) "  DIC         ",DIC  *1.d6, 2105.208,"mu mol / kg-sol"
write(6,1) "  HCO3        ",HCO3 *1.d6, 1812.654,"mu mol / kg-sol"
write(6,1) "   CO3        ", CO3 *1.d6,  283.106,"mu mol / kg-sol"
write(6,1) "   CO2        ", CO2 *1.d6,    9.448,"mu mol / kg-sol"
write(6,1) "  Alk B       ",  AB *1.d6,  111.033,"mu mol / kg-sol"
write(6,1) "  Omega_cal   ", Omega_cal,    6.792," 1"
write(6,1) "  Omega_arg   ", Omega_arg,    4.526," 1"


S   = 37
T   = 22 + 273.15
AT  = 2350.d-6
DIC = 2100.d-6
call CC_solve_DIC_Talk       ( DIC, AT, T, S, CO2, H,   &
                        HCO3=HCO3, CO3=CO3, pH=pH,&
                        Omega_arg=Omega_arg,        &
                        Omega_cal=Omega_cal, AB=AB )
fCO2 =  CO2 / K0_CO2(T,S)
pCO2 = fCO2 / CO2_fugacity_const(T)

write(6,*)
write(6,*) "TEST : solve using Total alkalinity and partial pressure CO2"
write(6,*)
write(6,'(34(" -"))')
write(6,1) "* Salinity    ",S         , 37.,     "psu"
write(6,1) "* Temperature ",T-273.15  , 22.,     "deg C"
write(6,1) "  fCO2        ",fCO2 *1.d6, 493.455, "mu atm"
write(6,1) "  pCO2        ",pCO2 *1.d6, 495.093, "mu atm"
write(6,1) "* AT          ",AT   *1.d6, 2350.,   "mu mol / kg-sol"
write(6,'(34(" -"))')
write(6,1) "  pH          ",pH        , 7.971,    "pH SWS"
write(6,1) "* DIC         ",DIC  *1.d6, 2100.,   "mu mol / kg-sol"
write(6,1) "  HCO3        ",HCO3 *1.d6, 1904.703,"mu mol / kg-sol"
write(6,1) "   CO3        ", CO3 *1.d6,  180.313,"mu mol / kg-sol"
write(6,1) "   CO2        ", CO2 *1.d6,   14.985,"mu mol / kg-sol"
write(6,1) "  Alk B       ",  AB *1.d6,   80.241,"mu mol / kg-sol"
write(6,1) "  Omega_cal   ", Omega_cal,    4.254," 1"
write(6,1) "  Omega_arg   ", Omega_arg,    2.787," 1"
write(6,'(34(" -"))')


S   = 34
T   = 18 + 273.15
AT  = 2350.d-6
TP  =  105.d-6
DIC = 2100.d-6
call CC_solve_DIC_Talk( DIC, AT, T, S, CO2, H,     &
                        TP = TP,                   &
                        HCO3=HCO3, CO3=CO3, OH=OH, pH=pH, &
                        Omega_arg=Omega_arg,       &
                        Omega_cal=Omega_cal,       &
                        AB=AB, AP=AP )
fCO2 =  CO2 / K0_CO2(T,S)
pCO2 = fCO2 / CO2_fugacity_const(T)

write(6,*)
write(6,*) "TEST : solve using Total alkalinity and DIC and TP"
write(6,*)
write(6,'(34(" -"))')
write(6,1) "* Salinity    ",S         , 34.,     "psu"
write(6,1) "* Temperature ",T-273.15  , 18.,     "deg C"
write(6,1) "  fCO2        ",fCO2 *1.d6, 679.938, "mu atm"
write(6,1) "  pCO2        ",pCO2 *1.d6, 682.310, "mu atm"
write(6,1) "* AT          ",AT   *1.d6, 2350.,   "mu mol / kg-sol"
write(6,1) "* TP          ",TP   *1.d6,  105.,   "mu mol / kg-sol"
write(6,'(34(" -"))')
write(6,1) "  pH          ",pH        , 7.838,   "pH SWS"
write(6,1) "* DIC         ",DIC  *1.d6, 2100.,   "mu mol / kg-sol"
write(6,1) "  HCO3        ",HCO3 *1.d6, 1966.432,"mu mol / kg-sol"
write(6,1) "   CO3        ", CO3 *1.d6,  110.127,"mu mol / kg-sol"
write(6,1) "   CO2        ", CO2 *1.d6,   23.442,"mu mol / kg-sol"
write(6,1) "  Alk B       ",  AB *1.d6,   50.181,"mu mol / kg-sol"
write(6,1) "   OH         ",  OH *1.d6,    2.137,"mu mol / kg-sol"
write(6,1) "  Alk P       ",  AP *1.d6 , 111.004,"mu mol / kg-sol" 
write(6,1) "  Omega_cal   ", Omega_cal,    2.648," 1"
write(6,1) "  Omega_arg   ", Omega_arg,    1.709," 1"
write(6,'(34(" -"))')


S   = 34
T   = 18 + 273.15
AT  = 2350.d-6
TSi =  110.d-6
DIC = 2100.d-6
call CC_solve_DIC_Talk( DIC, AT, T, S, CO2, H,     &
                         TSi=TSi, TP=0.d0,         &
                        HCO3=HCO3, CO3=CO3, OH=OH, pH=pH, &
                        Omega_arg=Omega_arg,       &
                        Omega_cal=Omega_cal,       &
                        AB=AB, AP=AP, ASi=ASi )
fCO2 =  CO2 / K0_CO2(T,S)
pCO2 = fCO2 / CO2_fugacity_const(T)

write(6,*)
write(6,*) "TEST : solve using Total alkalinity and DIC and TP"
write(6,*)
write(6,'(34(" -"))')
write(6,1) "* Salinity    ",S         , 34.,     "psu"
write(6,1) "* Temperature ",T-273.15  , 18.,     "deg C"
write(6,1) "  fCO2        ",fCO2 *1.d6, 394.214, "mu atm"
write(6,1) "  pCO2        ",pCO2 *1.d6, 395.589, "mu atm"
write(6,1) "* AT          ",AT   *1.d6, 2350.,   "mu mol / kg-sol"
write(6,1) "* TSi         ",TSi  *1.d6,  110.,   "mu mol / kg-sol"
write(6,'(34(" -"))')
write(6,1) "  pH          ",pH        , 8.061,   "pH SWS"
write(6,1) "* DIC         ",DIC  *1.d6, 2100.,   "mu mol / kg-sol"
write(6,1) "  HCO3        ",HCO3 *1.d6, 1907.658,"mu mol / kg-sol"
write(6,1) "   CO3        ", CO3 *1.d6,  178.751,"mu mol / kg-sol"
write(6,1) "   CO2        ", CO2 *1.d6,   13.591,"mu mol / kg-sol"
write(6,1) "  Alk B       ",  AB *1.d6,   77.483,"mu mol / kg-sol"
write(6,1) "   OH         ",  OH *1.d6,    3.584,"mu mol / kg-sol"
write(6,1) "  Alk Si      ",  ASi*1.d6 ,   3.783,"mu mol / kg-sol"
write(6,1) "  Omega_cal   ", Omega_cal,    4.298," 1"
write(6,1) "  Omega_arg   ", Omega_arg,    2.775," 1"
write(6,'(34(" -"))')

S   = 33.63
T   = 278
AT  = 2340.d-6
DIC = 2162.d-6
call CC_solve_DIC_Talk( DIC, AT, T, S, CO2, H,     &
                        TSi=Tsi,                   &
                        HCO3=HCO3, CO3=CO3, OH=OH, pH=pH, &
                        Omega_arg=Omega_arg,       &
                        Omega_cal=Omega_cal,       &
                        AB=AB, ASi=ASi  )
fCO2 =  CO2 / K0_CO2(T,S)
pCO2 = fCO2 / CO2_fugacity_const(T)

write(6,*)
write(6,*) "TEST : solve using Total alkalinity and DIC and TP"
write(6,*)
write(6,'(34(" -"))')
write(6,1) "* Salinity    ",S         , 34.,     "psu"
write(6,1) "* Temperature ",T-273.15  , 18.,     "deg C"
write(6,1) "  fCO2        ",fCO2 *1.d6, 394.214, "mu atm"
write(6,1) "  pCO2        ",pCO2 *1.d6, 395.589, "mu atm"
write(6,1) "* AT          ",AT   *1.d6, 2350.,   "mu mol / kg-sol"
write(6,'(34(" -"))')
write(6,1) "  pH          ",pH        , 8.061,   "pH SWS"
write(6,1) "* DIC         ",DIC  *1.d6, 2100.,   "mu mol / kg-sol"
write(6,1) "  HCO3        ",HCO3 *1.d6, 1907.658,"mu mol / kg-sol"
write(6,1) "   CO3        ", CO3 *1.d6,  178.751,"mu mol / kg-sol"
write(6,1) "   CO2        ", CO2 *1.d6,   13.591,"mu mol / kg-sol"
write(6,1) "  Alk B       ",  AB *1.d6,   77.483,"mu mol / kg-sol"
write(6,1) "   OH         ",  OH *1.d6,    3.584,"mu mol / kg-sol"
write(6,1) "  Omega_cal   ", Omega_cal,    4.298," 1"
write(6,1) "  Omega_arg   ", Omega_arg,    2.775," 1"
write(6,'(34(" -"))')

return
end subroutine CC_run_checks
end module carbonate_chemistry_mod

