!--------------------------------------------------------------------------------
module py_interface_mod
!--------------------------------------------------------------------------------
use carbonate_chemistry_mod

implicit none

contains

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function CC_solve( input, const ) result( output)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
double precision,intent(in) :: input(:)
integer,intent(in)          :: const
double precision            :: output(14)

double precision :: H 

output=0

                        ! DIC       TA      T         S
call CC_solve_DIC_Talk( input(1), input(2), input(3),input(4),         &
                       output(8), H,                                   &
                       TP=input(5), TSi=input(6),                      &
                       pH=output( 3), HCO3=output( 6), CO3=output( 7),  &
                       pCO2=output(5), &
                       OH=output(10),   AP=output(11), ASi=output(12),  & 
                       const=const)

output(1) = input(2)
output(2) = input(1)
return
end function  CC_solve

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function calc_curve( theta, vdata, aux_data ) result(pH_curve)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
double precision,intent(in) :: theta(:), vdata(:,:), aux_data(:)
double precision            :: pH_curve( size(vdata,2))

double precision :: pHi, vol, S,  DIC,  TA,    mass_s, conc_a, K1_f, K2_f, &
                   mass_a,TK,  S_i,DIC_i,TA_i, H, CO2, pH, rho_a=1.d0, a_pH,b_pH,&
                   TP=0., TSi=0.
                   

integer :: i,np, const=10, naux

a_pH = theta(1)
b_pH = theta(2)
K1_f = theta(3)
K2_f = theta(4)
DIC  = theta(5)*1d-3
TA   = theta(6)*1d-3

naux = size(aux_data)
mass_s = aux_data(1)
conc_a = aux_data(2)
S      = aux_data(3)
const=10 ; if(4<= naux) const = int(aux_data(4))
TP   =0. ; if(5<= naux) TP    =     aux_data(5)*1.d-6
TSi  =0. ; if(6<= naux) TSi   =     aux_data(6)*1.d-6

pHi = 4.0
np  = size(vdata,2)
do i = np,1,-1
   vol= vdata(1,i)
   TK = vdata(3,i)+273.15
   mass_a      = vol*rho_a
   DIC_i       =    mass_s*DIC                 /(mass_s + mass_a )
   TA_i        =  ( mass_s*TA - mass_a*conc_a )/(mass_s + mass_a )
   S_i         =  ( mass_s*S  + mass_a*0.     )/(mass_s + mass_a )
   call CC_solve_DIC_Talk( DIC_i, TA_i, TK, S_i, CO2,H, &
                  pH=pH, &
                  pHi=pHi, &
                  TP=TP, TSi=TSi, &
                  K1_f=K1_f, K2_f=K2_f,          const=const)
   pH_curve(i) = a_pH*(pH - 7.d0)+7.d0 + b_pH
   pHi = pH
enddo

return
end function calc_curve

end module py_interface_mod
