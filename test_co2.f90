!--------------------------------------------------------------------------------
program ctest
!--------------------------------------------------------------------------------
use carbonate_chemistry_mod
implicit none
character*(300) :: line
double precision :: T,S,  P,TP,TSi, TA,DIC, CO2,CO3, HCO3,pH,H,OH,AP,Asi


S   = 33.0d0
DIC = 2200
TA  = 2300


1 format(22(1x,f14.6))
do 
   read(5,'(a300)',end=5) line
   line = adjustl(line)
   if( line(1:1)=="#" ) cycle
   read(line,*) S,T, P,TP,TSi, TA,DIC
   call CC_solve_DIC_Talk( DIC*1d-6, TA*1d-6, T+273.15, S, CO2, H, pH=pH, CO3=CO3, HCO3=HCO3,&
          TP=TP*1d-6, TSi=TSi*1d-6,OH=OH,AP=AP,ASi=ASi,pHi=7.d0, const=10 )
   write(6,1) T,S, 0.,0.,0.,TA,DIC, pH, 0.,0.,HCO3*1d6, CO3*1d6, CO2*1d6,0.d0,OH*1d6,AP*1d6,ASi*1d6
enddo
5 continue

end program ctest
