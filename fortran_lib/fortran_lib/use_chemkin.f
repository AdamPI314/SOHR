      Subroutine chemkinInitialize()
      implicit none
      include "chemkin/ckvariables.f"
      character*106 ctmp

      OPEN (LINKCK,FORM='UNFORMATTED',STATUS='UNKNOWN',
     >FILE='./input/chem.bin')
      CALL CKLEN(LINKCK,LOUT,LENICK,LENRCK,LENCCK)
      CALL CKINIT(LENICK,LENRCK,LENCCK,LINKCK,LOUT,ICKWRK,
     >RCKWRK,CCKWRK)
      CLOSE (LINKCK)

      open(unit=51,file="./output/general_output.out")
      write(51,"(A)") "General output:"
      write(51,*) "   "
      write(51,*) " Number of Elements (NMM) = ", nmm
      write(51,*) " Number of Species (NKK) = ", nkk
      write(51,*) " Number of Reactions (NII) = ", nii
      write(51,*) ""

      ctmp="  Please make sure that in the file 'global_extern_vars.h"
     >//"', the value of the constant 'NumberOfReactions'"
      write(51,"(A)") ctmp
      ctmp="  is  the same as 'Number of Reactions (NII)', and the"
     >//" value of the constant 'NumberOfSpecies' is the"
      write(51,"(A)") ctmp
      ctmp="  same as 'Number of Species (NKK)'."
      write(51,"(A)") ctmp

      write(51,*) " LENICK = ", LENICK
      write(51,*) " LENRCK = ", LENRCK
      write(51,*) " LENCCK = ", LENCCK
      write(51,*) "   "

      ctmp="  Please make sure that the values of constants 'LENIWK'"
     >//", 'LENRWK' and 'LENCWK' in the file "
      write(51,"(A)") ctmp
      ctmp="  'global_extern_vars.h"
     >//"' are the same as those in the file"
     >//" 'ckvariables.f'. Also, they must"
      write(51,"(A)") ctmp
      ctmp="  be larger than"
     >//"  'LENICK', 'LENRCK' and 'LENCCK' respectively."
      write(51,"(A)") ctmp

      close(unit=51);

      Return
      End


      Subroutine ChemkinUpdate(neq,t,xgst,derstr)
      ! at constant volume
      ! 'neq' is number of first order ode-s.
      ! 't' is a dummy variable here.
      ! 'xgst' is array of initial values, of length neq. The last value in this array is temperature. The others are mass fraction.
      ! 'derstr' is the first derivative of 'xgst', d[xgst]/dt.
      implicit none
      include "chemkin/ckvariables.f"

      double precision rhomass
      COMMON /ckstore/ rhomass

      integer neq, i
      double precision t, xgst(neq), derstr(neq)
      double precision wdot(nkk), yy(nkk), wt(nkk), h(nkk)
      double precision ttr, p, cvbs, sumh

      ttr=xgst(neq)

      call ckwt(ICKWRK,RCKWRK,wt)
      call ckwyr(rhomass,ttr,xgst,ickwrk,rckwrk,wdot)

      do i=1,nkk
          derstr(i)=wt(i)*wdot(i)/rhomass
      end do

      call ckums(ttr,ickwrk,rckwrk,h)
      call ckcvbs(ttr,xgst,ickwrk,rckwrk,cvbs)
      sumh=0.d0

      do i=1,nkk
          sumh=sumh+wt(i)*wdot(i)*h(i)
      end do

      derstr(neq)=-sumh/(cvbs*rhomass)

      Return
      End

      Subroutine ChemkinUpdateVT(neq,t,xgst,derstr)
      ! at constant volume
      ! at constant temperature
      ! 'neq' is number of first order ode-s.
      ! 't' is a dummy variable here.
      ! 'xgst' is array of initial values, of length neq. The last value in this array is temperature. The others are mass fraction.
      ! 'derstr' is the first derivative of 'xgst', d[xgst]/dt.
      implicit none
      include "chemkin/ckvariables.f"

      double precision rhomass
      COMMON /ckstore/ rhomass

      integer neq, i
      double precision t, xgst(neq), derstr(neq)
      double precision wdot(nkk), yy(nkk), wt(nkk), h(nkk)
      double precision ttr, p, cvbs, sumh

      ttr=xgst(neq)

      call ckwt(ICKWRK,RCKWRK,wt)
      call ckwyr(rhomass,ttr,xgst,ickwrk,rckwrk,wdot)

      do i=1,nkk
          derstr(i)=wt(i)*wdot(i)/rhomass
      end do

      derstr(neq)=0.0

      Return
      End

      Subroutine ChemkinUpdateP(neq,t,xgst,derstr)
      ! at constant pressure
      ! 'neq' is number of first order ode-s.
      ! 't' is a dummy variable here.
      ! 'xgst' is array of initial values, of length neq. The last value in this array is temperature. The others are mass fraction.
      ! 'derstr' is the first derivative of 'xgst', d[xgst]/dt.
      implicit none
      include "chemkin/ckvariables.f"

      double precision rhomass, pressure
      COMMON /ckstore/ rhomass, pressure

      integer neq, i
      double precision t, xgst(neq), derstr(neq)
      double precision wdot(nkk), yy(nkk), wt(nkk), h(nkk)
      double precision ttr, cpbs, sumh, rhomass_inverse

      ttr=xgst(neq)

      ! define mass fractions
      call ckwt(ICKWRK,RCKWRK,wt)
      ! call ckwyr(rhomass,ttr,xgst,ickwrk,rckwrk,wdot)
      call ckwyp(pressure,ttr,xgst,ickwrk,rckwrk,wdot)
      call ckrhoy(pressure, ttr, xgst, ickwrk, rckwrk, rhomass)

      ! multiplication is faster than division
      rhomass_inverse= 1.d0/rhomass
      do i=1,nkk
          derstr(i)=wt(i)*wdot(i)*rhomass_inverse
      end do

      ! temperature equation, get enthalpies in mass unit
      !call ckums(ttr,ickwrk,rckwrk,h)
      call ckhms(ttr,ickwrk,rckwrk,h)

      ! mean specific heat
      !call ckcpbs(ttr,xgst,ickwrk,rckwrk,cpbs)
      call ckcpbs(ttr,xgst,ickwrk,rckwrk,cpbs)
      sumh=0.d0

      do i=1,nkk
          sumh=sumh+wt(i)*wdot(i)*h(i)
      end do

      derstr(neq)=-sumh/(cpbs*rhomass)

      Return
      End

      Subroutine ChemkinUpdatePT(neq,t,xgst,derstr)
      ! at constant pressure and constant temperature
      ! 'neq' is number of first order ode-s.
      ! 't' is a dummy variable here.
      ! 'xgst' is array of initial values, of length neq. The last value in this array is temperature. The others are mass fraction.
      ! 'derstr' is the first derivative of 'xgst', d[xgst]/dt.
      implicit none
      include "chemkin/ckvariables.f"

      double precision rhomass, pressure
      COMMON /ckstore/ rhomass, pressure

      integer neq, i
      double precision t, xgst(neq), derstr(neq)
      double precision wdot(nkk), yy(nkk), wt(nkk), h(nkk)
      double precision ttr, cpbs, sumh, rhomass_inverse

      ttr=xgst(neq)

      ! define mass fractions
      call ckwt(ICKWRK,RCKWRK,wt)
      ! call ckwyr(rhomass,ttr,xgst,ickwrk,rckwrk,wdot)
      call ckwyp(pressure,ttr,xgst,ickwrk,rckwrk,wdot)
      call ckrhoy(pressure, ttr, xgst, ickwrk, rckwrk, rhomass)

      ! multiplication is faster than division
      rhomass_inverse= 1.d0/rhomass
      do i=1,nkk
          derstr(i)=wt(i)*wdot(i)*rhomass_inverse
      end do

      ! temperature is constant
      derstr(neq)=0.0

      Return
      End

      Subroutine ChemkinUpdateST(neq,t,ct,derstr)
      ! surface reaction
      ! at constant temperature
      ! 'neq' is number of first order ode-s.
      ! 't' is a dummy variable here.
      ! 'ct' is array of initial values, of length neq. The last value in this array is temperature. The others are molar concentration
      ! 'derstr' is the first derivative of 'ct', d[ct]/dt.
      implicit none
      include "chemkin/ckvariables.f"

      double precision rhomass, pressure
      COMMON /ckstore/ rhomass, pressure

      integer neq, i
      double precision t, ct(neq), derstr(neq)
      double precision wdot(nkk), yy(nkk), wt(nkk), h(nkk)
      double precision ttr, cpbs, sumh, rhomass_inverse

      ttr=ct(neq)

      call ckwc(ttr,ct,ickwrk,rckwrk,wdot)

      do i=1,nkk
          derstr(i)=wdot(i)
      end do

      ! temperature is constant
      derstr(neq)=0.0

      Return
      End

      Subroutine ChemkinUpdateSTCC1(neq,t,ct,derstr)
      ! hold the concentration of the first species to be constant
      ! surface reaction
      ! at constant temperature
      ! 'neq' is number of first order ode-s.
      ! 't' is a dummy variable here.
      ! 'ct' is array of initial values, of length neq. The last value in this array is temperature. The others are molar concentration
      ! 'derstr' is the first derivative of 'ct', d[ct]/dt.
      implicit none
      include "chemkin/ckvariables.f"

      double precision rhomass, pressure
      COMMON /ckstore/ rhomass, pressure

      integer neq, i
      double precision t, ct(neq), derstr(neq)
      double precision wdot(nkk), yy(nkk), wt(nkk), h(nkk)
      double precision ttr, cpbs, sumh, rhomass_inverse

      ttr=ct(neq)

      call ckwc(ttr,ct,ickwrk,rckwrk,wdot)

      do i=1,nkk
          derstr(i)=wdot(i)
      end do

      ! temperature is constant
      derstr(neq)=0.0
      ! the concentration of the first species to be constant
      derstr(1)=0.0

      Return
      End

      Subroutine ChemkinUpdateSTCC2(neq,t,ct,derstr)
      ! hold the concentration of the first two species to be constant
      ! surface reaction
      ! at constant temperature
      ! 'neq' is number of first order ode-s.
      ! 't' is a dummy variable here.
      ! 'ct' is array of initial values, of length neq. The last value in this array is temperature. The others are molar concentration
      ! 'derstr' is the first derivative of 'ct', d[ct]/dt.
      implicit none
      include "chemkin/ckvariables.f"

      double precision rhomass, pressure
      COMMON /ckstore/ rhomass, pressure

      integer neq, i
      double precision t, ct(neq), derstr(neq)
      double precision wdot(nkk), yy(nkk), wt(nkk), h(nkk)
      double precision ttr, cpbs, sumh, rhomass_inverse

      ttr=ct(neq)

      call ckwc(ttr,ct,ickwrk,rckwrk,wdot)

      do i=1,nkk
          derstr(i)=wdot(i)
      end do

      ! temperature is constant
      derstr(neq)=0.0
      ! the concentration of the first two species to be constant
      derstr(1)=0.0
      derstr(2)=0.0

      Return
      End

      Subroutine calculateTdotV(y, t, tdot)
      ! constant volume
      ! calculate temperature derivative at a initial temperature and initial mass fraction

      implicit none
      include "chemkin/ckvariables.f"

      double precision y(nkk)
      double precision t, tdot

      double precision rhomass
      COMMON /ckstore/ rhomass

      integer i
      double precision wdot(nkk), wt(nkk), h(nkk)
      double precision cvbs, sumh

      call ckwt(ICKWRK,RCKWRK,wt)
      call ckwyr(rhomass,t,y,ickwrk,rckwrk,wdot)

      call ckums(t,ickwrk,rckwrk,h)
      call ckcvbs(t,y,ickwrk,rckwrk,cvbs)

      sumh=0.d0
      do i=1,nkk
          sumh=sumh+wt(i)*wdot(i)*h(i)
      end do

      tdot=-sumh/(cvbs*rhomass)
      Return
      End
