#ifdef __GFORTRAN__

      Subroutine canteraInitialize()


      Return
      End


      Subroutine CanteraUpdate(neq,t,xgst,derstr)
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

      Subroutine CanteraUpdateVT(neq,t,xgst,derstr)
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

      Subroutine CanteraUpdateP(neq,t,xgst,derstr)
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

      Subroutine CanteraUpdatePT(neq,t,xgst,derstr)
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

      Subroutine CanteraUpdateST(neq,t,ct,derstr)
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

      Subroutine CanteraUpdateSTCC1(neq,t,ct,derstr)
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

      Subroutine CanteraUpdateSTCC2(neq,t,ct,derstr)
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

#endif
