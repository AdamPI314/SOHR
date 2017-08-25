#ifdef __GFORTRAN__

      Subroutine canteraInitialize()
      ! Read in the reaction mechanism. Since this is done differently
      ! than in Chemkin, this function does not correspond to any CKLIB
      ! subroutine.
      call newIdealGasMix('./input/chem.xml','gas','')

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
      !only change this routine, minimize change from chemkin --> cantera
      call ctwyr(rhomass,ttr,xgst,ickwrk,rckwrk,wdot)

      do i=1,nkk
          derstr(i)=wt(i)*wdot(i)/rhomass
          write(*,*) wdot(i)
      end do
      stop 1

      call ckums(ttr,ickwrk,rckwrk,h)
      call ckcvbs(ttr,xgst,ickwrk,rckwrk,cvbs)
      sumh=0.d0

      do i=1,nkk
          sumh=sumh+wt(i)*wdot(i)*h(i)
      end do

      derstr(neq)=-sumh/(cvbs*rhomass)

      Return
      End


#endif
