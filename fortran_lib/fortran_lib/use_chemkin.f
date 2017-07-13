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
     >// "', the value of the constant 'NumberOfReactions'"
      write(51,"(A)") ctmp
      ctmp="  is  the same as 'Number of Reactions (NII)', and the"
     >// " value of the constant 'NumberOfSpecies' is the"
      write(51,"(A)") ctmp
      ctmp="  same as 'Number of Species (NKK)'."
      write(51,"(A)") ctmp

      write(51,*) " LENICK = ", LENICK
      write(51,*) " LENRCK = ", LENRCK
      write(51,*) " LENCCK = ", LENCCK
      write(51,*) "   "

      ctmp="  Please make sure that the values of constants 'LENIWK'"
     >// ", 'LENRWK' and 'LENCWK' in the file "
      write(51,"(A)") ctmp
      ctmp="  'global_extern_vars.h"
     >// "' are the same as those in the file"
     >// " 'ckvariables.f'. Also, they must"
      write(51,"(A)") ctmp
      ctmp="  be larger than"
     >//"  'LENICK', 'LENRCK' and 'LENCCK' respectively."
      write(51,"(A)") ctmp

      close(unit=51);

      Return
      End


      Subroutine cppdlsodeV(ti,tout,neq,xgst)
      ! at constant volume
      ! 'ti' is the initial value of the independent variable.
      ! 'tout' is the first point where output is desired (.ne. ti).
      ! 'neq' is number of first order ode-s.
      ! 'xgst' is array of initial values, of length neq. The last value in this array is temperature. The others are mass fraction.
      implicit none
      include "chemkin/ckvariables.f"
      external Fcn11r

      integer itol,mf,itask,istate,iopt,lrw,liw
      double precision ti,tout,dt,atol,rtol,rjac
      !parameter (lrw=500,liw=500)
      !Have to make sure lrw, liw big enough, or it will have error like:
      !*** glibc detected *** ./a.out: double free or corruption (top):
      !*0x08901d70 ***
      !======= Backtrace: =========
      !/lib/libc.so.6(+0x6c501)[0x17c501]
      !/lib/libc.so.6(+0x6dd70)[0x17dd70]
      !/lib/libc.so.6(cfree+0x6d)[0x180e5d]
      !/lib/libc.so.6(fclose+0x14a)[0x16c81a]
      !./a.out[0x8048998]
      !/lib/libpthread.so.0(+0x5cc9)[0xc1fcc9]
      !/lib/libc.so.6(clone+0x5e)[0x1e069e]
      !======= Memory map: ========
      !parameter (lrw=500000,liw=500000)
      parameter (lrw=50000,liw=50000)
      double precision rwork(lrw),iwork(liw)
      integer neq
      double precision xgst(neq)
      common /lsodestore/ atol,rtol,mf,itask,istate,iopt,itol

      call dlsode (Fcn11r,neq,xgst,ti,tout,itol,rtol,atol,
     >itask,istate,iopt,rwork,lrw,iwork,liw,rjac,mf)

      Return
      End
      
      Subroutine cppdlsodeVT(ti,tout,neq,xgst)
      ! at constant volume and constant temperature
      ! 'ti' is the initial value of the independent variable.
      ! 'tout' is the first point where output is desired (.ne. ti).
      ! 'neq' is number of first order ode-s.
      ! 'xgst' is array of initial values, of length neq. The last value in this array is temperature. The others are mass fraction.
      implicit none
      include "chemkin/ckvariables.f"
      external Fcn11rVT

      integer itol,mf,itask,istate,iopt,lrw,liw
      double precision ti,tout,dt,atol,rtol,rjac
      parameter (lrw=50000,liw=50000)
      double precision rwork(lrw),iwork(liw)
      integer neq
      double precision xgst(neq)
      common /lsodestore/ atol,rtol,mf,itask,istate,iopt,itol

      call dlsode (Fcn11rVT,neq,xgst,ti,tout,itol,rtol,atol,
     >itask,istate,iopt,rwork,lrw,iwork,liw,rjac,mf)

      Return
      End
      
      Subroutine cppdlsodeP(ti,tout,neq,xgst)
      ! at constant pressure
      ! 'ti' is the initial value of the independent variable.
      ! 'tout' is the first point where output is desired (.ne. ti).
      ! 'neq' is number of first order ode-s.
      ! 'xgst' is array of initial values, of length neq. The last value in this array is temperature. The others are mass fraction.
      implicit none
      include "chemkin/ckvariables.f"
      external Fcn11rP

      integer itol,mf,itask,istate,iopt,lrw,liw
      double precision ti,tout,dt,atol,rtol,rjac
      parameter (lrw=50000,liw=50000)
      double precision rwork(lrw),iwork(liw)
      integer neq
      double precision xgst(neq)
      common /lsodestore/ atol,rtol,mf,itask,istate,iopt,itol

      call dlsode (Fcn11rP,neq,xgst,ti,tout,itol,rtol,atol,
     >itask,istate,iopt,rwork,lrw,iwork,liw,rjac,mf)

      Return
      End

      Subroutine cppdlsodePT(ti,tout,neq,xgst)
      ! at constant pressure and constant temperature
      ! 'ti' is the initial value of the independent variable.
      ! 'tout' is the first point where output is desired (.ne. ti).
      ! 'neq' is number of first order ode-s.
      ! 'xgst' is array of initial values, of length neq. The last value in this array is temperature. The others are mass fraction.
      implicit none
      include "chemkin/ckvariables.f"
      external Fcn11rPT

      integer itol,mf,itask,istate,iopt,lrw,liw
      double precision ti,tout,dt,atol,rtol,rjac
      parameter (lrw=50000,liw=50000)
      double precision rwork(lrw),iwork(liw)
      integer neq
      double precision xgst(neq)
      common /lsodestore/ atol,rtol,mf,itask,istate,iopt,itol

      call dlsode (Fcn11rPT,neq,xgst,ti,tout,itol,rtol,atol,
     >itask,istate,iopt,rwork,lrw,iwork,liw,rjac,mf)

      Return
      End

      Subroutine cppdlsodeST(ti,tout,neq,ct)
      ! surface reaction
      ! at constant temperature
      ! 'ti' is the initial value of the independent variable.
      ! 'tout' is the first point where output is desired (.ne. ti).
      ! 'neq' is number of first order ode-s.
      ! 'ct' is array of initial values, of length neq. The last value in this array is temperature. The others are molar concentration
      implicit none
      include "chemkin/ckvariables.f"
      external Fcn11rST

      integer itol,mf,itask,istate,iopt,lrw,liw
      double precision ti,tout,dt,atol,rtol,rjac
      parameter (lrw=50000,liw=50000)
      double precision rwork(lrw),iwork(liw)
      integer neq
      double precision ct(neq)
      common /lsodestore/ atol,rtol,mf,itask,istate,iopt,itol

      call dlsode (Fcn11rST,neq,ct,ti,tout,itol,rtol,atol,
     >itask,istate,iopt,rwork,lrw,iwork,liw,rjac,mf)

      Return
      End

      Subroutine cppdlsodaV(ti,tout,neq,xgst)
      ! at constant volume
      ! 'ti' is the initial value of the independent variable.
      ! 'tout' is the first point where output is desired (.ne. ti).
      ! 'neq' is number of first order ode-s.
      ! 'xgst' is array of initial values, of length neq. The last value in this array is temperature. The others are mass fraction.
      ! lsoda, mf should be jt
      implicit none
      include "chemkin/ckvariables.f"
      external Fcn11r

      integer itol,mf,itask,istate,iopt,lrw,liw,jt
      double precision ti,tout,dt,atol,rtol,jdum

      parameter (lrw=50000,liw=50000)
      double precision rwork(lrw),iwork(liw)
      integer neq
      double precision xgst(neq)
      common /lsodestore/ atol,rtol,mf,itask,istate,iopt,itol,jt

      call dlsoda (Fcn11r,neq,xgst,ti,tout,itol,rtol,atol,
     >itask,istate,iopt,rwork,lrw,iwork,liw,jdum,jt)

      Return
      End

      Subroutine cppdlsodaVT(ti,tout,neq,xgst)
      ! at constant volume and temperature
      ! 'ti' is the initial value of the independent variable.
      ! 'tout' is the first point where output is desired (.ne. ti).
      ! 'neq' is number of first order ode-s.
      ! 'xgst' is array of initial values, of length neq. The last value in this array is temperature. The others are mass fraction.
      ! lsoda, mf should be jt
      implicit none
      include "chemkin/ckvariables.f"
      external Fcn11rVT

      integer itol,mf,itask,istate,iopt,lrw,liw,jt
      double precision ti,tout,dt,atol,rtol,jdum

      parameter (lrw=50000,liw=50000)
      double precision rwork(lrw),iwork(liw)
      integer neq
      double precision xgst(neq)
      common /lsodestore/ atol,rtol,mf,itask,istate,iopt,itol,jt

      call dlsoda (Fcn11rVT,neq,xgst,ti,tout,itol,rtol,atol,
     >itask,istate,iopt,rwork,lrw,iwork,liw,jdum,jt)

      Return
      End      
      
      Subroutine cppdlsodaP(ti,tout,neq,xgst)
      ! at constant pressure
      ! 'ti' is the initial value of the independent variable.
      ! 'tout' is the first point where output is desired (.ne. ti).
      ! 'neq' is number of first order ode-s.
      ! 'xgst' is array of initial values, of length neq. The last value in this array is temperature. The others are mass fraction.
      ! lsoda, mf should be jt
      implicit none
      include "chemkin/ckvariables.f"
      external Fcn11rP

      integer itol,mf,itask,istate,iopt,lrw,liw,jt
      double precision ti,tout,dt,atol,rtol,jdum

      parameter (lrw=50000,liw=50000)
      double precision rwork(lrw),iwork(liw)
      integer neq
      double precision xgst(neq)
      common /lsodestore/ atol,rtol,mf,itask,istate,iopt,itol,jt

      call dlsoda (Fcn11rP,neq,xgst,ti,tout,itol,rtol,atol,
     >itask,istate,iopt,rwork,lrw,iwork,liw,jdum,jt)

      Return
      End

      Subroutine cppdlsodaPT(ti,tout,neq,xgst)
      ! at constant pressure and constant temperature
      ! 'ti' is the initial value of the independent variable.
      ! 'tout' is the first point where output is desired (.ne. ti).
      ! 'neq' is number of first order ode-s.
      ! 'xgst' is array of initial values, of length neq. The last value in this array is temperature. The others are mass fraction.
      ! lsoda, mf should be jt
      implicit none
      include "chemkin/ckvariables.f"
      external Fcn11rPT

      integer itol,mf,itask,istate,iopt,lrw,liw,jt
      double precision ti,tout,dt,atol,rtol,jdum

      parameter (lrw=50000,liw=50000)
      double precision rwork(lrw),iwork(liw)
      integer neq
      double precision xgst(neq)
      common /lsodestore/ atol,rtol,mf,itask,istate,iopt,itol,jt

      call dlsoda (Fcn11rPT,neq,xgst,ti,tout,itol,rtol,atol,
     >itask,istate,iopt,rwork,lrw,iwork,liw,jdum,jt)

      Return
      End

      Subroutine cppdlsodaST(ti,tout,neq,ct)
      ! surface reaction
      ! at constant temperature
      ! 'ti' is the initial value of the independent variable.
      ! 'tout' is the first point where output is desired (.ne. ti).
      ! 'neq' is number of first order ode-s.
      ! 'ct' is array of initial values, of length neq. The last value in this array is temperature. The others are molar concentration
      ! lsoda, mf should be jt
      implicit none
      include "chemkin/ckvariables.f"
      external Fcn11rST

      integer itol,mf,itask,istate,iopt,lrw,liw,jt
      double precision ti,tout,dt,atol,rtol,jdum

      parameter (lrw=50000,liw=50000)
      double precision rwork(lrw),iwork(liw)
      integer neq
      double precision ct(neq)
      common /lsodestore/ atol,rtol,mf,itask,istate,iopt,itol,jt

      call dlsoda (Fcn11rST,neq,ct,ti,tout,itol,rtol,atol,
     >itask,istate,iopt,rwork,lrw,iwork,liw,jdum,jt)

      Return
      End
 
      Subroutine cppdlsodaSTCC1(ti,tout,neq,ct)
	  ! hold the concentration of the first species to be constant
      ! surface reaction
      ! at constant temperature
      ! 'ti' is the initial value of the independent variable.
      ! 'tout' is the first point where output is desired (.ne. ti).
      ! 'neq' is number of first order ode-s.
      ! 'ct' is array of initial values, of length neq. The last value in this array is temperature. The others are molar concentration
      ! lsoda, mf should be jt
      implicit none
      include "chemkin/ckvariables.f"
      external Fcn11rSTCC1

      integer itol,mf,itask,istate,iopt,lrw,liw,jt
      double precision ti,tout,dt,atol,rtol,jdum

      parameter (lrw=50000,liw=50000)
      double precision rwork(lrw),iwork(liw)
      integer neq
      double precision ct(neq)
      common /lsodestore/ atol,rtol,mf,itask,istate,iopt,itol,jt

      call dlsoda (Fcn11rSTCC1,neq,ct,ti,tout,itol,rtol,atol,
     >itask,istate,iopt,rwork,lrw,iwork,liw,jdum,jt)

      Return
      End 

      Subroutine cppdlsodaSTCC2(ti,tout,neq,ct)
	  ! hold the concentration of the first species to be constant
      ! surface reaction
      ! at constant temperature
      ! 'ti' is the initial value of the independent variable.
      ! 'tout' is the first point where output is desired (.ne. ti).
      ! 'neq' is number of first order ode-s.
      ! 'ct' is array of initial values, of length neq. The last value in this array is temperature. The others are molar concentration
      ! lsoda, mf should be jt
      implicit none
      include "chemkin/ckvariables.f"
      external Fcn11rSTCC2

      integer itol,mf,itask,istate,iopt,lrw,liw,jt
      double precision ti,tout,dt,atol,rtol,jdum

      parameter (lrw=50000,liw=50000)
      double precision rwork(lrw),iwork(liw)
      integer neq
      double precision ct(neq)
      common /lsodestore/ atol,rtol,mf,itask,istate,iopt,itol,jt

      call dlsoda (Fcn11rSTCC2,neq,ct,ti,tout,itol,rtol,atol,
     >itask,istate,iopt,rwork,lrw,iwork,liw,jdum,jt)

      Return
      End 	  
	  
      Subroutine Fcn11r(neq,t,xgst,derstr)
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

      Subroutine Fcn11rVT(neq,t,xgst,derstr)
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
      
      Subroutine Fcn11rP(neq,t,xgst,derstr)
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

      Subroutine Fcn11rPT(neq,t,xgst,derstr)
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

      Subroutine Fcn11rST(neq,t,ct,derstr)
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

      Subroutine Fcn11rSTCC1(neq,t,ct,derstr)
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

      Subroutine Fcn11rSTCC2(neq,t,ct,derstr)
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
