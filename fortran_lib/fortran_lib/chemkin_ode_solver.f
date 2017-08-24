      Subroutine ckcppdlsodeV(ti,tout,neq,xgst)
      ! at constant volume
      ! 'ti' is the initial value of the independent variable.
      ! 'tout' is the first point where output is desired (.ne. ti).
      ! 'neq' is number of first order ode-s.
      ! 'xgst' is array of initial values, of length neq. The last value in this array is temperature. The others are mass fraction.
      implicit none
      include "chemkin/ckvariables.f"
      external ChemkinUpdate

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

      call dlsode (ChemkinUpdate,neq,xgst,ti,tout,itol,rtol,atol,
     >itask,istate,iopt,rwork,lrw,iwork,liw,rjac,mf)

      Return
      End

      Subroutine ckcppdlsodeVT(ti,tout,neq,xgst)
      ! at constant volume and constant temperature
      ! 'ti' is the initial value of the independent variable.
      ! 'tout' is the first point where output is desired (.ne. ti).
      ! 'neq' is number of first order ode-s.
      ! 'xgst' is array of initial values, of length neq. The last value in this array is temperature. The others are mass fraction.
      implicit none
      include "chemkin/ckvariables.f"
      external ChemkinUpdateVT

      integer itol,mf,itask,istate,iopt,lrw,liw
      double precision ti,tout,dt,atol,rtol,rjac
      parameter (lrw=50000,liw=50000)
      double precision rwork(lrw),iwork(liw)
      integer neq
      double precision xgst(neq)
      common /lsodestore/ atol,rtol,mf,itask,istate,iopt,itol

      call dlsode (ChemkinUpdateVT,neq,xgst,ti,tout,itol,rtol,atol,
     >itask,istate,iopt,rwork,lrw,iwork,liw,rjac,mf)

      Return
      End

      Subroutine ckcppdlsodeP(ti,tout,neq,xgst)
      ! at constant pressure
      ! 'ti' is the initial value of the independent variable.
      ! 'tout' is the first point where output is desired (.ne. ti).
      ! 'neq' is number of first order ode-s.
      ! 'xgst' is array of initial values, of length neq. The last value in this array is temperature. The others are mass fraction.
      implicit none
      include "chemkin/ckvariables.f"
      external ChemkinUpdateP

      integer itol,mf,itask,istate,iopt,lrw,liw
      double precision ti,tout,dt,atol,rtol,rjac
      parameter (lrw=50000,liw=50000)
      double precision rwork(lrw),iwork(liw)
      integer neq
      double precision xgst(neq)
      common /lsodestore/ atol,rtol,mf,itask,istate,iopt,itol

      call dlsode (ChemkinUpdateP,neq,xgst,ti,tout,itol,rtol,atol,
     >itask,istate,iopt,rwork,lrw,iwork,liw,rjac,mf)

      Return
      End

      Subroutine ckcppdlsodePT(ti,tout,neq,xgst)
      ! at constant pressure and constant temperature
      ! 'ti' is the initial value of the independent variable.
      ! 'tout' is the first point where output is desired (.ne. ti).
      ! 'neq' is number of first order ode-s.
      ! 'xgst' is array of initial values, of length neq. The last value in this array is temperature. The others are mass fraction.
      implicit none
      include "chemkin/ckvariables.f"
      external ChemkinUpdatePT

      integer itol,mf,itask,istate,iopt,lrw,liw
      double precision ti,tout,dt,atol,rtol,rjac
      parameter (lrw=50000,liw=50000)
      double precision rwork(lrw),iwork(liw)
      integer neq
      double precision xgst(neq)
      common /lsodestore/ atol,rtol,mf,itask,istate,iopt,itol

      call dlsode (ChemkinUpdatePT,neq,xgst,ti,tout,itol,rtol,atol,
     >itask,istate,iopt,rwork,lrw,iwork,liw,rjac,mf)

      Return
      End

      Subroutine ckcppdlsodeST(ti,tout,neq,ct)
      ! surface reaction
      ! at constant temperature
      ! 'ti' is the initial value of the independent variable.
      ! 'tout' is the first point where output is desired (.ne. ti).
      ! 'neq' is number of first order ode-s.
      ! 'ct' is array of initial values, of length neq. The last value in this array is temperature. The others are molar concentration
      implicit none
      include "chemkin/ckvariables.f"
      external ChemkinUpdateST

      integer itol,mf,itask,istate,iopt,lrw,liw
      double precision ti,tout,dt,atol,rtol,rjac
      parameter (lrw=50000,liw=50000)
      double precision rwork(lrw),iwork(liw)
      integer neq
      double precision ct(neq)
      common /lsodestore/ atol,rtol,mf,itask,istate,iopt,itol

      call dlsode (ChemkinUpdateST,neq,ct,ti,tout,itol,rtol,atol,
     >itask,istate,iopt,rwork,lrw,iwork,liw,rjac,mf)

      Return
      End

      Subroutine ckcppdlsodaV(ti,tout,neq,xgst)
      ! at constant volume
      ! 'ti' is the initial value of the independent variable.
      ! 'tout' is the first point where output is desired (.ne. ti).
      ! 'neq' is number of first order ode-s.
      ! 'xgst' is array of initial values, of length neq. The last value in this array is temperature. The others are mass fraction.
      ! lsoda, mf should be jt
      implicit none
      include "chemkin/ckvariables.f"
      external ChemkinUpdate

      integer itol,mf,itask,istate,iopt,lrw,liw,jt
      double precision ti,tout,dt,atol,rtol,jdum

      parameter (lrw=50000,liw=50000)
      double precision rwork(lrw),iwork(liw)
      integer neq
      double precision xgst(neq)
      common /lsodestore/ atol,rtol,mf,itask,istate,iopt,itol,jt

      call dlsoda (ChemkinUpdate,neq,xgst,ti,tout,itol,rtol,atol,
     >itask,istate,iopt,rwork,lrw,iwork,liw,jdum,jt)

      Return
      End

      Subroutine ckcppdlsodaVT(ti,tout,neq,xgst)
      ! at constant volume and temperature
      ! 'ti' is the initial value of the independent variable.
      ! 'tout' is the first point where output is desired (.ne. ti).
      ! 'neq' is number of first order ode-s.
      ! 'xgst' is array of initial values, of length neq. The last value in this array is temperature. The others are mass fraction.
      ! lsoda, mf should be jt
      implicit none
      include "chemkin/ckvariables.f"
      external ChemkinUpdateVT

      integer itol,mf,itask,istate,iopt,lrw,liw,jt
      double precision ti,tout,dt,atol,rtol,jdum

      parameter (lrw=50000,liw=50000)
      double precision rwork(lrw),iwork(liw)
      integer neq
      double precision xgst(neq)
      common /lsodestore/ atol,rtol,mf,itask,istate,iopt,itol,jt

      call dlsoda (ChemkinUpdateVT,neq,xgst,ti,tout,itol,rtol,atol,
     >itask,istate,iopt,rwork,lrw,iwork,liw,jdum,jt)

      Return
      End

      Subroutine ckcppdlsodaP(ti,tout,neq,xgst)
      ! at constant pressure
      ! 'ti' is the initial value of the independent variable.
      ! 'tout' is the first point where output is desired (.ne. ti).
      ! 'neq' is number of first order ode-s.
      ! 'xgst' is array of initial values, of length neq. The last value in this array is temperature. The others are mass fraction.
      ! lsoda, mf should be jt
      implicit none
      include "chemkin/ckvariables.f"
      external ChemkinUpdateP

      integer itol,mf,itask,istate,iopt,lrw,liw,jt
      double precision ti,tout,dt,atol,rtol,jdum

      parameter (lrw=50000,liw=50000)
      double precision rwork(lrw),iwork(liw)
      integer neq
      double precision xgst(neq)
      common /lsodestore/ atol,rtol,mf,itask,istate,iopt,itol,jt

      call dlsoda (ChemkinUpdateP,neq,xgst,ti,tout,itol,rtol,atol,
     >itask,istate,iopt,rwork,lrw,iwork,liw,jdum,jt)

      Return
      End

      Subroutine ckcppdlsodaPT(ti,tout,neq,xgst)
      ! at constant pressure and constant temperature
      ! 'ti' is the initial value of the independent variable.
      ! 'tout' is the first point where output is desired (.ne. ti).
      ! 'neq' is number of first order ode-s.
      ! 'xgst' is array of initial values, of length neq. The last value in this array is temperature. The others are mass fraction.
      ! lsoda, mf should be jt
      implicit none
      include "chemkin/ckvariables.f"
      external ChemkinUpdatePT

      integer itol,mf,itask,istate,iopt,lrw,liw,jt
      double precision ti,tout,dt,atol,rtol,jdum

      parameter (lrw=50000,liw=50000)
      double precision rwork(lrw),iwork(liw)
      integer neq
      double precision xgst(neq)
      common /lsodestore/ atol,rtol,mf,itask,istate,iopt,itol,jt

      call dlsoda (ChemkinUpdatePT,neq,xgst,ti,tout,itol,rtol,atol,
     >itask,istate,iopt,rwork,lrw,iwork,liw,jdum,jt)

      Return
      End

      Subroutine ckcppdlsodaST(ti,tout,neq,ct)
      ! surface reaction
      ! at constant temperature
      ! 'ti' is the initial value of the independent variable.
      ! 'tout' is the first point where output is desired (.ne. ti).
      ! 'neq' is number of first order ode-s.
      ! 'ct' is array of initial values, of length neq. The last value in this array is temperature. The others are molar concentration
      ! lsoda, mf should be jt
      implicit none
      include "chemkin/ckvariables.f"
      external ChemkinUpdateST

      integer itol,mf,itask,istate,iopt,lrw,liw,jt
      double precision ti,tout,dt,atol,rtol,jdum

      parameter (lrw=50000,liw=50000)
      double precision rwork(lrw),iwork(liw)
      integer neq
      double precision ct(neq)
      common /lsodestore/ atol,rtol,mf,itask,istate,iopt,itol,jt

      call dlsoda (ChemkinUpdateST,neq,ct,ti,tout,itol,rtol,atol,
     >itask,istate,iopt,rwork,lrw,iwork,liw,jdum,jt)

      Return
      End

      Subroutine ckcppdlsodaSTCC1(ti,tout,neq,ct)
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
      external ChemkinUpdateSTCC1

      integer itol,mf,itask,istate,iopt,lrw,liw,jt
      double precision ti,tout,dt,atol,rtol,jdum

      parameter (lrw=50000,liw=50000)
      double precision rwork(lrw),iwork(liw)
      integer neq
      double precision ct(neq)
      common /lsodestore/ atol,rtol,mf,itask,istate,iopt,itol,jt

      call dlsoda (ChemkinUpdateSTCC1,neq,ct,ti,tout,itol,rtol,atol,
     >itask,istate,iopt,rwork,lrw,iwork,liw,jdum,jt)

      Return
      End

      Subroutine ckcppdlsodaSTCC2(ti,tout,neq,ct)
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
      external ChemkinUpdateSTCC2

      integer itol,mf,itask,istate,iopt,lrw,liw,jt
      double precision ti,tout,dt,atol,rtol,jdum

      parameter (lrw=50000,liw=50000)
      double precision rwork(lrw),iwork(liw)
      integer neq
      double precision ct(neq)
      common /lsodestore/ atol,rtol,mf,itask,istate,iopt,itol,jt

      call dlsoda (ChemkinUpdateSTCC2,neq,ct,ti,tout,itol,rtol,atol,
     >itask,istate,iopt,rwork,lrw,iwork,liw,jdum,jt)

      Return
      End