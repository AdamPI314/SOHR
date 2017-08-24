#ifdef __GFORTRAN__


      Subroutine canteracppdlsodeV(ti,tout,neq,xgst)
      ! at constant volume
      ! 'ti' is the initial value of the independent variable.
      ! 'tout' is the first point where output is desired (.ne. ti).
      ! 'neq' is number of first order ode-s.
      ! 'xgst' is array of initial values, of length neq. The last value in this array is temperature. The others are mass fraction.
      implicit none
      include "chemkin/ckvariables.f"
      external CanteraUpdate

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

      call dlsode (CanteraUpdate,neq,xgst,ti,tout,itol,rtol,atol,
     >itask,istate,iopt,rwork,lrw,iwork,liw,rjac,mf)

      Return
      End


      Subroutine canteracppdlsodaV(ti,tout,neq,xgst)
      ! at constant volume
      ! 'ti' is the initial value of the independent variable.
      ! 'tout' is the first point where output is desired (.ne. ti).
      ! 'neq' is number of first order ode-s.
      ! 'xgst' is array of initial values, of length neq. The last value in this array is temperature. The others are mass fraction.
      ! lsoda, mf should be jt
      implicit none
      include "chemkin/ckvariables.f"
      external CanteraUpdate

      integer itol,mf,itask,istate,iopt,lrw,liw,jt
      double precision ti,tout,dt,atol,rtol,jdum

      parameter (lrw=50000,liw=50000)
      double precision rwork(lrw),iwork(liw)
      integer neq
      double precision xgst(neq)
      common /lsodestore/ atol,rtol,mf,itask,istate,iopt,itol,jt

      call dlsoda (CanteraUpdate,neq,xgst,ti,tout,itol,rtol,atol,
     >itask,istate,iopt,rwork,lrw,iwork,liw,jdum,jt)

      Return
      End

#endif