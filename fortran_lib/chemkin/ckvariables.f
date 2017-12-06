C     VARIABLES IN 'CKVARIABLES.F'
C
C        LENIWK allocates the integer working space
C        LENRWK allocates the real working space
C        LENCWK allocates the character working space
C        LENSYM is the length of a character string
C
C        LENICK is the length of the integer work array
C        LENRCK is the length of the real work array
C        LENCCK is the length of the character work array
C
C        ICKWRK(*) is the integer workspace array; dimension at least LENICK
C        RCKWRK(*) is the real workspace array; dimension at least LENRCK
C        CCKWRK(*) is the character string workspace array; dimension at least LENCCK
C
C        LINKCK is unit from which the Chemkin binary file is read
C        LOUT is the unit to which printed output is written
C
C        OTHER USEFUL INTEGER VARIABLES IN 'CKVARIABLES.F'
C           (Their meaning can be found in 'CHEMKIN III')
C           NMM , NKK , NII , MXSP, MXTB, MXTP, NCP , NCP1,
C           NCP2, NCP2T,NPAR, NLAR, NFAR, NLAN, NFAL, NREV,
C           NTHB, NRLT, NWL,  NRNU, NORD, MXORD,IcMM, IcKK, 
C           IcNC, IcPH, IcCH, IcNT, IcNU, IcNK, IcNS, IcNR, 
C           IcLT, IcRL, IcRV, IcWL, IcFL, IcFO, IcKF, IcTB, 
C           IcKN, IcKT, IcRNU,IcORD,IcKOR,NcAW, NcWT, NcTT, 
C           NcAA, NcCO, NcRV, NcLT, NcRL, NcFL, NcKT, NcWL, 
C           NcRU, NcRC, NcPA, NcKF, NcKR, NcRNU,NcKOR,NcK1, 
C           NcK2, NcK3, NcK4, NcI1, NcI2, NcI3, NcI4



      INTEGER LENIWK, LENRWK, LENCWK, LENSYM
      INTEGER LENICK, LENRCK, LENCCK
      INTEGER LINKCK, LOUT
      !PARAMETER (LENIWK=1000, LENRWK=1000, LENCWK=100,
      PARAMETER (LENIWK=16000, LENRWK=16000, LENCWK=200,
      !PARAMETER (LENIWK=300, LENRWK=600, LENCWK=20,
     >           LENSYM=16, LOUT=6, LINKCK=25)
      INTEGER ICKWRK(LENIWK)
      DOUBLE PRECISION RCKWRK(LENRWK)
      CHARACTER CCKWRK(LENCWK)*(LENSYM)

      INTEGER  NMM , NKK , NII , MXSP, MXTB, MXTP, NCP , NCP1,
     >         NCP2, NCP2T,NPAR, NLAR, NFAR, NLAN, NFAL, NREV,
     >         NTHB, NRLT, NWL,  NRNU, NORD, MXORD,IcMM, IcKK, 
     >         IcNC, IcPH, IcCH, IcNT, IcNU, IcNK, IcNS, IcNR, 
     >         IcLT, IcRL, IcRV, IcWL, IcFL, IcFO, IcKF, IcTB, 
     >         IcKN, IcKT, IcRNU,IcORD,IcKOR,NcAW, NcWT, NcTT, 
     >         NcAA, NcCO, NcRV, NcLT, NcRL, NcFL, NcKT, NcWL, 
     >         NcRU, NcRC, NcPA, NcKF, NcKR, NcRNU,NcKOR,NcK1, 
     >         NcK2, NcK3, NcK4, NcI1, NcI2, NcI3, NcI4



      COMMON /CHEMKIN/ ICKWRK, RCKWRK, CCKWRK, LENICK, LENRCK, LENCCK

      COMMON /CKSTRT/ NMM , NKK , NII , MXSP, MXTB, MXTP, NCP , NCP1,
     >                NCP2, NCP2T,NPAR, NLAR, NFAR, NLAN, NFAL, NREV,
     >                NTHB, NRLT, NWL,  NRNU, NORD, MXORD,IcMM, IcKK, 
     >                IcNC, IcPH, IcCH, IcNT, IcNU, IcNK, IcNS, IcNR, 
     >                IcLT, IcRL, IcRV, IcWL, IcFL, IcFO, IcKF, IcTB, 
     >                IcKN, IcKT, IcRNU,IcORD,IcKOR,NcAW, NcWT, NcTT, 
     >                NcAA, NcCO, NcRV, NcLT, NcRL, NcFL, NcKT, NcWL, 
     >                NcRU, NcRC, NcPA, NcKF, NcKR, NcRNU,NcKOR,NcK1, 
     >                NcK2, NcK3, NcK4, NcI1, NcI2, NcI3, NcI4

