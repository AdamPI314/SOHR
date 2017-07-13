      PARAMETER (LENSYM=16)

      COMMON /CKSTRT/ NMM , NKK , NII , MXSP, MXTB, MXTP, NCP , NCP1,
     1                NCP2, NCP2T,NPAR, NLAR, NFAR, NLAN, NFAL, NREV,
     2                NTHB, NRLT, NWL,  NRNU, NORD, MXORD,IcMM, IcKK, 
     3                IcNC, IcPH, IcCH, IcNT, IcNU, IcNK, IcNS, IcNR, 
     4                IcLT, IcRL, IcRV, IcWL, IcFL, IcFO, IcKF, IcTB, 
     5                IcKN, IcKT, IcRNU,IcORD,IcKOR,NcAW, NcWT, NcTT, 
     6                NcAA, NcCO, NcRV, NcLT, NcRL, NcFL, NcKT, NcWL, 
     7                NcRU, NcRC, NcPA, NcKF, NcKR, NcRNU,NcKOR,NcK1, 
     8                NcK2, NcK3, NcK4, NcI1, NcI2, NcI3, NcI4
C
C
C     VARIABLES in CHEMKIN II
C
C        LENSYM is the length of a character string
C
C        LENIWK is the length of the integer work array
C        LENRWK is the length of the real work array
C        LENCWK is the length of the character work array
C
C        ICKWRK(*) is the integer workspace array; dimension at least LENICK
C        RCKWRK(*) is the real workspace array; dimension at least LENRCK
C        CCKWRK(*) is the character string workspace array; dimension at least LENCCK
C
C        LINKCK is unit from which the Chemkin binary file is read
C        LOUT is the unit to which printed output is written
C
C        OTHER USEFUL INTEGER VARIABLES
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
C
