#ifndef __GLOBAL_EXTERN_VARS_H_
#define __GLOBAL_EXTERN_VARS_H_

#include <iostream>
#include <fstream>

#include "fortran_routine_block_alias.h"

//Global variables.
//const int LENIWK=1000, LENRWK=1000, LENCWK=100, LENSYM=16;
const int LENIWK=16000, LENRWK=16000, LENCWK=200, LENSYM=16;
//const int LENIWK=300, LENRWK=600, LENCWK=20, LENSYM=16;

// C++ and Fortran mixed programming.
extern "C"
{
	// Fortran 'COMMON' data struct.
	extern struct
	{
		int nmm, nkk, nii, mxsp, mxtb, mxtp, ncp, ncp1, 
			ncp2, ncp2t, npar, nlar, nfar, nlan, nfal, nrev, 
			nthb, nrlt, nwl, nrnu, nord, mxord, icmm, ickk, 
			icnc, icph, icch, icnt, icnu, icnk, icns, icnr, 
			iclt, icrl, icrv, icwl, icfl, icfo, ickf, ictb, 
			ickn, ickt, icrnu, icord, ickor, ncaw, ncwt, nctt, 
			ncaa, ncco, ncrv, nclt, ncrl, ncfl, nckt, ncwl, 
			ncru, ncrc, ncpa, nckf, nckr, ncrnu, nckor, nck1, 
			nck2, nck3, nck4, nci1, nci2, nci3, nci4;
	} ckstrt;

	extern struct
	{
		int ickwrk[LENIWK];
		double rckwrk[LENRWK];
		char cckwrk[LENCWK][LENSYM];
		int lenick, lenrck, lencck;
	} chemkin;

	extern struct
	{
		double rhomass;
		double pressure;
	} ckstore;

	extern struct
	{
		double atol, rtol;
		int mf, itask, istate, iopt, itol;
		int jt;
		int deltaN1, deltaN2;
	} lsodestore;


	// Fortran Subroutine.
	void chemkininitialize();
	//Returns universal gas constants and the pressure of one standard atmosphere
	void ckrp(const int *ICKWRK, const double *RCKWRK, double *ru, double *ruc, double *pa);
	void cknu(const int *KDIM,const int *ICKWRK,const double *RCKWRK,int *NUKI);
	void ckitr(const int *ICKWRK,const double *RCKWRK,int *ITHB,int *IREV);
	void ckncf(const int *MDIM,const int *ICKWRK,const double *RCKWRK,int *NCF);
	void ckrhoy(const double *P,const double *T,const double *Y,const int *ICKWRK,const double *RCKWRK,double *RHO);
	void ckxty(const double *X,const int *ICKWRK,const double *RCKWRK,double *Y);
	void ckytx(const double *Y,const int *ICKWRK,const double *RCKWRK,double *X);
	void ckctx(const double *C,const int *ICKWRK,const double *RCKWRK,double *X);
	void ckcty(const double *C, const int *ICKWRK, const double *RCKWRK, double *Y);
	void ckxtcp(const double *P, const double *T, const double *X, const int *ICKWRK,const double *RCKWRK, double *C);

	void ckytcr(const double *RHO,const double *T,const double *Y,const int *ICKWRK,const double *RCKWRK,double *C);
	void cppdlsodev(const double *ti,const double *tout,const int *neq,double *xgst);
	void cppdlsodav(const double *ti,const double *tout,const int *neq,double *xgst);
	void cppdlsodevt(const double *ti, const double *tout, const int *neq, double *xgst);
	void cppdlsodavt(const double *ti, const double *tout, const int *neq, double *xgst);
	void cppdlsodep(const double *ti,const double *tout,const int *neq,double *xgst);
	void cppdlsodap(const double *ti,const double *tout,const int *neq,double *xgst);
	void cppdlsodept(const double *ti,const double *tout,const int *neq,double *xgst);
	void cppdlsodapt(const double *ti,const double *tout,const int *neq,double *xgst);
	void cppdlsodest(const double *ti,const double *tout,const int *neq,double *ct);
	void cppdlsodast(const double *ti,const double *tout,const int *neq,double *ct);
	void cppdlsodastcc1(const double *ti, const double *tout, const int *neq, double *ct);
	void cppdlsodastcc2(const double *ti, const double *tout, const int *neq, double *ct);
	void ckkfkr(const double *P,const double *T,const double *X,const int *ICKWRK,const double *RCKWRK,double *FWDK,double *REVK);
	// Shirong Bai wrote a fortron subroutine to calculate the reaction rates given temperature and molar concentration
	// Applicable for reactions with rate constant independent of pressure
	// where sr stands for Shirong
	void ckkfkrsr(const double *T,const double *C,const int *ICKWRK,const double *RCKWRK,double *FWDK,double *REVK);

	void ckpy(const double *RHO,const double *T,const double *Y,const int *ICKWRK,const double *RCKWRK,double *P);
	void ckcdyr(const double *RHO,const double *T,const double *Y,const int *ICKWRK,const double *RCKWRK,double *CDOT,double *DDOT);
	void ckkfrt(const double *P,const double *T,const int *ICKWRK,const double *RCKWRK, double *RKFT, double *RKRT);
	void ckraex(const int *I, const double *RCKWRK, double *RA);
	void ckindx(const int *ICKWRK, const double *RCKWRK, int *MM, int *KK, int *II, int *NFIT);
	void ckabe(const int *ICKWRK, const double *RCKWRK, double *RA, double *RB, double *RE);
	void ckcdc(const double *T,const double *C,const int *ICKWRK,const double *RCKWRK,double *CDOT,double *DDOT);

	void calculatetdotv(const double* Y, const double *t, double *tdot);

}

#endif