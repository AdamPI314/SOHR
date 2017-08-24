/*
 * chemkincpp.cpp
 *
 *  Created on: Jul 9, 2014
 *      Author: Shirong
 */

#ifndef __CHEMKINCPP_CPP_
#define __CHEMKINCPP_CPP_

#include "mechanism.h"

namespace mechanism {
	kinetics::kinetics()
	{
	}

	kinetics::~kinetics()
	{
	}


	void kinetics::ckrp(double * ru, double * ruc, double * pa)
	{
		::ckrp(::chemkin.ickwrk, ::chemkin.rckwrk, ru, ruc, pa);
	}

	void kinetics::ckindx(int* MM, int* KK, int* II, int* NFIT) {
		::ckindx(::chemkin.ickwrk, ::chemkin.rckwrk, MM, KK, II, NFIT);
	}

	void kinetics::ckxty(const double *X, double *Y) {
		::ckxty(X, ::chemkin.ickwrk, ::chemkin.rckwrk, Y);
	}

	void kinetics::ckytx(const double *Y, double *X) {
		::ckytx(Y, ::chemkin.ickwrk, ::chemkin.rckwrk, X);
	}

	void kinetics::ckctx(const double *C, double *X) {
		::ckctx(C, ::chemkin.ickwrk, ::chemkin.rckwrk, X);
	}

	void kinetics::ckcty(const double * C, double * Y)
	{
		::ckcty(C, ::chemkin.ickwrk, ::chemkin.rckwrk, Y);
	}

	void kinetics::ckxtcp(const double *P, const double *T, const double *X, double *C) {
		::ckxtcp(P, T, X, ::chemkin.ickwrk, ::chemkin.rckwrk, C);
	}

	void kinetics::ckrhoy(const double *P, const double *T, const double *Y, double *RHO) {
		::ckrhoy(P, T, Y, ::chemkin.ickwrk, ::chemkin.rckwrk, RHO);
	}

	void kinetics::ckpy(const double *RHO, const double *T, const double *Y, double *P) {
		::ckpy(RHO, T, Y, ::chemkin.ickwrk, ::chemkin.rckwrk, P);
	}

	void kinetics::ckytcr(const double *RHO, const double *T, const double *Y, double *C) {
		::ckytcr(RHO, T, Y, ::chemkin.ickwrk, ::chemkin.rckwrk, C);
	}

	void kinetics::ckcdyr(const double *RHO, const double *T, const double *Y, double *CDOT, double *DDOT) {
		::ckcdyr(RHO, T, Y, ::chemkin.ickwrk, ::chemkin.rckwrk, CDOT, DDOT);
	}

	void kinetics::ckkfkr(const double *P, const double *T, const double *X, double *FWDK, double *REVK) {
		::ckkfkr(P, T, X, ::chemkin.ickwrk, ::chemkin.rckwrk, FWDK, REVK);
	}

	// Shirong Bai wrote a fortron subroutine to calculate the reaction rates given temperature and molar concentration
	// Applicable for reactions with rate constant independent of pressure
	// where sr stands for Shirong
	void kinetics::ckkfkrsr(const double *T, const double *C, double *FWDK, double *REVK) {
		::ckkfkrsr(T, C, ::chemkin.ickwrk, ::chemkin.rckwrk, FWDK, REVK);
	}

	void kinetics::ckkfrt(const double *P, const double *T, double *RKFT, double *RKRT) {
		::ckkfrt(P, T, ::chemkin.ickwrk, ::chemkin.rckwrk, RKFT, RKRT);
	}

	void kinetics::ckraex(const int *I, double *RA) {
		::ckraex(I, ::chemkin.rckwrk, RA);
	}

	void kinetics::ckabe(double *RA, double *RB, double *RE) {
		::ckabe(::chemkin.ickwrk, ::chemkin.rckwrk, RA, RB, RE);
	}

	void kinetics::ckcdc(const double *T, const double *C, double *CDOT, double *DDOT) {
		::ckcdc(T, C, ::chemkin.ickwrk, ::chemkin.rckwrk, CDOT, DDOT);
	}


} /*namespace chemkincpp_sr */

#endif /* __CHEMKINCPP_CPP_ */



