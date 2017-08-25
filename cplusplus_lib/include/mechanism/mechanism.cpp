/*
 * chemkincpp.cpp
 *
 *  Created on: Jul 9, 2014
 *      Author: Shirong
 */

#ifndef __CHEMKINCPP_CPP_
#define __CHEMKINCPP_CPP_

#include "mechanism.h"
#include "../../include/tools/misc/global_macros.h"

namespace mechanism {
	kinetics::kinetics()
	{
	}

	kinetics::~kinetics()
	{
	}


	void kinetics::rp(double * ru, double * ruc, double * pa)
	{
		::ckrp(::chemkin.ickwrk, ::chemkin.rckwrk, ru, ruc, pa);
	}

	void kinetics::indx(int* MM, int* KK, int* II, int* NFIT) {
		::ckindx(::chemkin.ickwrk, ::chemkin.rckwrk, MM, KK, II, NFIT);
	}

	void kinetics::xty(const double *X, double *Y) {
		::ckxty(X, ::chemkin.ickwrk, ::chemkin.rckwrk, Y);
	}

	void kinetics::ytx(const double *Y, double *X) {
		::ckytx(Y, ::chemkin.ickwrk, ::chemkin.rckwrk, X);
	}

	void kinetics::ctx(const double *C, double *X) {
		::ckctx(C, ::chemkin.ickwrk, ::chemkin.rckwrk, X);
	}

	void kinetics::cty(const double * C, double * Y)
	{
		::ckcty(C, ::chemkin.ickwrk, ::chemkin.rckwrk, Y);
	}

	void kinetics::xtcp(const double *P, const double *T, const double *X, double *C) {
		::ckxtcp(P, T, X, ::chemkin.ickwrk, ::chemkin.rckwrk, C);
	}

	void kinetics::rhoy(const double *P, const double *T, const double *Y, double *RHO) {
		::ckrhoy(P, T, Y, ::chemkin.ickwrk, ::chemkin.rckwrk, RHO);
	}

	void kinetics::py(const double *RHO, const double *T, const double *Y, double *P) {
		::ckpy(RHO, T, Y, ::chemkin.ickwrk, ::chemkin.rckwrk, P);
	}

	void kinetics::ytcr(const double *RHO, const double *T, const double *Y, double *C) {
		::ckytcr(RHO, T, Y, ::chemkin.ickwrk, ::chemkin.rckwrk, C);
	}

	void kinetics::cdyr(const double *RHO, const double *T, const double *Y, double *CDOT, double *DDOT) {
#ifdef __USE_CANTERA_
		::ctcdyr(RHO, T, Y, ::chemkin.ickwrk, ::chemkin.rckwrk, CDOT, DDOT);
#else
		::ckcdyr(RHO, T, Y, ::chemkin.ickwrk, ::chemkin.rckwrk, CDOT, DDOT);
#endif // __USE_CANTERA_
	}

	void kinetics::kfkr(const double *P, const double *T, const double *X, double *FWDK, double *REVK) {
		::ckkfkr(P, T, X, ::chemkin.ickwrk, ::chemkin.rckwrk, FWDK, REVK);
#ifdef __USE_CANTERA_
		::ctkfkr(P, T, X, ::chemkin.ickwrk, ::chemkin.rckwrk, FWDK, REVK);
#else
		::ckkfkr(P, T, X, ::chemkin.ickwrk, ::chemkin.rckwrk, FWDK, REVK);
#endif // __USE_CANTERA_
	}

	// Shirong Bai wrote a fortron subroutine to calculate the reaction rates given temperature and molar concentration
	// Applicable for reactions with rate constant independent of pressure
	// where sr stands for Shirong
	void kinetics::kfkrsr(const double *T, const double *C, double *FWDK, double *REVK) {
		::ckkfkrsr(T, C, ::chemkin.ickwrk, ::chemkin.rckwrk, FWDK, REVK);
	}

	void kinetics::kfrt(const double *P, const double *T, double *RKFT, double *RKRT) {
		::ckkfrt(P, T, ::chemkin.ickwrk, ::chemkin.rckwrk, RKFT, RKRT);
	}

	void kinetics::raex(const int *I, double *RA) {
		::ckraex(I, ::chemkin.rckwrk, RA);
	}

	void kinetics::abe(double *RA, double *RB, double *RE) {
		::ckabe(::chemkin.ickwrk, ::chemkin.rckwrk, RA, RB, RE);
	}

	void kinetics::cdc(const double *T, const double *C, double *CDOT, double *DDOT) {
		::ckcdc(T, C, ::chemkin.ickwrk, ::chemkin.rckwrk, CDOT, DDOT);
	}


} /*namespace chemkincpp_sr */

#endif /* __CHEMKINCPP_CPP_ */



