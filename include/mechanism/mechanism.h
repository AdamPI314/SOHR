/*
 * chemkincpp.h
 *
 *  Created on: Jul 9, 2014
 *      Author: Shirong
 */

#ifndef __CHEMKINCPP_H_
#define __CHEMKINCPP_H_

#include <string>
#include <iostream>
#include <vector>
#include <algorithm>
#include "../tools/misc/global_extern_vars.h"
#include "../tools/misc/misc_template.h"

namespace mechanism {
	/*
	 * reaction network reaction index and ChemKin reaction index lookup table, the reason the ChemKin index is a vector is there might exist duplicated reaction
	 * ChemKin has its own index of reaction, our reaction network redefine reaction index, here is the corresponding relation between them
	 * The former is C/C++ style index, The later is Fortran Style index, actually they are exact index in ChemKin space
	 * negative index represents backward reaction
	 */

	class kinetics {
		/*
		 * some C++ wrappers of chemkin functions which are originally in Fortran language
		 * have to initialize the system, call function chemkininitialize_() somewhere at least once
		 */
	public:
		kinetics();
		~kinetics();

#ifdef __CHEMKIN_AVAILABLE_

	public:
		//:: represent global namespace
		static const int nkk(void) { return ::ckstrt.nkk; }
		static const int nii(void) { return ::ckstrt.nii; }


	public:
		static void chemkin_init(const char* infile, const char* outfile) { ::chemkininitialize(infile, outfile); }
		static void cantera_init();
		static void rp(double *ru, double *ruc, double *pa);
		static void indx(int* MM, int* KK, int* II, int* NFIT);
		static void xty(const double *X, double *Y);
		static void ytx(const double *Y, double *X);
		static void ctx(const double *C, double *X);
		static void cty(const double *C, double *Y);
		static void xtcp(const double *P, const double *T, const double *X, double *C);

		static void rhoy(const double *P, const double *T, const double *Y, double *RHO);
		static void py(const double *RHO, const double *T, const double *Y, double *P);
		static void ytcr(const double *RHO, const double *T, const double *Y, double *C);
		static void cdyr(const double *RHO, const double *T, const double *Y, double *CDOT, double *DDOT);
		static void kfkr(const double *P, const double *T, const double *X, double *FWDK, double *REVK);
		// Shirong Bai wrote a fortron subroutine to calculate the reaction rates given temperature and molar concentration
		// Applicable for reactions with rate constant independent of pressure
		// where sr stands for Shirong
		static void kfkrsr(const double *T, const double *C, double *FWDK, double *REVK);

		static void kfrt(const double *P, const double *T, double *RKFT, double *RKRT);
		static void raex(const int *I, double *RA);
		static void abe(double *RA, double *RB, double *RE);
		static void cdc(const double *T, const double *C, double *CDOT, double *DDOT);
		//calculate temperature derivative
		static void calculate_t_dot_cv(const double *Y, const double *t, double *tdot) { ::calculatetdotv(Y, t, tdot); }
#endif // __CHEMKIN_AVAILABLE_


	};

} /*namespace chemkincpp_sr */



#endif /* __CHEMKINCPP_H_ */
