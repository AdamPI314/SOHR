/*
 * ODESOLVER.h
 *
 *  Created on: Jul 9, 2014
 *      Author: Shirong
 */
#ifndef __ODESOLVER_H_
#define __ODESOLVER_H_

#include <string>
#include <iostream>
#include <vector>
#include <algorithm>
#include "../tools/misc/global_extern_vars.h"
#include "../tools/misc/misc_template.h"

namespace ODE {
	class solver {
	public:
		solver();
		~solver();

	public:
		static void cppdlsodev(const double *ti, const double *tout, const int *neq, double *xgst);
		static void cppdlsodav(const double *ti, const double *tout, const int *neq, double *xgst);
		static void cppdlsodevt(const double *ti, const double *tout, const int *neq, double *xgst);
		static void cppdlsodavt(const double *ti, const double *tout, const int *neq, double *xgst);
		static void cppdlsodep(const double *ti, const double *tout, const int *neq, double *xgst);
		static void cppdlsodap(const double *ti, const double *tout, const int *neq, double *xgst);
		static void cppdlsodept(const double *ti, const double *tout, const int *neq, double *xgst);
		static void cppdlsodapt(const double *ti, const double *tout, const int *neq, double *xgst);
		static void cppdlsodest(const double *ti, const double *tout, const int *neq, double *ct);
		static void cppdlsodast(const double *ti, const double *tout, const int *neq, double *ct);
		static void cppdlsodastcc1(const double *ti, const double *tout, const int *neq, double *ct);
		static void cppdlsodastcc2(const double *ti, const double *tout, const int *neq, double *ct);

	};

} /*namespace ODE */



#endif /* __ODESOLVER_H_ */
