/*
 * ODESOLVER.cpp
 *
 *  Created on: Aug 24, 2017
 *      Author: Shirong
 */
#ifndef __ODESOLVER_CPP_
#define __ODESOLVER_CPP_

#include "../../include/odeSolver/odesolver.h"
#include "../../include/tools/misc/global_macros.h"

namespace ODE {
	solver::solver()
	{
	}

	solver::~solver()
	{
	}

	void solver::cppdlsodev(const double * ti, const double * tout, const int * neq, double * xgst)
	{
#ifdef __USE_CANTERA_
		::canteracppdlsodev(ti, tout, neq, xgst);
#else
		::ckcppdlsodev(ti, tout, neq, xgst);
#endif
	}

	void solver::cppdlsodav(const double * ti, const double * tout, const int * neq, double * xgst)
	{
#ifdef __USE_CANTERA_
		::canteracppdlsodav(ti, tout, neq, xgst);
#else
		::ckcppdlsodav(ti, tout, neq, xgst);
#endif
	}

	void solver::cppdlsodevt(const double * ti, const double * tout, const int * neq, double * xgst)
	{
		::ckcppdlsodevt(ti, tout, neq, xgst);
	}

	void solver::cppdlsodavt(const double * ti, const double * tout, const int * neq, double * xgst)
	{
		::ckcppdlsodavt(ti, tout, neq, xgst);
	}

	void solver::cppdlsodep(const double * ti, const double * tout, const int * neq, double * xgst)
	{
		::ckcppdlsodep(ti, tout, neq, xgst);
	}

	void solver::cppdlsodap(const double * ti, const double * tout, const int * neq, double * xgst)
	{
		::ckcppdlsodap(ti, tout, neq, xgst);
	}

	void solver::cppdlsodept(const double * ti, const double * tout, const int * neq, double * xgst)
	{
		::ckcppdlsodept(ti, tout, neq, xgst);
	}

	void solver::cppdlsodapt(const double * ti, const double * tout, const int * neq, double * xgst)
	{
		::ckcppdlsodapt(ti, tout, neq, xgst);
	}

	void solver::cppdlsodest(const double * ti, const double * tout, const int * neq, double * ct)
	{
		::ckcppdlsodest(ti, tout, neq, ct);
	}

	void solver::cppdlsodast(const double * ti, const double * tout, const int * neq, double * ct)
	{
		::ckcppdlsodast(ti, tout, neq, ct);
	}

	void solver::cppdlsodastcc1(const double * ti, const double * tout, const int * neq, double * ct)
	{
		::ckcppdlsodastcc1(ti, tout, neq, ct);
	}

	void solver::cppdlsodastcc2(const double * ti, const double * tout, const int * neq, double * ct)
	{
		::ckcppdlsodastcc2(ti, tout, neq, ct);
	}




} /*namespace ODE */

#endif /* __ODESOLVER_CPP_ */



