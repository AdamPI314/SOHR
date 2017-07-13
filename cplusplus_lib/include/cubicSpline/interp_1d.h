#ifndef _INTERP_1D_H_
#define _INTERP_1D_H_
#include "nr3.h"

struct Base_interp
{
	Int n, mm, jsav, cor, dj;
	const Doub *xx, *yy;
	Base_interp(VecDoub_I &x, const Doub *y, Int m)
		: n(x.size()), mm(m), jsav(0), cor(0), xx(&x[0]), yy(y) {
		dj = MIN(1, (int)pow((Doub)n, 0.25));
	}

	Doub interp(Doub x) {
		Int jlo = cor ? hunt(x) : locate(x);
		return rawinterp(jlo, x);
	}

	Int locate(const Doub x);

	Int hunt(const Doub x);

	Doub virtual rawinterp(Int jlo, Doub x) = 0;
	virtual ~Base_interp() {};

};

struct Poly_interp : Base_interp
{
	Doub dy;
	Poly_interp(VecDoub_I &xv, VecDoub_I &yv, Int m)
		: Base_interp(xv, &yv[0], m), dy(0.) {}
	Doub rawinterp(Int jl, Doub x);
};


struct Rat_interp : Base_interp
{
	Doub dy;
	Rat_interp(VecDoub_I &xv, VecDoub_I &yv, Int m)
		: Base_interp(xv, &yv[0], m), dy(0.) {}
	Doub rawinterp(Int jl, Doub x);
};

struct Spline_interp : Base_interp
{
	VecDoub y2;

	Spline_interp(VecDoub_I &xv, VecDoub_I &yv, Doub yp1 = 1.e99, Doub ypn = 1.e99)
		: Base_interp(xv, &yv[0], 2), y2(xv.size())
	{
		sety2(&xv[0], &yv[0], yp1, ypn);
	}

	Spline_interp(VecDoub_I &xv, const Doub *yv, Doub yp1 = 1.e99, Doub ypn = 1.e99)
		: Base_interp(xv, yv, 2), y2(xv.size())
	{
		sety2(&xv[0], yv, yp1, ypn);
	}

	void sety2(const Doub *xv, const Doub *yv, Doub yp1, Doub ypn);
	Doub rawinterp(Int jl, Doub xv);
};

struct BaryRat_interp : Base_interp
{
	VecDoub w;
	Int d;
	BaryRat_interp(VecDoub_I &xv, VecDoub_I &yv, Int dd);
	Doub rawinterp(Int jl, Doub x);
	Doub interp(Doub x);
};

#endif /*_INTERP_1D_H_*/
