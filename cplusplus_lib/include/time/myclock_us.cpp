#ifndef __MYCLOCK_US_CPP_
#define __MYCLOCK_US_CPP_

#include "myclock_us.h"
#include "timing.h"


////////////////////////////////////////////////////////////////////////////////////
// class 
////////////////////////////////////////////////////////////////////////////////////


void MyClock_us::begin()
{
	if (fBegin == false)
	{
		get_walltime(&dBegin);
		fBegin = true;
	}
}

void MyClock_us::end()
{
	if (fBegin == true)
	{
		get_walltime(&dEnd);
		fBegin = false;
		Time = dEnd - dBegin;
	}
}

double MyClock_us::GetHowLong()
{
	if (fBegin == false)
		return Time;
	else
		return -1.0;
}


#endif
