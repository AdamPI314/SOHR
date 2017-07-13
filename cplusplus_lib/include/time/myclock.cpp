#ifndef _MYCLOCK_CPP_
#define _MYCLOCK_CPP_

#include <ctime>
#include "myclock.h"

#define CLK_TCK CLOCKS_PER_SEC //There is no CLK_TCK in Linux

MyClock::MyClock()
{
	dBegin = dEnd = Time = 0.0;
	fBegin = false;
}

void MyClock::begin()
{
	if (fBegin == false)
	{
		dBegin = double(clock()) / CLK_TCK;
		fBegin = true;
	}
}

void MyClock::end()
{
	if (fBegin == true)
	{
		dEnd = double(clock()) / CLK_TCK;
		fBegin = false;
		Time = dEnd - dBegin;
	}
}

double MyClock::GetHowLong(void)
{
	if (fBegin == false)
		return Time;
	else
		return -1.0;
}

#endif
