#ifndef __TIMING_CPP_
#define __TIMING_CPP_
#include <stdlib.h>
#include "timing.h"

//#define __WINDOWS_
#define __LINUX_

#ifdef __LINUX_
#include <sys/time.h>
void get_walltime_(double* wcTime) {
	struct timeval tp;
	gettimeofday(&tp, NULL);
	*wcTime = (double)(tp.tv_sec + tp.tv_usec / 1000000.0);
}

#endif

#ifdef __WINDOWS_
#include <windows.h>

//https://msdn.microsoft.com/en-us/library/windows/desktop/dn553408(v=vs.85).aspx

void get_walltime_(double* wcTime) {
	LARGE_INTEGER StartingTime;
	LARGE_INTEGER Frequency;

	QueryPerformanceFrequency(&Frequency);
	QueryPerformanceCounter(&StartingTime);

	//microsecond
	//StartingTime.QuadPart *= 1000000;
	//second
	StartingTime.QuadPart *= static_cast<LONGLONG>(1.0);
	StartingTime.QuadPart /= Frequency.QuadPart;

	*wcTime = static_cast<double>(StartingTime.QuadPart);

}

#endif


void get_walltime(double* wcTime) {
	get_walltime_(wcTime);
}



#endif
