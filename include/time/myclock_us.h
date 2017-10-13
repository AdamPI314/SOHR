#ifndef __MYCLOCK_US_H_
#define __MYCLOCK_US_H_

#include <iostream>
#include <fstream>
#include <ostream>

class MyClock_us {
public:
	double dBegin;
	double dEnd;
	double Time;
	bool   fBegin;
public:
	MyClock_us() :dBegin(0), dEnd(0), Time(0), fBegin(false) {};
	~MyClock_us() {};
public:
	void   begin();
	void   end();
	double GetHowLong();
};

#endif
