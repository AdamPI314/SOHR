#ifndef _MY_CLOCK_H_
#define _MY_CLOCK_H_

class MyClock
{
private:
	double dBegin;
	double dEnd;
	double Time;
	bool   fBegin;
public:
	MyClock();
	void   begin(void);
	void   end(void);
	double GetHowLong(void);
};

#endif
