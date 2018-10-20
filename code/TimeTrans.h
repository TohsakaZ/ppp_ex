/*************************************************
Beta1.1
    修改时间：2017.07.18
    修改内容：定义结构体、常数

*************************************************/
#ifndef TimeTrans_H
#define TimeTrans_H


#include "P_Struct.h"




class TimeTrans
{
public:
	TimeTrans();
	~TimeTrans();

    static GPSTime CalendarTime2GPSTime(const CalendarTime &A);
	static JulianDay CalendarTime2JulianDay(const CalendarTime &A);//通用时转换儒略日
	static CalendarTime JulianDay2CalendarTime(const JulianDay &A);//儒略日转换通用时
	static GPSTime  JulianDay2GPSTime(const JulianDay &A);//儒略日转换GPS时
	static JulianDay GPSTime2JulianDay(const GPSTime &A);//GPS时转换儒略日
	static DayofYear CalendarTime2DayofYear(const CalendarTime &A);//通用时转换年积日(整数天)
	
};

#endif;
