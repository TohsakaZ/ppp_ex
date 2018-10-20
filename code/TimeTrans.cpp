/*************************************************
Beta1.1
    修改时间：2017.07.18
    修改内容：实现了不同时间的转换函数

*************************************************/
#include<math.h>
#include "TimeTrans.h"


/**************通用时间转 GPS 时***************************/
GPSTime TimeTrans::CalendarTime2GPSTime(const CalendarTime &A)
{
    const int doy[]={1,32,60,91,121,152,182,213,244,274,305,335};
    int days;
    GPSTime time;
    
    days = (A.year - 1980)*365+(A.year-1977)/4 + doy[A.month-1]+A.day-2+(A.year%4==0&&A.month>=3?1:0);
    days -= 5;
    
    time.wn = days / 7;
    time.tow.sn = (days - time.wn* 7) * 86400 + A.hour * 3600 + A.minute*60 + (int)A.second;
    time.tow.tos = A.second - (int)A.second;
    return time;
    
}

/***********通用时转儒略日*********************************/
/*        输入：通用时   A

          输出：儒略日  B  */
JulianDay TimeTrans::CalendarTime2JulianDay(const CalendarTime &A)//
{
	int y, m;//对年月进行改正
	if (A.month <= 2)
	{
		y = A.year - 1;
		m = A.month + 12;
	}
	else
	{
		y = A.year;
		m = A.month;
	}
	JulianDay B;
	TimeofDay C;
	B.day = int(365.25*y) + int(30.6001*(m + 1)) + A.day + int(1720981.5 + A.hour / 24.0);
	double n = A.hour / 24.0 + 0.5 - int(A.hour / 24.0 + 0.5) + A.minute / 24.0 / 60.0 + A.second / 24.0 / 60.0 / 60.0;
	C.sn = int(n * 24 * 60 * 60);
	C.tos = n * 24 * 60 * 60 - C.sn;
	B.tod = C;
	return B;
}



/***********儒略日转通用时*********************************/
/*        输入：儒略日   A

          输出：通用时   B */

CalendarTime TimeTrans::JulianDay2CalendarTime(const JulianDay &A)
{
	int a, b, c, d, e;
	double f, g, h;
	a = A.day + int((A.tod.sn + A.tod.tos)/ 60.0/60.0/24.0+ 0.5);
	b = a + 1537;
	c = int((b - 122.1) / 365.25);
	d = int(365.25*c);
	e = int((int(b - d) / 30.6001));
	f = ((A.tod.sn + A.tod.tos) / 60.0 / 60.0 / 24.0 + 0.5);
	g = 24 * (f - int(f));
	h = 60 * (g - int(g));
	CalendarTime B;
	B.day = b - d - int(30.6001*e);
	B.month = e - 1 - 12 * int(e / 14.0);
	B.year = c - 4715 - int((7 + B.month) / 10.0);
	B.hour = int(g);
	B.minute = int(h);
	B.second = 60 * (h - int(h));
	return B;

}





/***********儒略日转GPS时*********************************/
/*        输入：儒略日   A

          输出：GPS时   B */
GPSTime TimeTrans::JulianDay2GPSTime(const JulianDay &A)
{
	GPSTime B;
	B.wn = floor((A.day - 44244 - 2400000.5 + (A.tod.sn + A.tod.tos) / 60.0 / 60.0 / 24.0) / 7.0); 
	double m = (A.day + (A.tod.sn + A.tod.tos) / 60.0 / 60.0 / 24.0 - 2400000.5-44244 - B.wn * 7)/*- floor(A.day + (A.tod.sn + A.tod.tos) / 60.0 / 60.0 / 24.0 - 2400000.5 - B.wn * 7)*/ * 86400.0;
	B.tow.sn =floor(m);
	B.tow.tos = m - floor(m);
	return B;
}


/***********GPS时转儒略日*********************************/
/*        输入：GPS时  A

         输出：儒略日   B */
JulianDay TimeTrans::GPSTime2JulianDay(const GPSTime &A)
{
	JulianDay B;
	B.day = floor (44244 + A.wn * 7 + (A.tow.sn + A.tow.tos) / 86400.0 + 2400000.5);
	double d = 44244 + A.wn * 7 +(A.tow.sn +A.tow.tos) / 86400.0 + 2400000.5-B.day;//为0.n天
	double e = floor(d*60.0*60.0*24.0);
	B.tod.sn = e;
	B.tod.tos = 60.0*60.0*24.0*d - e;
	return B;

}




/***********通用时转换年积日*********************************/
/*        输入：通用时  A

          输出：年积日   B（整数） */
DayofYear TimeTrans::CalendarTime2DayofYear(const CalendarTime &A)
{
	DayofYear B;
	B.year = A.year;
	B.day = floor(A.month * 275 / 9.0) - floor((A.month + 9) / 12.0)*(floor((A.year - 4 * floor(A.year / 4.0) + 2) / 3) + 1) + A.day - 30;
	return B;
}


