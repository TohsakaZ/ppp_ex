/************************************************************************
Beta1.0
   修改时间：2017.7.17
   修改内容：1）完成了各种结构体的定义。
   
Beta1.1
   修改时间：2017.7.18
   修改内容：1）增加了NavRecord中的wn变量。（用于表示周数）




*************************************************************************************/

#ifndef P_STRUCT_H
#define P_STRUCT_H

#include <string>
#include <iostream>
#include <iomanip>
using namespace std;


/******************************************定义常量*************************************************************/
const double a = 6378137;
const double e2 = 0.00669437999013;
const double pi = 3.14159265358979324;

const int maxtypeofobs = 30;
const int maxsatnume = 64;

const double c = 2.99792458E8; //m/s
const double GM = 3.986005E14; //m^3/s^2
const double omega_e = 0.000072921151467;//rad/s

const double ura_eph[]={         /* ura values (ref [3] 20.3.3.3.1.1) */
    2.4,3.4,4.85,6.85,9.65,13.65,24.0,48.0,96.0,192.0,384.0,768.0,1536.0,
    3072.0,6144.0,0.0
};

/***************************************卫星系统常量*************************/
#define MAXGPS 50    //GPS卫星系统最大数量



#define IB(s,opt)    ((s) + 3 + 1  + (opt.ion_corr_mode>=3?0:0) -1) //确定卫星数索引

#define SQR(x)      ((x)*(x))

#define EFACT_GPS   1.0                 /* error factor: GPS */

#define ERR_SAAS    0.3             /* saastamoinen model error std (m) */
#define ERR_BRDCI   0.5             /* broadcast iono model error factor */
#define ERR_CBIAS   0.3             /* code bias error std (m) */


#define VAR_POS     SQR( 30.0) /* initial variance of receiver position (m^2) */
#define VAR_CLK     SQR(100.0) /* initial variance of receiver clock (m^2) */
#define VAR_ZTD     SQR(  0.3) /* initial variance of ztd (m^2) */
#define VAR_GRA     SQR(0.001) /* initial variance of gradient (m^2) */
#define VAR_BIAS    SQR( 30.0) /* initial variance of phase-bias (m^2) */

const double lamda0 = c/(1575.42E6);  //GPS上L1载波波长（m）
const double lamda1 = c/120.0/10.23E6;  //GPS上L2载波波长（m）

/*****************************************定义数据结构体*******************************************************/

//时间结构体

/*****当前时刻*****/
struct TimeofDay
{
	int sn; //秒
	double tos; //秒的小数部分
};

/*****当前周*******/
struct TimeofWeek
{
	int sn;
	double tos;
};

/****通用时******/
struct CalendarTime
{
	int year;
	int month;
	int day;
	int hour;
	int minute;
	double second;

	friend ostream& operator<<(ostream& out, const CalendarTime &time) 
	{
		out << setw(6) << setiosflags(ios::right) << time.year;
		out << setw(2) << setiosflags(ios::left) << time.month;
		out << setw(4) << setiosflags(ios::left) << time.day;
		out << setw(3) << setiosflags(ios::left) << time.hour;
		out << setw(3) << setiosflags(ios::left) << time.minute;
		out << setw(11) << setiosflags(ios::right) <<setiosflags(ios::fixed) << setprecision(7) << setiosflags(ios::left) << time.second;
		return out;
	}
};

/****儒略日*******/
struct JulianDay
{
	int day;
	TimeofDay tod; //当前时刻
};

/****年积日*******/
struct DayofYear
{
	int  year;
	int  day;
};

/****GPS时********/
struct GPSTime
{
	int wn; //周数
	TimeofWeek tow;//时刻

	friend bool operator < (const GPSTime &t1, const GPSTime &t2)//定义了GPSTime的"<"运算，这样的话GPSTime结构就可以作为map容器的键值
	{
		return (t1.wn + (t1.tow.sn+t1.tow.tos) / 604800) < (t2.wn + (t2.tow.sn+t2.tow.tos) / 604800);
	}
    
    double operator -(const GPSTime &t2) const
    {
        double tt;
        tt  = tow.sn+tow.tos -(t2.tow.sn + t2.tow.tos);
        tt += (wn-t2.wn) *604800;
        return tt;
    }

	GPSTime operator +(const double &t2) const//还没测试
	{
		GPSTime new_t = { wn,tow };
		new_t.tow.tos += t2;
		if (new_t.tow.tos >= 1)
		{
			new_t.tow.sn += (int)new_t.tow.tos;
			new_t.tow.tos -= (int)new_t.tow.tos;
		}
		else if (new_t.tow.tos < 0)
		{
			new_t.tow.sn += (int)new_t.tow.tos-1;
			new_t.tow.tos -= (int)new_t.tow.tos - 1;
		}

		if (new_t.tow.sn >= 604800)
		{
			new_t.wn += (int)(new_t.tow.sn / 604800.0);
			new_t.tow.sn -= (int)(new_t.tow.sn / 604800.0) * 604800;
		}
		else if (new_t.tow.sn < 0)
		{
			new_t.wn += (int)(new_t.tow.sn / 604800.0)-1;
			new_t.tow.sn -= ((int)(new_t.tow.sn / 604800.0)-1) * 604800;
		}

		return new_t;
	}

};



//坐标系结构

/*****笛卡尔坐标系******/
struct Cartesian
{
	double X;
	double Y;
	double Z;
    
    Cartesian operator+(Cartesian &value)
    {
        Cartesian new_value ;
        new_value.X = X + value.X;
        new_value.Y = Y + value.Y;
        new_value.Z = Z + value.Z;
        return new_value;
    }
    
    Cartesian operator-(Cartesian &value)
    {
        Cartesian new_value ;
        new_value.X = X - value.X;
        new_value.Y = Y - value.Y;
        new_value.Z = Z - value.Z;
        return new_value;
    }

};

/*******大地坐标系******/
struct Geodetic
{
	double B;//大地经度（度）
	double L;//大地纬度（度）
	double H;//大地高（米）
};

/*******站心地平坐标系（线坐标）***/
struct Topocentric
{
	double N;//北方向（米)
	double E;//东方向（米）
	double U;//天顶方向(米)
};

/*******站心地平坐标系（极坐标）***/
struct Topopolar
{
	double S;//到站心距离（米）
	double E;//星视仰角（度）
	double A;//星视方向角（度）
};

//卫星观测数据结构体
/*******观测数据头文件*************/
struct GPSObsHdr
{
	double RINEX_version;//版本号
	char RINEX_type;
	string mark_name;//天线标志的名称（点名）
	string mark_number;//天线标志编号；
	string receiverNumber; //接收机序列号
	string receiverType;   //接收机类型
	string recerverVersion; //接收机版本号
	string antennaNumber;//天线序列号
	string antennaType;//天线类型
	Cartesian approxPos; //测站标志的近似位置（WGS-84）
	Topocentric antennaDelta; //天线中心相对于测站标志的位置
	int obsTypeNumber;//观测值类型数目
	char obsType[maxtypeofobs][2]; //具体观测类型
	double interval; //观测值的历元间隔（秒）
	int leap_sec; //自1980年1月6日以来的跳秒数
	GPSTime startTime; //数据文件第一个观测记录的时刻
	GPSTime endTime; //数据文件最后一个观测记录的时刻
	int headLineNumber; //头文件最后一行行号
};
/******每条观测数据文件的结构体*****/
struct GpsObservation
{
	char PRN[3];//卫星名
    double P[3]; //伪距观测值
    double L[3]; //载波相位观测值
    double D[3]; //多普勒观测值
	double obs[maxtypeofobs];    //观测值
	int LLI[maxtypeofobs];       //LLI(失锁标识符）
	int signal_str[maxtypeofobs]; //信号强度
};

/******整条观测记录的结构体**********/
struct ObsRecord
{
	GPSTime obstime;      //观测历元时刻GPS时
	CalendarTime obstime_c;  //观测历元时刻通用时
	int epoch_mark;       // 历元标志
	int sat_num;          // 当前历元所观测到的卫星数
	GpsObservation obsdata[maxsatnume];
};



//导航电文观测数据结构体
/*******导航电文头文件***********/
struct GPSNavHdr
{
	double RINEX_version;    //RINEX版本
	char RINEX_type;         //文件类型
	double ion_alpha[4], ion_beta[4]; //电离层参数
	double utc_a0, utc_a1;   //用于计算UTC时间的历书参数
	int utc_t, utc_w;       //用于计算UTC时间的历书参数
	int leap_sec;           //跳秒数
	int nheader_len, nrec_len;    //文件头行数，数据行数
};

/******导航电文记录*************/
struct  NavRecord
{
	char PRN[3];   //卫星的PRN号
	//第一行
	CalendarTime CalendarTime_0;   //日历时
	GPSTime TOC;       //卫星钟的参考时间
	double SClockBias;//卫星钟的偏差(s)
	double SClockDri;//卫星钟的漂移(s/s)
	double SClockDriV;//卫星钟的漂移速度(s/s^2)
	//第二行
	double IODE;//数据、星历发布时间
	double Crs;//(m)
	double deltan;//(rad/s)
	double M0;//(rad)
	//第三行
	double Cuc;//(rad)
	double e;//轨道偏心率
	double Cus;//(rad)
	double sqrtA;//(m^0.5)
	//第四行
	GPSTime TOE;//星历的参考时刻（GPS周内的秒数）
	double Cic;//(rad)
	double omega;//(rad)
	double Cis;//(rad)
	//第五行
	double i0;//(rad)
	double Crc;//(m)
	double omega1;//(rad)
	double dOmega;//(rad/s)
	//第六行
	double di;
	double l2_code;
	double wn;//GPS周数（与TOE一同表示时间）
	double l2_P_tag;
	//第七行
	double SatAccur;//卫星精度(m)
	double SatState;//卫星健康状态
	double TGD;//(sec)
	double IODC;//钟的数据龄期
	//第八行
	double SendTime;//电文发送时刻（单位为GPS周的秒，通过交接字（HOW）中的Z计数得出）
	double FitRange;
	double spare1;
	double spare2;
};

/***********卫星位置计算数据****************/
//struct SatInfo
//{
//	Cartesian SatPos;//卫星位置
//	
//	double delta_t;//钟差
//};

//对流层延迟改正
/*******对流层数据*********/
struct meteodata
{
	double temp;//温度（K）
	double pres;//压强（hpa）
	double RH;
	double height;//高度（km）
};

//平差方案
struct  Adjust_Scheme
{
    int mode;                 //positioning model : 0 SPP 1 Static_PPP 2 Dynamic_PPP
    
    double cut_off_angle;     // 截止高度角
	int trop_corr_mode;       // 对流层改正选项
	int ion_corr_mode;        // 电离层改正选项 0 不改正 1 简化改正模型 2 SAAS模型 3 ZTD模型 4 ZTD+Grad 模型
    int niter;                // 滤波迭代次数
    
    double thresslip;         // 无几何距离组合周跳探测阈值
};

//卫星信息和接收机信息结构体
/*******卫星信息********/
struct  SatInfo
{
	string prn;  //卫星的PRN号
	GPSTime t;   //所在的时刻
	Cartesian pos;  //该时刻卫星的位置
    Cartesian velocity; //
	double delta_t;  //卫星的改正后钟差（？）
    double drift;
    double var;
    bool flag;
};

/*******单点定位得到的接收机信息*******/
struct ssat_t
{
    int slip[2] ;    //周跳标识 1 有周跳
    double azel[2]; //表示方位角/仰角 {az,el} (rad)
    
    double gf;    //无组合距离
    double phw;
    bool vs;    //valid of satllite 
};

struct RcvInfo
{
	GPSTime rcvtime; //接收时刻
	Cartesian pos;   //接收机位置的笛卡尔坐标系
    Cartesian velocity; //接收机速度笛卡尔坐标系
    Topocentric diff; //计算得到的接收机位置与观测文件给出的近似位置的差
	double delta_t;  //接收机钟差
    double clock_drift;
	int validSatNum; //有效卫星数
	double sigma0;
	double  GDOP;
	double  PDOP;
	double  TDOP;
	double  HDOP;
    double  VDOP;
    int  nx;
    double  x[50] = {0.0};       //浮点解的状态向量
    double  P[50][50] = {0.0};   // 状态向量的方差
    double tt;  //距离上个时刻的时间差
    ssat_t ssat[MAXGPS] = {{0}};
    Adjust_Scheme opt;
};



#endif
