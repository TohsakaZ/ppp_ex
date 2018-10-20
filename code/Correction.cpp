#include "Correction.h"
#include "CoordTrans.h"
#include <math.h>

Correction::Correction()
{

}

Correction::~Correction()
{

}

/*****************定义对流层延迟改正简化模型函数*******************************/
/*
   函数名：tropSimple()
   输入：E 卫星相对于测站的高度角
   输出： 对流层折射延迟改正，单位（m）
*/
double Correction::tropSimple(const double E)
{
	double ans;
	ans = 2.47*(1 / (sin(E*pi / 180) + 0.0121));
	return (ans);
}

/*****************定义对流层折射改正SAAS模型函数*******************************/
/*
   函数名：tropSaas();
   输入：X1           测站笛卡尔坐标
        X2           卫星笛卡尔坐标
		meteoData_0  测得的气象资料
   输出：对流层折射延迟
*/
double Correction::tropSaas(const Cartesian &X1,const Cartesian &X2)
{
	Topopolar T2p;
    Geodetic pos;
    meteodata meteodata_0; //标准条件下的气象数据
    const double temp0 = 15.0;

    double e,z,trph,trpw;

    pos = CoordTrans::Cart2Geod(X1);
    T2p = CoordTrans::Topo2Topop(CoordTrans::Cart2Topo(X1, X2));
    if(pos.H <-100 || pos.H > 1E4 ||T2p.E<=0.0) return 0.0;
    meteodata_0.height = pos.H;
    meteodata_0.pres =1013.25*pow(1.0-2.2557E-5*meteodata_0.height,5.2568);
    meteodata_0.temp = temp0-6.5E-3*meteodata_0.height+273.16;
    meteodata_0.RH = 0.7;
    
    e=6.108*meteodata_0.RH*exp((17.15*meteodata_0.temp-4684.0)/(meteodata_0.temp-38.45));
    
    /* saastamoninen model */
    z=pi/2.0-(T2p.E /180.0 *pi);
    trph=0.0022768*meteodata_0.pres/(1.0-0.00266*cos(2.0*(pos.B/180.0 *pi))-0.00028*meteodata_0.height/1E3)/cos(z);
    trpw=0.002277*(1255.0/meteodata_0.temp+0.05)*e/cos(z);
    return trph+trpw;
    
//    es = meteodata_0.RH*exp(-37.2465 + 0.213166*meteodata_0.temp -
//        0.000256908*meteodata_0.temp*meteodata_0.temp);
//
//
//    aa = 1.16 - 0.15E-3*meteodata_0.height + 0.716E-8*pow(meteodata_0.height, 2);
//    dlta_E = 16 / meteodata_0.temp *(meteodata_0.pres + 4810 / meteodata_0.temp*es)*(1 / tan(T2p.E*pi/180.0));
//    Eq = T2p.E + dlta_E;
//
//
//
//    ans = 0.002277 / sin(T2p.E*pi / 180.0)*(meteodata_0.pres + (1255/meteodata_0.temp +0.05)*es
//        - aa / (tan(T2p.E*pi / 180.0)*tan(T2p.E*pi / 180.0)));
}

/*****************定义电离层延迟Klobuchar模型改正函数*************************/
/*
   函数名：ionKlobuchar()
   输入：a[4],b[4]   Klobuchar模型系数，从星历文件中的头文件header获取
        gpstime     观测时刻的GPS时，单位（S）
		X1          测站的笛卡尔坐标
		X2          卫星的笛卡尔坐标
   输出：用Klobuchar改正模型计算得出电离层延迟改正，单位（m）
*/
double Correction::ionKlobuchar(const double a[4],const double b[4],const int &gpstime,
   const	Cartesian &X1,const Cartesian &X2)
{
	double Psi, Phi_i, lambda_i, Phi_m;  //单位是semi—circls
	Topopolar T2p;
	Geodetic g1;
	double t, F, PER, AMP, x;
	
	T2p = CoordTrans::Topo2Topop(CoordTrans::Cart2Topo(X1, X2));
    if (T2p.E<=0.0) return 0.0;
    
	Psi = 0.0137 / (T2p.E / 180.0 + 0.11) - 0.022;
	g1 = CoordTrans::Cart2Geod(X1);

	Phi_i = g1.B / 180.0 + Psi*cos(T2p.A*pi / 180.0);
	if (Phi_i > 0.416)
	{
		Phi_i = 0.416;
	}
	else if (Phi_i < -0.416)
	{
		Phi_i = -0.416;
	}

	lambda_i = g1.L / 180.0 + Psi*sin(T2p.A*pi / 180.0)/cos(Phi_i*pi);

	t = 4.32E4*lambda_i + gpstime;
	if (t>86400)
	{
		t = t - 86400;
	}
	else if (t < 0)
	{
		t = t + 86400;
	}

	Phi_m = Phi_i + 0.064*cos((lambda_i - 1.617)*pi);
	F = 1.0 + 16.0*pow((0.53 - T2p.E / 180.0),3);
	
	//Phi_m的单位是否是semi——circles
	PER = b[0] + b[1] * (Phi_m) + b[2] * pow(Phi_m, 2)
		+ b[3] * pow(Phi_m, 3);
	if (PER < 72000) PER = 72000;

	AMP = a[0] + a[1] * (Phi_m) + a[2] * pow(Phi_m, 2)
		+ a[3] * pow(Phi_m, 3);
	if (AMP < 0) AMP = 0;

	x = 2.0 * pi*(t - 50400) / PER;

	if (fabs(x) < 1.57 / 648000 * 2 * pi)
		return c*(F*(5.0E-9 + AMP*(1 - x*x / 2 + pow(x,4) / 24)));
	else return c*(F*5.0E-9);
}

/*****************地球自转改正函数*************************/
/*
函数名：EarthRotatCorr()
输入：t       卫星信号真实传播时间
     SatPos   卫星位置的笛卡尔坐标
输出：地球自转改正后的卫星坐标
*/
Cartesian Correction::EarthRotatCorr(double t, Cartesian SatPos)
{
	return {sin(omega_e*t)*SatPos.Y+cos(omega_e*t)*SatPos.X,-sin(omega_e*t)*SatPos.X +cos(omega_e*t)* SatPos.Y ,SatPos.Z};
}
