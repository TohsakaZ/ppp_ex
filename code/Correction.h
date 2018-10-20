/***************************************************





******************************************************/

#ifndef CORRECTION_H
#define CORRECTION_H

#include "P_Struct.h"

class Correction
{

public:
	Correction();
	~Correction();

	/************卫星钟差改正******************************/

	/************相对论效应改正****************************/

	/************对流层延迟改正****************************/
	static double tropSimple(const double E);    //对流层折射延迟的简化模型
	static double tropSaas(const Cartesian &X1, const Cartesian &X2);      //对流层折射延迟的SAAS模型 


	/************电离层延迟改正***************************/
	static double ionKlobuchar(const double a[4],const double b[4],const int &gpstime,
		const Cartesian &X1,const Cartesian &X2);  //Klobuchar模型对对电离层延迟

	/************地球自转改正*****************************/
	static Cartesian EarthRotatCorr(double t, Cartesian SatPos);//t为卫星信号真实传播时间，（接收机钟面时+接收机钟差/c）-(卫星发射钟面时)



};




#endif;
