/*****************************************************************************
Beta1.0
    修改时间：2017.7.18
    修改内容：定义坐标转换函数




*******************************************************************************/
#ifndef COORDINATRANS_H
#define COORDINATRANS_H

#include "P_Struct.h"

class CoordTrans
{
public:
	CoordTrans();
	~CoordTrans();

public:
	/*******************************定义坐标转换函数**************************************/
	static Geodetic Cart2Geod(const Cartesian &X); //笛卡尔坐标系->大地坐标系
	static Cartesian Geod2Cart(const Geodetic &G); //大地坐标系->笛卡尔坐标系
	static Topocentric Cart2Topo(const Cartesian &X1, const Cartesian &X2); //笛卡尔坐标系->站心线坐标系
	static Cartesian Topo2Cart(const Cartesian &X1, const Topocentric &T2); //站心线坐标系->笛卡尔坐标系
	static Topopolar Topo2Topop(const Topocentric &Tc); //站心线坐标系->站心极坐标系
	static Topocentric Topop2Topo(const Topopolar &Tp); //站心极坐标系->站心线坐标系

};

#endif;