/*****************************************************************************
Beta1.0
    ??????2017.7.18
    ???????????????????????




*******************************************************************************/

#include "CoordTrans.h"
#include "CMatrix.h"
#include <math.h>
CoordTrans::CoordTrans()
{
}


CoordTrans::~CoordTrans()
{
}

/********************??????????????????????***********************************/
/*
   ????????Cart2Geod()
   ?????????????X
   ????????????
*/
Geodetic CoordTrans::Cart2Geod(const Cartesian &X)
{
	Geodetic ans;
	double B1, B2, N, W;
	ans.L = atan(X.Y / X.X);
	ans.L = ans.L / pi * 180;
	if (X.X<0 && X.Y>0) ans.L = ans.L + 180;
	if (X.X < 0 && X.Y<0) ans.L = ans.L - 180;
	B1 = atan(X.Z / sqrt(X.X*X.X + X.Y*X.Y));
	N = a / sqrt(1 - e2*sin(B1)*sin(B1));
	B2 = atan((X.Z + N*e2*sin(B1)) / sqrt(X.X*X.X + X.Y*X.Y));
	/**************??????????¦Ã??	************************/
	while (abs(B1 - B2) >(0.001 / 3600 * pi / 180))
	{
		B1 = B2;
		N = a / sqrt(1 - e2*sin(B1)*sin(B1));
		B2 = atan((X.Z + N*e2*sin(B1)) / sqrt(X.X*X.X + X.Y*X.Y));
	}
	ans.B = B2 / pi * 180;
	ans.H = X.Z / sin(B2) - N*(1 - e2);
	return (ans);
}


/**************??????????????????????***********************************/
/*
   ????????Geod2Cart()
   ???????????G
   ??????????????
*/
Cartesian CoordTrans::Geod2Cart(const Geodetic &G)
{
	Cartesian ans;
	double N, L, B, H;
	B = G.B / 180 * pi;
	L = G.L / 180 * pi;
	H = G.H;
	N = a / sqrt(1 - e2*sin(B)*sin(B));
	ans.X = (N + H)*cos(B)*cos(L);
	ans.Y = (N + H)*cos(B)*sin(L);
	ans.Z = (N*(1 - e2) + H)*sin(B);
	return ans;
}

/********************????????????????????????******************************/
/*
   ????????Cart2Topo()
   ?????¦Ï???????????X1????????????????X2
   ??????????????????
*/
Topocentric CoordTrans::Cart2Topo(const Cartesian & X1, const Cartesian & X2)
{
	Geodetic G1 = Cart2Geod(X1);
	CMatrix M1 = CMatrix(3, 3);
	CMatrix M2 = CMatrix(3, 1);
	double B = G1.B*pi / 180;
	double L = G1.L*pi / 180;

	M1.Set_number(0, 0, -sin(B)*cos(L));
	M1.Set_number(0, 1, -sin(B)*sin(L));
	M1.Set_number(0, 2, cos(B));
	M1.Set_number(1, 0, -sin(L));
	M1.Set_number(1, 1, cos(L));
	M1.Set_number(1, 2, 0);
	M1.Set_number(2, 0, cos(B)*cos(L));
	M1.Set_number(2, 1, cos(B)*sin(L));
	M1.Set_number(2, 2, sin(B));

	M2.Set_number(0, 0, X2.X - X1.X);
	M2.Set_number(1, 0, X2.Y - X1.Y);
	M2.Set_number(2, 0, X2.Z - X1.Z);

	CMatrix M3 = M1*M2;

	Topocentric Tc = { M3.Num(0, 0), M3.Num(1, 0), M3.Num(2, 0) };

	return Tc;
}

/**********************???????????????????????*****************************/
/*
   ????????Topo2Cart()
   ?????¦Ï???????????X1???????????????T2
   ???????????????????
*/
Cartesian CoordTrans::Topo2Cart(const Cartesian & X1, const Topocentric & T2)
{
	Geodetic G1 = Cart2Geod(X1);
	CMatrix M1(3, 3);
	CMatrix M2(3, 1);

	double B = G1.B*pi / 180;
	double L = G1.L*pi / 180;

	M1.Set_number(0, 0, -sin(B)*cos(L));
	M1.Set_number(0, 1, -sin(L));
	M1.Set_number(0, 2, cos(B)*cos(L));
	M1.Set_number(1, 0, -sin(B)*sin(L));
	M1.Set_number(1, 1, cos(L));
	M1.Set_number(1, 2, cos(B)*sin(L));
	M1.Set_number(2, 0, cos(B));
	M1.Set_number(2, 1, 0);
	M1.Set_number(2, 2, sin(B));

	M2.Set_number(0, 0, T2.N);
	M2.Set_number(1, 0, T2.E);
	M2.Set_number(2, 0, T2.U);

	CMatrix M3 = M1*M2;

	Cartesian cart = { M3.Num(0, 0) + X1.X, M3.Num(1, 0) + X1.Y, M3.Num(2, 0) + X1.Z };

	return cart;
}

/*************************??????????????????????????????*************************/
/*
   ????????Topo2Topop()
   ????????????????Tc
   ????????????????
*/
Topopolar CoordTrans::Topo2Topop(const Topocentric & Tc)
{

	double S = sqrt(Tc.N*Tc.N + Tc.E*Tc.E + Tc.U*Tc.U);


	double A;

	if (Tc.E >= 0 && Tc.N > 0)
		A = atan(Tc.E / Tc.N);
	else if (Tc.E >= 0 && Tc.N == 0)
		A = pi / 2;
	else if (Tc.N < 0)
		A = atan(Tc.E / Tc.N) + pi;
	else if (Tc.E < 0 && Tc.N == 0)
		A = pi * 3 / 2;
	else A = atan(Tc.E / Tc.N) + pi * 2;


	double E = asin(Tc.U / S);

	Topopolar Tp = { S, E * 180 / pi, A * 180 / pi };
	return Tp;
}

/*************************??????????????????????????????*************************/
/*
????????Topop2Topo()
???????????????Tp
?????????????????
*/
Topocentric CoordTrans::Topop2Topo(const Topopolar & Tp)
{
	double E0 = Tp.E / 180 * pi;
	double A0 = Tp.A / 180 * pi;

	double N = Tp.S*cos(E0)*cos(A0);
	double E = Tp.S*cos(E0)*sin(A0);
	double U = Tp.S*sin(E0);

	Topocentric Tc = { N, E, U };
	return Tc;
}

