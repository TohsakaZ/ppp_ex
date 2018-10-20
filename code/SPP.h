#ifndef  SPP_H
#define  SPP_H

#include <string>
#include "NavFile.h"
#include "ObsFile.h"
#include "CMatrix.h"
#include "CoordTrans.h"
#include "Correction.h"

class SPP
{
public:
	SPP();
	~SPP();

public:
    
    static void initx(RcvInfo &info,double xi,double var,int i);
    //进行处理
    static void propos(const NavFile &navfile,const ObsFile & obsfile,const Adjust_Scheme &scheme);
    
	static SatInfo Cal_SatInfo(const NavRecord &record, GPSTime t);
    
    static void Infoinit(RcvInfo &info,const Adjust_Scheme &scheme);
    
    static void pppos(RcvInfo &info,const ObsRecord &oRec,const NavFile &nav);
    
    static int filter(const ObsRecord &oRec,const RcvInfo &info, CMatrix &x,CMatrix &P,CMatrix &H,CMatrix &v,CMatrix &R,int n,int m);
    
    static void udstate_ppp(RcvInfo &info,const ObsRecord &oRec,const NavFile &nav);
    static void udpos_ppp(RcvInfo &info);
    static void udclk_pp(RcvInfo &info);
    static void udtrop_ppp(RcvInfo &info);
    static void udbias_ppp(RcvInfo &info,const ObsRecord &oRec,const NavFile &nav);
    
    static void detslp_gf(RcvInfo &info,const ObsRecord &oRec,const NavFile &nav);
    
    static double gfmeas(const GpsObservation &obs);
    
    static int prn2num(const GpsObservation &obs);
    
    static double Prange(const GPSTime &time, const GpsObservation &obs,const NavFile &nav,const Adjust_Scheme &info,double &var);
    
    static int corrmeas(const GPSTime  &time,const GpsObservation &obs,const NavFile &nav, const Cartesian &X1,const Cartesian &X2,double azel[],const Adjust_Scheme &opt,
                        double dantr[], double dants[],double phw,double meas[],double var[] );
    
    static double varerr(double el,int type,const Adjust_Scheme &opt);
    
    static int valsat(const CMatrix &B,const CMatrix &x ,const CMatrix &l,const Cartesian &pos,const int valid);
    
    static int  Cal_RcvT(const GPSObsHdr &oHeader, const ObsRecord &oRec, const NavFile &nav,RcvInfo &info);
    
	static int  Cal_RcvInfo(const GPSObsHdr &oHeader, const ObsRecord &oRec, const NavFile &nav,RcvInfo &info);
    //计算卫星位置
    static int Satposs(const ObsRecord &oRec,const NavFile &nav,SatInfo sat[]);
    //计算接收机位置
    static int Estpos(const ObsRecord &oRec,const NavFile &nav,SatInfo sat[],Cartesian &pos,const Adjust_Scheme &scheme,RcvInfo &info,bool vsat[]);
    //计算接收机速度
    static int Estevl(const ObsRecord &oRec,const NavFile &nav,SatInfo sat[],RcvInfo &info);
    //计算伪距残差
    static int Rescode(const ObsRecord &oRec, const NavFile &nav,SatInfo sat[],Cartesian &pos,double &dtr,const Adjust_Scheme &scheme,CMatrix &B,CMatrix &L,double var[],bool vsat[]);
    //计算多普勒残差
    static int Resdop(const ObsRecord &oRec, const NavFile &nav, SatInfo sat[],RcvInfo &info,Cartesian &rec_v,double c_drift,CMatrix &H,CMatrix &v);
    //计算PPP残差
    static int Resppp(const ObsRecord &oRec, const NavFile &nav,SatInfo sat[],RcvInfo &info,CMatrix &x,CMatrix &H,CMatrix &v,CMatrix &R);
    //输出文件头
    static void outsolhead();
    //输出该状态下RcvInfo中的信息
    static void outsolstat(const int &count,const RcvInfo &info,std::string output);
    
};

#endif // ! SPP_H
