//增加计算接收机速度的动能
//        修改时间2018.04.03

#include "SPP.h"

SPP::SPP()
{
}


SPP::~SPP()
{
}

void SPP::initx(RcvInfo &info, double xi, double var, int i)
{
    info.x[i] = xi;
    for (int j=0;j<info.nx;j++)
    {
        if (i == j)
        {
            info.P[i][j] = var;
        }
        else
        {
            info.P[i][j] = info.P[j][i] = 0.0;
        }
    }
}

void SPP::propos(const NavFile &navfile,const ObsFile &obsfile,const Adjust_Scheme &scheme)
{
    string output = "output_SPP_1.txt";
    ofstream os(output);
    RcvInfo info;
    GPSTime time;
    
    Infoinit(info, scheme);
    info.rcvtime = obsfile.Data(0).obstime;
    
    for (int i = 0; i< obsfile.getDataNum(); i++)
    {
        time = info.rcvtime;
        info.tt = obsfile.Data(i).obstime - time;
        std::cout <<" processing ::" << obsfile.Data(i).obstime_c <<std::endl;
        
        if (!Cal_RcvInfo(obsfile.Header(), obsfile.Data(i), navfile, info)) continue;
        
        //outsolstat(i,info,output);
        
        os <<obsfile.Data(i).obstime_c;
        os << setw(5) << setiosflags(ios::right) << info.validSatNum;
        os << setw(15) << setiosflags(ios::right) << fixed << setprecision(3) << info.pos.X;
        os << setw(15) << setiosflags(ios::right) << fixed << setprecision(3) << info.pos.Y;
        os << setw(15) << setiosflags(ios::right) << fixed << setprecision(3) << info.pos.Z;
        //    os << setw(15) << setiosflags(ios::right) << fixed << setprecision(8) << info.delta_t;
        os << setw(15) << setiosflags(ios::right) << fixed << setprecision(8) << info.velocity.X;
        os << setw(15) << setiosflags(ios::right) << fixed << setprecision(8) << info.velocity.Y;
        os << setw(15) << setiosflags(ios::right) << fixed << setprecision(8) << info.velocity.Z;
        //    os << setw(15) << setiosflags(ios::right) << fixed << setprecision(5) << info.sigma0;
        //    os << setw(15) << setiosflags(ios::right) << fixed << setprecision(5) << info.GDOP;
        //    os << setw(15) << setiosflags(ios::right) << fixed << setprecision(5) << info.PDOP;
        //    os << setw(15) << setiosflags(ios::right) << fixed << setprecision(5) << info.TDOP;
        //    os << setw(15) << setiosflags(ios::right) << fixed << setprecision(5) << info.HDOP;
        //    os << setw(15) << setiosflags(ios::right) << fixed << setprecision(5) << info.VDOP;
        //    os << setw(15) << setiosflags(ios::right) << fixed << setprecision(5) << info.diff.N;
        //    os << setw(15) << setiosflags(ios::right) << fixed << setprecision(5) << info.diff.E;
        //    os << setw(15) << setiosflags(ios::right) << fixed << setprecision(5) << info.diff.U;
        os << endl;
        //

    }
    os.close();
}

void SPP::outsolstat(const int &count, const RcvInfo &info, std::string output)
{
    ofstream os(output,ios::app);
    if (!count)
        os = ofstream(output,ios::out);
    
    
//    ofstream os(output);
//    os <<obsfile.Data(i).obstime_c;
    os << setw(5) << setiosflags(ios::right) << info.validSatNum;
    os << setw(15) << setiosflags(ios::right) << fixed << setprecision(3) << info.pos.X;
    os << setw(15) << setiosflags(ios::right) << fixed << setprecision(3) << info.pos.Y;
    os << setw(15) << setiosflags(ios::right) << fixed << setprecision(3) << info.pos.Z;
//    os << setw(15) << setiosflags(ios::right) << fixed << setprecision(8) << info.delta_t;
    os << setw(15) << setiosflags(ios::right) << fixed << setprecision(8) << info.velocity.X;
    os << setw(15) << setiosflags(ios::right) << fixed << setprecision(8) << info.velocity.Y;
    os << setw(15) << setiosflags(ios::right) << fixed << setprecision(8) << info.velocity.Z;
//    os << setw(15) << setiosflags(ios::right) << fixed << setprecision(5) << info.sigma0;
//    os << setw(15) << setiosflags(ios::right) << fixed << setprecision(5) << info.GDOP;
//    os << setw(15) << setiosflags(ios::right) << fixed << setprecision(5) << info.PDOP;
//    os << setw(15) << setiosflags(ios::right) << fixed << setprecision(5) << info.TDOP;
//    os << setw(15) << setiosflags(ios::right) << fixed << setprecision(5) << info.HDOP;
//    os << setw(15) << setiosflags(ios::right) << fixed << setprecision(5) << info.VDOP;
//    os << setw(15) << setiosflags(ios::right) << fixed << setprecision(5) << info.diff.N;
//    os << setw(15) << setiosflags(ios::right) << fixed << setprecision(5) << info.diff.E;
//    os << setw(15) << setiosflags(ios::right) << fixed << setprecision(5) << info.diff.U;
    os << endl;
//
   os.close();
}


void SPP::Infoinit(RcvInfo &info, const Adjust_Scheme &scheme)
{
    info.pos = {0};
    info.velocity = {0};
    info.nx = 50;
    info.tt = 0.0;
    
    info.opt = scheme;
    
}

void SPP::pppos(RcvInfo &info, const ObsRecord &oRec, const NavFile &nav)
{
    CMatrix xp,Pp,H,v,R;
    Adjust_Scheme opt = info.opt;
    SatInfo sat[oRec.sat_num];
    int iter =1 ,nv ; //滤波迭代次数
    
    udstate_ppp(info, oRec, nav);
    
    Satposs(oRec, nav, sat);
    
    //分配矩阵
    xp = CMatrix(info.nx,1);  Pp = CMatrix(info.nx,info.nx);
    for (int i = 0;i<info.nx;i++) xp.Set_number(i, 0, info.x[i]);
    nv = oRec.sat_num * 2; //观测值个数
    v =  CMatrix(nv,1); H = CMatrix(nv,info.nx); R =CMatrix(nv,nv);
    
    for (int i= 0;i < iter;i++){
        //计算残差 更新设计矩阵
        if (!(nv = Resppp(oRec, nav,sat, info, xp, H, v, R))) break;
        
        for (int j = 0;j<info.nx;j++)
            for (int k =0;k<info.nx;k++)
                Pp.Set_number(j, k, info.P[j][k]);
        
        // 矩阵计算
        if (!filter(oRec,info, xp, Pp, H, v, R, info.nx,nv)) break;

//        for (int j = 0;j<info.nx;j++) {
//            info.x[j] = xp[j][0];}
    }
    
    //更新结果
    for (int i = 0;i<info.nx;i++) {
        info.x[i] = xp[i][0];
        for (int j =0;j<info.nx;j++)
            info.P[i][j] = Pp[i][j];
    }
    
    info.pos.X = info.x[0]; info.pos.Y = info.x[1]; info.pos.Z = info.x[2];
    info.validSatNum = nv/2;
    //存储精度信息的
    
}

int SPP::
filter(const ObsRecord &oRec,const RcvInfo &info, CMatrix &x, CMatrix &P, CMatrix &H, CMatrix &v, CMatrix &R, int n, int m)
{
    int ix[n];
    int i,j,k,kk,ll;
    int prn;
    CMatrix x_,Pp_,H_,xp_,P_;
    
    //获取其中需要计算的
    for (i = k=0; i<n;i++) if (x[i][0]!= 0.0  && P[i][i] > 0.0) ix[k++] = i;
//    for (k = 0;k<4;k++)
//        ix[k] = k;
//    for (i = 0;i<oRec.sat_num; i++)
//    {
//        prn = prn2num(oRec.obsdata[i]);
//        if (prn == 0) continue;
//        ix[k] = IB(prn,info.opt);
//        k++;
//    }
    

    x_=CMatrix(k,1); xp_ =CMatrix(k,1);  H_ = CMatrix(m,k); Pp_ = CMatrix(k,k); P_ = CMatrix(k,k);
    for (i = 0;i <k;i++){
        x_.Set_number(i, 0, x[ix[i]][0]);
        for (j = 0;j<m;j++ ) H_.Set_number(j, i, H[j][ix[i]]);
        for (j = 0;j<k;j++ ) P_.Set_number(j, i, P[ix[j]][ix[i]]);
    }
    
    ofstream log;
    log.open("v.txt",ios::app);
    log << m << ":  " ;
    for (kk = 0; kk < m; kk++)
        log << v[kk][0] << "  ";
    log << endl;
    log.close();
    
    log.open("x.txt",ios::app);
    log << m << ":  " ;
    for (kk = 0; kk < k; kk++)
        log << x_[kk][0] << "  ";
    log << endl;
    log.close();
    
//    log.open("P.txt",ios::app);
//    log << k << ": " << endl;
//    for (kk = 0;kk <k; kk++){
//        for (ll = 0; ll<k;ll++)
//            log << P_[kk][ll] << "  ";
//        log << endl;
//    }
//    log << endl;
//    log.close();
//
//    log.open("H.txt",ios::app);
//    log << m << " " << k << ": " <<endl;
//    for (kk = 0;kk<m;kk++){
//        for (ll = 0; ll<k; ll++ )
//            log << H_[kk][ll] << "  ";
//        log <<endl;
//    }
//    log << endl;
//    log.close();
    
    CMatrix Q = CMatrix::Inverse(H_*P_*CMatrix::Transpose(H_) + R);
    CMatrix Kk = P_*CMatrix::Transpose(H_) *Q;
    xp_ = x_ + Kk *v;
    CMatrix I(k);
    Pp_ = (I-Kk*H_)*P_;
    //
    for (i = 0;i<k;i++)
    {
        x.Set_number((ix[i]), 0, xp_[i][0]);
        for (j = 0;j<k;j++)
            P.Set_number(ix[i], ix[j], Pp_[i][j]);
    }
    
    return 1;
}


void SPP::udstate_ppp(RcvInfo &info, const ObsRecord &oRec, const NavFile &nav)
{
    udpos_ppp(info);
    
    udclk_pp(info);
    
    if (info.opt.ion_corr_mode >= 1)
    {
        udtrop_ppp(info);
    }
    
    udbias_ppp(info, oRec, nav);
    
}
void SPP::udpos_ppp(RcvInfo &info)
{
    //如果是第一个时刻 则直接初始化为SPP结果
    if (info.x[0] == 0.0)
    {
        initx(info, info.pos.X, VAR_POS, 0);
        initx(info, info.pos.Y, VAR_POS, 1);
        initx(info, info.pos.Z, VAR_POS, 2);
    }
    
    if (info.opt.mode ==1) return ;
    
    if (info.opt.mode ==2){
        initx(info, info.pos.X, VAR_POS, 0);
        initx(info, info.pos.Y, VAR_POS, 1);
        initx(info, info.pos.Z, VAR_POS, 2);
//        initx(info, -2262920.934293, VAR_POS, 0);
//        initx(info, 5018933.289840 , VAR_POS, 1);
//        initx(info, 3209511.301931, VAR_POS, 2);
    }

}
void SPP::udclk_pp(RcvInfo &info)
{
    double dtr;
    dtr = info.delta_t;
    
    initx(info, c*dtr, VAR_CLK, 3);
}
void SPP::udtrop_ppp(RcvInfo &info)
{
    
}
void SPP::udbias_ppp(RcvInfo &info, const ObsRecord &oRec, const NavFile &nav)
{
    double meas[2],var[2],bias[MAXGPS] ={0} ,offset = 0.0;
    Geodetic pos;
    int i,j,k,prn;
    
    for ( i = 0;i < MAXGPS;i++)
        for ( j =0;j< 2;j++)
            info.ssat[i].slip[j] = 0;
    
    detslp_gf(info, oRec, nav);
    
    pos = CoordTrans::Cart2Geod(info.pos);
    
    /*  首先利用伪距方程与载波相位房产等之差确定模糊度初值  即bias
        然后统计所有没有发生周跳的卫星 bias 与前一时刻模糊度之差 即offset
        然后将取均值 加到上一时刻的模糊度估计量中
     */
    for (i = k =0; i<oRec.sat_num ; i++)
    {
        prn = prn2num(oRec.obsdata[i]);
        if (!prn) continue;
        j = IB(prn,info.opt);
        if ( ! corrmeas(oRec.obstime, oRec.obsdata[i],nav, info.pos,{0,0,0}, info.ssat[prn].azel, info.opt, NULL,
                        NULL, 0.0, meas, var))  continue;
        bias[i] = meas[0] - meas[1];
        if (info.x[j] == 0.0 || info.ssat[prn].slip[0] || info.ssat[prn].slip[1]) continue;
        offset += bias[i]-info.x[j];
        k++;
    }
    
    if (k >=2 && fabs(offset / k) > 0.005*c) {
        for (i = 0;i <MAXGPS ;i ++){
            j = IB(i+1,info.opt);
            if (info.x[j] != 0.0 ) info.x[j] += offset/k;
        }
    }
    
    for (i=0 ;i <oRec.sat_num;i++){
        prn = prn2num(oRec.obsdata[i]);
        if (!prn) continue;
        j = IB(prn,info.opt);
        //有问题
        // info.P[j][j] += SQR(0.0001*fabs(info.tt));
        //如果 没有发生周跳 就认为下一时刻的周跳等于上一时刻的周跳
        if (info.x[j] !=0.0 && !info.ssat[prn].slip[0] && !info.ssat[prn].slip[1] ) continue;
        //如果发生了周跳 就重新为周跳赋初值
        //if (bias[i] == 0) continue;
        initx(info, bias[i], VAR_BIAS, j);
    }
    
}

void SPP::detslp_gf(RcvInfo &info, const ObsRecord &oRec, const NavFile &nav)
{
    double g0,g1;
    int prn;
    for (int i = 0 ;i <oRec.sat_num;i++)
    {
        prn = prn2num(oRec.obsdata[i]);
        if (!prn) continue;
        if ((g1 = gfmeas(oRec.obsdata[i])) == 0.0 ) continue;
        
        g0 = info.ssat[prn2num(oRec.obsdata[i])].gf;
        info.ssat[prn2num(oRec.obsdata[i])].gf = g1;
        
        if (g0 !=0 && fabs(g1-g0)> info.opt.thresslip)
        {
            for (int j =0;j<2;j++)
                info.ssat[prn2num(oRec.obsdata[i])].slip[j] = 1;
        }
    }
}

double SPP::gfmeas(const GpsObservation &obs)
{
    if(obs.L[0] == 0.0 || obs.L[1] == 0.0 ) return 0.0;
    
    return (lamda0 *obs.L[0] - lamda1 *obs.L[1]);
}

int SPP::corrmeas(const GPSTime &time, const GpsObservation &obs,const NavFile &nav, const Cartesian &X1,const Cartesian &X2 ,double azel[],const Adjust_Scheme &opt,
                    double dantr[],double dants[],double phw,double meas[],double var[] )
{

    double L1,P1,ion,vari;
    meas[0] = meas[1] = var[0] = var[1] =0.0;
    
    if (obs.L[0] == 0.0 || obs.P[0] == 0.0) return 0;
    L1 = obs.L[0]*lamda0;
    P1 = obs.P[0];
    
    //这里应该有相应的电离层 改正
//    ion = Correction::ionKlobuchar(nav.Header().ion_alpha, nav.Header().ion_beta,time.tow.sn , X1, X2);
//    vari = SQR(ion*ERR_BRDCI);
    ion = 0.0;
    vari = 0.0;
    
    
    meas[0] = L1 + ion ;
    meas[1] = P1 - ion;
    var[0] = vari;
    var[1] = vari+SQR(ERR_CBIAS);
    //对卫星和接收机天线相位中心进行改正
    for (int i = 0;i < 2; i++)
    {
        if (dants) meas[i] -=dants[0];
        if (dantr) meas[i] -=dantr[0];
    }
    
    return 1;
}

double SPP::varerr(double el, int type, const Adjust_Scheme &opt)
{
    double a,b,c,f;
    a = 0.003; b= 0.003;
    if (type ==0)
        c =1.0;
    else
        c= 100;

    f = EFACT_GPS;
    if (opt.ion_corr_mode == 3)
        f *= 3.0;
    return f*f*c*c*(a*a+b*b/sin(el));
//    if (type)
//        return 100.0;
//    else
//        return 1.0;
}

int SPP::prn2num(const GpsObservation &obs)
{
    string prn;
    prn = obs.PRN;
    if (prn[0] != 'G')
        return 0;
    prn.erase(0,1);
    return (atoi(prn.c_str()));
}

double SPP::Prange(const GPSTime &time, const GpsObservation &obs, const NavFile &nav, const Adjust_Scheme &info,double &var)
{
    var = 0.0;
    double PC,P1,P2,P1_P2 = 0.0,gamma,tgd;
    P1 = obs.P[0];
    P2 = obs.P[1];
    
    gamma = SQR(lamda1) / SQR(lamda0);
    tgd = nav.GetRecord(time, obs.PRN).TGD;
    P1_P2 = c*tgd;
    //get_TGD
    if (info.ion_corr_mode == 3 && fabs(P1) > 0.0 && fabs(P2) >0.0 )
    {
        //无电离层组合改正模型
        PC = (gamma*P1-P2) / (gamma-1.0);
    }
    else //采用单频伪距模型
    {
        PC = P1 -P1_P2;
    }
    
    var = SQR(ERR_CBIAS);
    return PC;
}


SatInfo SPP::Cal_SatInfo(const NavRecord & record, GPSTime t)
{
	//NavRecord record = nav.GetRecord(t, prn);//获取某一时刻该卫星的星历记录
	SatInfo sat;
	double a = record.sqrtA*record.sqrtA;//长半轴a
	double n0 = sqrt(GM / a / a / a);//计算卫星运动的平均角速度n0
	double n = n0 + record.deltan;//平均角速度改正后的值
	double tk = (t.tow.sn + t.tow.tos) - (record.TOE.tow.sn + record.TOE.tow.tos);//
	if (tk>302400)
	{
		tk -= 604800;
	}
	else if (tk<-302400)
	{
		tk += 604800;
	}
	else
	{
		tk = tk;
	}

	double Mk = record.M0 + n*tk;//计算观测时刻平近点角M rad

	double Ek = 0;//计算观测时刻偏近点角Ek rad
	double _E , E1 ;
    double sinE;
//    while (_E>1e-12)
//    {
//        E1 = Mk + record.e*sin(Ek);
//        _E = abs(Ek - E1);
//        Ek = E1;
//    }
//
    for (E1 = Mk,sinE = Ek = 0.0; fabs(E1-Ek)>1E-12; ){
        Ek = E1; sinE = sin(Ek); E1 = Mk +record.e *sinE;
    }

	/************计算卫星位置*****************************************************/
	double f = atan2(sqrt(1 - record.e*record.e)*sin(Ek) / (1 - record.e*cos(Ek)), (cos(Ek) - record.e) / (1 - record.e*cos(Ek)));//真近点角f
	double u = record.omega1 + f;//升交角距 u rad

								 //进行摄动改正

	double _u = record.Cus*sin(2 * u) + record.Cuc*cos(2 * u);//升交角距的改正数
	double _r = record.Crs*sin(2 * u) + record.Crc*cos(2 * u);//向径的改正数
	double _i = record.Cis*sin(2 * u) + record.Cic*cos(2 * u);// 轨道倾角改正数
	u = u + _u;
	double R = a*(1 - record.e*cos(Ek)) + _r;
	double i = record.i0 + _i + record.di*tk;

	//计算卫星在轨道面坐标系中的位置
	double x = R*cos(u);
	double y = R*sin(u);

	//计算观测瞬间升交点的经度
	double L = record.omega + (record.dOmega - omega_e)*tk - omega_e*(record.TOE.tow.sn + record.TOE.tow.tos);

	//计算卫星在瞬时地球坐标系中的位置
	sat.pos.X = x*cos(L) - y*cos(i)*sin(L);
	sat.pos.Y = x*sin(L) + y*cos(i)*cos(L);
	sat.pos.Z = y*sin(i);
 
    sat.delta_t = record.SClockBias + record.SClockDri*tk + record.SClockDriV*tk*tk;
    sat.delta_t -= 2.0 * sqrt(GM)*record.sqrtA*sin(Ek)*record.e/SQR(c);
    
    //计算卫星星历误差
    if (record.SatAccur<0 || record.SatAccur>15)
        sat.var = 6144.0;
    else
        sat.var = SQR(ura_eph[int(record.SatAccur)]);
    
	////_delta_t = record.SClockBias + record.SClockDri*((t.tow.sn + t.tow.tos) - (record.TOC.tow.sn + record.TOC.tow.tos)) + record.SClockDriV*((t.tow.sn + t.tow.tos) - (record.TOC.tow.sn + record.TOC.tow.tos))*((t.tow.sn + t.tow.tos) - (record.TOC.tow.sn + record.TOC.tow.tos)) + x - record.TGD;
	//sat.delta_t = record.SClockBias + record.SClockDri*tk + record.SClockDriV*tk*tk + F - record.TGD;
	sat.prn = record.PRN;
	sat.t = t;
	return sat;
}

int SPP::Satposs( const ObsRecord &oRec, const NavFile &nav,SatInfo sat[])
{
    double tt = 1E-3;
    SatInfo sat_test;
    
    for (int i =0;i<oRec.sat_num;i++)
    {
        string prn = oRec.obsdata[i].PRN;
        if (prn[1] == ' ')
            prn[1] = '0';
        double P ;
        double t,dts;
        GPSTime time =oRec.obstime;
        NavRecord nRec;
        
        for (int j = 0;j<3;j++) if ((P=oRec.obsdata[i].P[j]) > 0.0) break;
        if (P<=0.0) {
            sat[i].flag = false;
            continue;
        }
    
        //确定初次卫星时刻
        time = time + (-P / c);
    
        //计算初次卫星钟差
        nRec  =  nav.GetRecord(time, prn);
        if (nRec.PRN[0]=='\0')
        {
            sat[i].flag = false;
            continue;
        }
        t = (time.tow.sn + time.tow.tos) - (nRec.TOE.tow.sn +nRec.TOE.tow.tos);
        if (t>302400)
        {
            t -= 604800;
        }
        else if (t<-302400)
        {
            t += 604800;
        }
        else
        {
            t = t;
        }
        
        for (int k =0;k<2;k++)
        {
            t -=nRec.SClockBias + nRec.SClockDri*t + nRec.SClockDriV*t*t;
        }
        dts = nRec.SClockBias + nRec.SClockDri*t + nRec.SClockDriV*t*t;
    
        // 确定最终发射时间并计算卫星位置及钟差
        time = time +(-dts);
        nRec = nav.GetRecord(time,prn);
        sat[i] = SPP::Cal_SatInfo(nRec, time);
        sat[i].pos = Correction::EarthRotatCorr( P/c +dts, sat[i].pos);
        
        //计算1E-3秒后卫星位置 用于计算卫星速度
        
        time  = time + tt;
        sat_test = SPP::Cal_SatInfo(nRec, time);
        sat_test.pos = Correction::EarthRotatCorr(P/c + dts -tt, sat_test.pos);
        //计算卫星速度
        sat[i].velocity.X = (sat_test.pos.X - sat[i].pos.X) /tt;
        sat[i].velocity.Y = (sat_test.pos.Y - sat[i].pos.Y) /tt;
        sat[i].velocity.Z = (sat_test.pos.Z - sat[i].pos.Z) /tt;
        sat[i].drift = (sat_test.delta_t -sat[i].delta_t) /tt;
        sat[i].flag = true;
    }
    return 1;
}

int SPP::Rescode(const ObsRecord &oRec, const NavFile &nav, SatInfo *sat, Cartesian &pos, double &dtr,const Adjust_Scheme &scheme, CMatrix &B, CMatrix &L,double var[],bool vsat[])
{
    int valid_sat_num = 0;
    
    
    for (int i = 0; i<oRec.sat_num; i++)
    {
        vsat[i] = false;
        if (!sat[i].flag)
            continue;
        double vmeas;
        double dion = 0,vion=0;
        double dtrop = 0,vtrop=0;
        double P;
        
        P = Prange(oRec.obstime,oRec.obsdata[i],nav, scheme,vmeas);
        
        if (P <=0.0) continue;
        GPSTime time = oRec.obstime + (-P / c) ;
        
        Topopolar sat_tp;
        //判断高度角
        if (pos.X != 0.0){
            sat_tp = CoordTrans::Topo2Topop(CoordTrans::Cart2Topo(pos, sat[i].pos));
            if (sat_tp.E <= scheme.cut_off_angle)  continue;
        }
        
        if (pos.X==0.0)
        {
            sat_tp.E = pi /2.0;
        }
        if(fabs(pos.X) >0.0)
        {
            dtrop = Correction::tropSaas(pos, sat[i].pos);
             vtrop = SQR(ERR_SAAS / (sin(sat_tp.E / 180.0 *pi)+0.1));
        }
        else
        {
            dtrop =0;
            vtrop = 0;
        }
 
        if (scheme.ion_corr_mode != 3 || oRec.obsdata[i].P[0] == 0.0 || oRec.obsdata[i].P[1] ==0.0 )
        {
            if(fabs(pos.X) >0.0)
            {
                dion = Correction::ionKlobuchar(nav.Header().ion_alpha, nav.Header().ion_beta, time.tow.sn, pos, sat[i].pos);
                vion = SQR(dion * ERR_BRDCI);
            }
            else
            {
                dion =0;
                vion =0;
            }
        }
        
        double dist0 = sqrt((sat[i].pos.X - pos.X)*(sat[i].pos.X - pos.X) + (sat[i].pos.Y - pos.Y)* (sat[i].pos.Y - pos.Y) + (sat[i].pos.Z - pos.Z)* (sat[i].pos.Z - pos.Z));
        //dist0 += omega_e/c *(sat[i].pos.X * pos.Y - sat[i].pos.Y*pos.X);
        
        double l = (sat[i].pos.X - pos.X) / dist0;
        double m = (sat[i].pos.Y - pos.Y) / dist0;
        double n = (sat[i].pos.Z - pos.Z) / dist0;
        

        
        CMatrix B_row(1,4);
        B_row.Set_number(0, 0, -l);
        B_row.Set_number(0, 1, -m);
        B_row.Set_number(0, 2, -n);
        B_row.Set_number(0, 3, 1);
        
        B = (B , B_row);
        
        vsat[i] = true;
        double Li = P - (dist0 + dtr - c * sat[i].delta_t + dion + dtrop);
        //double Li = dist - c*(-V_Ts) + Vion + Vtrop - dist0;
        
        CMatrix L_row(1, 1);
        L_row.Set_number(0, 0, Li);
        
        L = (L, L_row);
        
        var[valid_sat_num] = varerr((sat_tp.E/180.0 *pi), 1,scheme) + sat[i].var + vmeas +vion +vtrop;
        valid_sat_num++;
    }
    
    if (valid_sat_num < 4 )
    {
       // log << "有效卫星数小于或等于4，单点定位失败。"<<endl;
        return 0;
    }
    else return valid_sat_num;
    
}

int SPP::Resdop( const ObsRecord &oRec, const NavFile &nav, SatInfo sat[],RcvInfo &info,Cartesian &rec_v,double c_drfit,CMatrix &H, CMatrix &v)
{
    double valid= 0;
    double lam,rate,vi,rx,ry; //lam记录频率D1 波长
    Cartesian vs; //接收机与卫星的相对速度
    for (int i=0; i<oRec.sat_num; i++)
    {
        lam = 0.19029367279836487; //λ1的波长
        double D1 = oRec.obsdata[i].D[0];
        if (D1 == 0.0) continue;
        if (!sat[i].flag) continue;
        
        rx = omega_e*info.pos.Y;
        ry = -omega_e*info.pos.X;
        
        vs = sat[i].velocity - rec_v;
        vs.X += rx;
        vs.Y += ry;
        
        //计算单位向量
        double l = (sat[i].pos.X - info.pos.X) ;
        double m = (sat[i].pos.Y - info.pos.Y);
        double n = (sat[i].pos.Z - info.pos.Z) ;
        double dist = sqrt(l*l + m*m + n*n);
        l  /= dist; m /= dist; n /= dist;  //确定接收机到卫星的单位向量
        
        //构造设计矩阵
        CMatrix B_row(1,4);
        B_row.Set_number(0, 0, -l);
        B_row.Set_number(0, 1, -m);
        B_row.Set_number(0, 2, -n);
        B_row.Set_number(0, 3, -1);
        H = (H,B_row);
        
        rate = l*vs.X + m*vs.Y +n*vs.Z; //+omega_e/c *(sat[i].velocity.Y*info.pos.X+sat[i].pos.Y*rec_v.X-
                                       // sat[i].velocity.X*info.pos.Y-sat[i].pos.X*rec_v.Y); //应该这里不需要地球自转改正 因为之前给出的卫星速度就是已经经过了地球自转改正
        
        vi = -lam *D1 -(rate + c_drfit -c*sat[i].drift);
        
        //构造v矩阵
        CMatrix v_row(1, 1);
        v_row.Set_number(0, 0, vi);
        v = (v, v_row);
        valid++;
    }
    
    if (valid <4 )
        return 0;
    return 1;
}

int SPP::Resppp(const ObsRecord &oRec,const NavFile &nav, SatInfo *sat, RcvInfo &info, CMatrix &x, CMatrix &H, CMatrix &v, CMatrix &R)
{
    ofstream log1,log2;
    log1.open("sat_pos.txt",ios::app);
    log2.open("sat_t.txt",ios::app);
    Cartesian rr,dr;
    Geodetic pos;
    int prn,nv = 0;
    double r,dtrp = 0.0,vi,vart = 0.0;
    double dantr[2] ={0}, dants[2] = {0}, var[MAXGPS*2],meas[2],azel[2],varm[2]={0};
    
    
    rr.X = x[0][0] ; rr.Y = x[1][0] ; rr.Z = x[2][0];
//    x.Set_number(3, 0, x[3][0]-12);
    
    
    // 固体潮改正
//    if (opt.tiedecorr >0)
//    {
//        ...........
//    }
    
    pos =  CoordTrans::Cart2Geod(rr);
    
    for (int i = 0;i < oRec.sat_num;i++)
    {
        prn = prn2num(oRec.obsdata[i]);
        if (!prn) continue;
        //判断卫星星历是否可用
        if ((! sat[i].flag) || (!info.ssat[prn].vs)) continue;
        
        dr = sat[i].pos - rr;
        r = sqrt(SQR(dr.X)+SQR(dr.Y)+SQR(dr.Z));
        dr.X /=r; dr.Y/=r; dr.Z/=r;
        
        azel[1] = (CoordTrans::Topo2Topop(CoordTrans::Cart2Topo(rr, sat[i].pos)).E) / 180.0 *pi;
        if (azel[1] < (info.opt.cut_off_angle /180.0 *pi) ) continue;
        
        if (r <= 0.0) continue ;
        
        //对流层延迟
        dtrp = 0.0;
        vart = 0.0;
//         dtrp = Correction::tropSaas(rr, sat[i].pos);
//         vart = SQR(ERR_SAAS);
        
        
        //相位饱和改正
        //...........
        
        if (!corrmeas(oRec.obstime, oRec.obsdata[i],nav,rr,sat[i].pos, info.ssat[prn].azel, info.opt, dantr, dants,
                      info.ssat[prn].phw, meas, varm))
            continue;
        

        log1.setf(ios::fixed);
        log2.setf(ios::fixed);
        log2 << c*sat[i].delta_t << " ";
        log1 << r << " ";
        

        r += -c*sat[i].delta_t + dtrp;
        
        for (int j = 0;j <2 ;j++)
        {

            if (meas[j] == 0.0) continue;
            
            for (int k = 0;k <info.nx;k++)
                H.Set_number(nv,k, 0.0);
            
            vi = meas[j] - r;
            
            H.Set_number(nv, 0, - dr.X);  H.Set_number(nv, 1, - dr.Y);
            H.Set_number(nv, 2, - dr.Z);

            vi -= x[3][0];
            H.Set_number(nv, 3, 1.0);
            
            if (info.opt.ion_corr_mode >=3)
            {
                //设置对流层天顶方向参数
                //............
            }
            
            //
            if (j == 0) {
                vi -= x[IB(prn,info.opt)][0];
                H.Set_number(nv,IB(prn,info.opt),1.0);
            }
            
            v.Set_number(nv, 0, vi);
            
            var[nv] = varerr(azel[1], j, info.opt) + varm[j] +vart;
            
            nv++;
        }
    }
    
    log1 << endl;
    log2 << endl;
    
    log1.close();
    log2.close();

    for (int i =0;i <nv; i++) for (int j = 0;j<nv;j++)
    {
        if (i==j)
            R.Set_number(i, j, var[i]);
        else
            R.Set_number(i, j, 0.0);
    }
    return  nv;
}

int SPP::Estevl(const ObsRecord &oRec, const NavFile &nav, SatInfo sat[], RcvInfo &info)
{
    CMatrix v,x,H;
    Cartesian rec_v ={0,0,0}; //记录接收机速度
    double clock_drift = 0; //记录接收机钟差
    double iter_count = 0;
    
    //进行迭代计算
    do
    {
        v = CMatrix(0,1);  H = CMatrix(0,4); x=CMatrix(3,1);
        
        if (!Resdop(oRec, nav, sat, info, rec_v,clock_drift,H, v))
            return 0;
        
        //进行最小二乘计算
        
        x =CMatrix::Inverse(CMatrix::Transpose(H)*H)*CMatrix::Transpose(H)*v;
        
        rec_v.X += x.Num(0, 0);
        rec_v.Y += x.Num(1, 0);
        rec_v.Z += x.Num(2, 0);
        
        clock_drift += x.Num(3, 0);
        iter_count++;
        
        if (iter_count >=10)
            return 0;
        
        
    } while (abs(x.Num(0, 0)*x.Num(0, 0)+x.Num(1, 0)*x.Num(1, 0)+x.Num(2, 0)*x.Num(2, 0))>1E-8);
    
    
    //填写info信息
    info.velocity.X = rec_v.X;
    info.velocity.Y = rec_v.Y;
    info.velocity.Z = rec_v.Z;
    
    info.clock_drift = clock_drift/c;
    return 1;
}

int SPP::Estpos(const ObsRecord &oRec, const NavFile &nav, SatInfo sat[], Cartesian &pos,const Adjust_Scheme &scheme,RcvInfo &info,bool vsat[])
{
    CMatrix B,L,x,P;
    double var[MAXGPS];
    double iter_count = 0;
    double dtr =0;
    double valid = 0;
    
    //进行迭代计算
    do
    {
        B = CMatrix(0,4); L = CMatrix(0,1); x = CMatrix(0,1);
        iter_count++;
        
        //计算设计矩阵
        if (!(valid = Rescode(oRec, nav, sat, pos,dtr, scheme, B, L,var,vsat)))
            return 0;
        
        P = CMatrix(valid,valid);
        for (int k=0;k<valid;k++)
        {
            P.Set_number(k, k, 1.0/sqrt(var[k]));
           // L.Set_number(k, 0, L[k][0]/sqrt(var[k]));
        }
        //对P进行赋值....
    
        //进行最小二乘
        x = CMatrix::Inverse(CMatrix::Transpose(B)*P*B)*CMatrix::Transpose(B)*P*L;
        
    
        pos.X += x.Num(0, 0);
        pos.Y += x.Num(1, 0);
        pos.Z += x.Num(2, 0);
        
        dtr += x.Num(3, 0) ;
        
        if (iter_count > 10)
            return 0;
        
    }
    while (abs(x.Num(0, 0)*x.Num(0, 0)+x.Num(1, 0)*x.Num(1, 0)+x.Num(2, 0)*x.Num(2, 0))>1e-6);
    
    if (!valsat(B,x,L,pos,valid))
        return 0;
    
    //以下计算相应误差矩阵
    CMatrix V = B*x - L;
    CMatrix Q = CMatrix::Inverse(CMatrix::Transpose(B)*B);
    
    Geodetic pos_g = CoordTrans::Cart2Geod(pos);
    
    CMatrix K(3,3);
    K.Set_number(0, 0, -sin(pos_g.B)*cos(pos_g.L));
    K.Set_number(0, 1, -sin(pos_g.B)*sin(pos_g.L));
    K.Set_number(0, 2, cos(pos_g.B));
    K.Set_number(1, 0, -sin(pos_g.L));
    K.Set_number(1, 1, cos(pos_g.B));
    K.Set_number(1, 2, 0);
    K.Set_number(2, 0, cos(pos_g.B)*cos(pos_g.L));
    K.Set_number(2, 1, cos(pos_g.B)*sin(pos_g.L));
    K.Set_number(2, 2, sin(pos_g.B));
    
    CMatrix Q_NEU = CMatrix::Transpose(K)*Q*K;
    info.rcvtime = oRec.obstime;
    info.validSatNum = valid;
    info.pos = pos;
    info.delta_t = dtr/c;
    info.sigma0 = sqrt((CMatrix::Transpose(V)*V).Num(0, 0) / (valid - 4));
    info.GDOP = sqrt(Q[0][0]+Q[1][1]+Q[2][2]+Q[3][3]);
    info.HDOP =sqrt(Q[0][0]+Q[1][1]+Q[2][2]);
    info.PDOP = sqrt(Q[3][3]);
    info.TDOP = sqrt(Q_NEU[0][0]+Q_NEU[1][1]);
    info.VDOP = sqrt(Q_NEU[2][2]);
    return 1;

}

int SPP::valsat(const CMatrix &B, const CMatrix &x, const CMatrix &L,const Cartesian &pos,const int valid)
{
    const double chisqr[100]={      /* chi-sqr(n) (alpha=0.001) */
        10.8,13.8,16.3,18.5,20.5,22.5,24.3,26.1,27.9,29.6,
        31.3,32.9,34.5,36.1,37.7,39.3,40.8,42.3,43.8,45.3,
        46.8,48.3,49.7,51.2,52.6,54.1,55.5,56.9,58.3,59.7,
        61.1,62.5,63.9,65.2,66.6,68.0,69.3,70.7,72.1,73.4,
        74.7,76.0,77.3,78.6,80.0,81.3,82.6,84.0,85.4,86.7,
        88.0,89.3,90.6,91.9,93.3,94.7,96.0,97.4,98.7,100 ,
        101 ,102 ,103 ,104 ,105 ,107 ,108 ,109 ,110 ,112 ,
        113 ,114 ,115 ,116 ,118 ,119 ,120 ,122 ,123 ,125 ,
        126 ,127 ,128 ,129 ,131 ,132 ,133 ,134 ,135 ,137 ,
        138 ,139 ,140 ,142 ,143 ,144 ,145 ,147 ,148 ,149
    };
    
    double vv,GDOP;
    CMatrix V = B*x - L;
    CMatrix Q = CMatrix::Inverse(CMatrix::Transpose(B)*B);
    
    Geodetic pos_g = CoordTrans::Cart2Geod(pos);
    
    CMatrix K(3,3);
    K.Set_number(0, 0, -sin(pos_g.B)*cos(pos_g.L));
    K.Set_number(0, 1, -sin(pos_g.B)*sin(pos_g.L));
    K.Set_number(0, 2, cos(pos_g.B));
    K.Set_number(1, 0, -sin(pos_g.L));
    K.Set_number(1, 1, cos(pos_g.B));
    K.Set_number(1, 2, 0);
    K.Set_number(2, 0, cos(pos_g.B)*cos(pos_g.L));
    K.Set_number(2, 1, cos(pos_g.B)*sin(pos_g.L));
    K.Set_number(2, 2, sin(pos_g.B));
    
    CMatrix Q_NEU = CMatrix::Transpose(K)*Q*K;

    vv = (CMatrix::Transpose(V)*V).Num(0, 0);
    GDOP = sqrt(Q[0][0]+Q[1][1]+Q[2][2]+Q[3][3]);
    if (vv>chisqr[valid-4-1])
        return 0;
    return 1;
}

int  SPP::Cal_RcvT(const GPSObsHdr & oHeader, const ObsRecord & oRec, const NavFile & nav,RcvInfo &info)
{
    return 0;
//    Adjust_Scheme scheme = info.opt;
//    //oRec.obstime;
//    //ofstream log;
//    SatInfo sat[MAXGPS];
//    //待估计得接收机位置与钟差
//    Cartesian pos = info.pos;
//    bool vsat[MAXGPS] = {false};
//
//    //log.open("Point_Positioning.log", ios::app);
//
//
//    //log << "开始单点定位" << "观测记录时间" << oRec.obstime_c<<endl;
//
//    if (!Satposs( oRec, nav, sat))
//    {
//        //log << "计算卫星位置有误 单点定位失败！" << endl;
//
//        return 0;
//    }
//
//    if (!Estpos( oRec, nav, sat, pos , scheme, info,vsat))
//    {
//        //log << "计算接收机位置时出错 单点定位失败！" <<endl;
//        for (int i = 0;i< MAXGPS ;i++)
//        {
//            info.ssat[i].vs = false;
//            info.sat[i] = sat[i];
//        }
//
//        int prn;
//        for (int i = 0;i< MAXGPS ;i++)
//        {
//            if (vsat[i]){
//                prn = prn2num(oRec.obsdata[i]);
//                info.ssat[prn].vs = true;
//            }
//        }
//
//        return 0;
//    }
//    else
//    {
//        for (int i = 0;i< MAXGPS ;i++)
//        {
//            info.ssat[i].vs = false;
//            info.sat[i] = sat[i];
//        }
//
//        int prn;
//        for (int i = 0;i< MAXGPS ;i++)
//        {
//            if (vsat[i]){
//                prn = prn2num(oRec.obsdata[i]);
//                info.ssat[prn].vs = true;
//            }
//        }
//        Estevl( oRec, nav, sat, info);
//        return 0;
//    }
}


int  SPP::Cal_RcvInfo(const GPSObsHdr & oHeader, const ObsRecord & oRec, const NavFile & nav,RcvInfo &info)
{
    
    Adjust_Scheme scheme = info.opt;
	//oRec.obstime;
	//ofstream log;
    SatInfo sat[MAXGPS];
    //待估计得接收机位置与钟差
    Cartesian pos = info.pos;
    bool vsat[MAXGPS] = {false};
    
	//log.open("Point_Positioning.log", ios::app);
    if (oRec.obstime_c.hour == 17 && oRec.obstime_c.minute == 0)
    {
        
    }

	//log << "开始单点定位" << "观测记录时间" << oRec.obstime_c<<endl;
	
    if (!Satposs( oRec, nav, sat))
    {
        //log << "计算卫星位置有误 单点定位失败！" << endl;
        return 0;
    }
    
    if (!Estpos( oRec, nav, sat, pos , scheme, info,vsat))
    {
        //log << "计算接收机位置时出错 单点定位失败！" <<endl;
        return 0;
    }
	info.diff = CoordTrans::Cart2Topo(oHeader.approxPos,pos);
	
    
    //计算接收机速度
    if (!Estevl( oRec, nav, sat, info))
    {
        //log << "接收机速度计算失败！"<<endl;
    }
    else
        ;
       // log << "接收机速度计算成功！"<<endl;

	//log << "单点定位成功！" << endl;
	//log.close();
    
    for (int i = 0;i< MAXGPS ;i++)
        info.ssat[i].vs = false;
    
    int prn;
    for (int i = 0;i< MAXGPS ;i++)
    {
        if (vsat[i]){
            prn = prn2num(oRec.obsdata[i]);
            if (!prn) continue;
            info.ssat[prn].vs = true;
        }
    }
    if (scheme.mode >=1)
    {
        //cout << "hello world!" << endl;
        pppos(info,oRec, nav);
    }
	return 1;
}
