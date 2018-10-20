/************************************************************************
Beta1.0
   修改时间：2017.7.17
   修改内容：1）完成了对导航电文文件的读取功能。
   
Beta1.1
   修改时间：11:30 2017.7.18
   修改内容：1）用stringstream替换了isstream，对相关语句做了修改，减少了stringstream的构造次数。
             2）增加了读入new_record.wn的语句。

Beta1.11
   修改时间：2017.7.18
   修改内容：1）修改了变量名。
			2）增加了const。


*************************************************************************************/

#include "NavFile.h"

NavFile::NavFile()
{
	
}

NavFile::NavFile(const string & name)
{
	ReadFile(name);
}


NavFile::~NavFile()
{
}

/**************定义读取导航电文文件函数*********/
/*
   函数名：ReadFile()
   输入：导航电文文件名
   输出：读取文件是否成功
*/

bool NavFile::ReadFile(const string & name)
{
    _filename = name;
    ifstream File;
    ofstream log;
    log.open("Point_Positioning.log", ios::app);
    File.open(_filename);
    log << "打开导航电文文件..." << _filename<<std::endl;
    int line_count=0; //line_count记录行数
    if (!File.good())
    {
        log << "打开失败！" << endl;
        return false;
    }

    /*********读取导航电文头文件***************************/
    string buffer;  //用于读取一整行的数据
    while (File.good())
    {
        line_count++;
        std::getline(File, buffer);

        if (buffer.find("RINEX VERSION / TYPE", 60) != string::npos)
        {
            _header.RINEX_version = std::atof(buffer.c_str());
            _header.RINEX_type = 'N';
        }
        else if (buffer.find("ION ALPHA ", 60) != string::npos)
        {
            _header.ion_alpha[0] = atof(num_std(buffer.substr(2,12)).c_str());
            _header.ion_alpha[1] = atof(num_std(buffer.substr(13,12)).c_str());
            _header.ion_alpha[2] = atof(num_std(buffer.substr(24,12)).c_str());
            _header.ion_alpha[3] = atof(num_std(buffer.substr(35,12)).c_str());
        }
        else if (buffer.find("ION BETA", 60) != string::npos)
        {
            _header.ion_beta[0] = atof(num_std(buffer.substr(2,12)).c_str());
            _header.ion_beta[1] = atof(num_std(buffer.substr(14,12)).c_str());
            _header.ion_beta[2] = atof(num_std(buffer.substr(26,12)).c_str());
            _header.ion_beta[3] = atof(num_std(buffer.substr(38,12)).c_str());

        }
        else if (buffer.find("DELTA-UTC: A0,A1,T,W", 60) != string::npos)
        {
            _header.utc_a0 = atof(num_std(buffer.substr(3,19)).c_str());
            _header.utc_a1 = atof(num_std(buffer.substr(22,19)).c_str());
            _header.utc_t = atoi(num_std(buffer.substr(41,9)).c_str());
            _header.utc_w = atoi(num_std(buffer.substr(50,9)).c_str());
        }
        else if (buffer.find("LEAP SECONDS" ,60) != string::npos)
        {
            _header.leap_sec = atoi(buffer.substr(0,6).c_str());
        }
        else if (buffer.find("END OF HEADER", 60) != string::npos)
        {
            _header.nheader_len = line_count;
            log << "头文件读取完毕!" << endl;
            break;
        }

    }

    line_count = 0;
    int count=0;
    int Rec_count = 0;
    NavRecord new_record;

    /******************读取导航电文数据记录*****************************/
    while (File.good())
    {
        std::getline(File, buffer);
        while (buffer.length() < 80)
        {
            buffer = buffer + " ";
        }


        if (buffer.length()<3) continue;  //判断buffer是否为空

        count = count % 8; //count记录一条记录中的行数
        count++;
        line_count++;     //line_count记录数据总行数

        if (count == 1)
        {
            string temp;

            new_record = NavRecord();   //new_record记录本条数据记录

            new_record.PRN[0] = 'G';       //读取卫星的PRN号
            if (buffer[0] == ' ') buffer[0] = '0';
            new_record.PRN[1] = buffer[0];
            new_record.PRN[2] = buffer[1];

            CalendarTime ctime = {};
            ctime.year = atoi(buffer.substr(3,2).c_str());
            ctime.month = atoi(buffer.substr(6,2).c_str());
            ctime.day = atoi(buffer.substr(9,2).c_str());
            ctime.hour = atoi(buffer.substr(12,2).c_str());
            ctime.minute = atoi(buffer.substr(15,2).c_str());
            ctime.second = atof(num_std(buffer.substr(17,5)).c_str());

            ctime.year += 2000;

            new_record.CalendarTime_0 = ctime;   //

            new_record.TOC = TimeTrans::CalendarTime2GPSTime(ctime);

            new_record.SClockBias = atof(num_std(buffer.substr(22,19)).c_str());
            new_record.SClockDri = atof(num_std(buffer.substr(41,19)).c_str());
            new_record.SClockDriV = atof(num_std(buffer.substr(60,19)).c_str());

        }
        if (count == 2)
        {
            new_record.IODE = atof(num_std(buffer.substr(3,19)).c_str());
            new_record.Crs = atof(num_std(buffer.substr(22,19)).c_str());
            new_record.deltan = atof(num_std(buffer.substr(41,19)).c_str());
            new_record.M0 = atof(num_std(buffer.substr(60,19)).c_str());
        }
        if (count == 3)
        {
            new_record.Cuc = atof(num_std(buffer.substr(3,19)).c_str());
            new_record.e = atof(num_std(buffer.substr(22,19)).c_str());
            new_record.Cus = atof(num_std(buffer.substr(41,19)).c_str());
            new_record.sqrtA = atof(num_std(buffer.substr(60,19)).c_str());
        }
        if (count == 4)
        {
            double secs;
            secs = atof(num_std(buffer.substr(3,19)).c_str());
            new_record.TOE = { 0,{ (int)secs,secs - (int)secs} };
            new_record.Cic = atof(num_std(buffer.substr(22,19)).c_str());
            new_record.omega = atof(num_std(buffer.substr(41,19)).c_str());
            new_record.Cis = atof(num_std(buffer.substr(60,19)).c_str());
        }
        if (count == 5)
        {
            new_record.i0 = atof(num_std(buffer.substr(3,19)).c_str());
            new_record.Crc = atof(num_std(buffer.substr(22,19)).c_str());
            new_record.omega1 = atof(num_std(buffer.substr(41,19)).c_str());
            new_record.dOmega= atof(num_std(buffer.substr(60,19)).c_str());
        }
        if (count == 6)
        {
            new_record.di = atof(num_std(buffer.substr(3,19)).c_str());
            new_record.l2_code = atof(num_std(buffer.substr(22,19)).c_str());
            new_record.wn= atof(num_std(buffer.substr(41,19)).c_str());
            new_record.l2_P_tag = atof(num_std(buffer.substr(60,19)).c_str());
        }
        if (count == 7)
        {
            new_record.SatAccur = atof(num_std(buffer.substr(3,19)).c_str());
            for (int i=0;i<15;i++) if (ura_eph[i]>=new_record.SatAccur)
            {
                new_record.SatAccur = i;
                break;
            }
            new_record.SatState = atof(num_std(buffer.substr(22,19)).c_str());
            new_record.TGD = atof(num_std(buffer.substr(41,19)).c_str());
            new_record.IODC = atof(num_std(buffer.substr(60,19)).c_str());
        }
        if (count == 8)
        {
            new_record.SendTime = atof(num_std(buffer.substr(3,19)).c_str());
            new_record.FitRange = atof(num_std(buffer.substr(22,19)).c_str());
            new_record.spare1 = atof(num_std(buffer.substr(41,19)).c_str());
            new_record.spare2 = atof(num_std(buffer.substr(60,19)).c_str());
            Rec_count++;
            new_record.TOC = TimeChange(new_record.TOC);
            _data[new_record.PRN].insert(make_pair(new_record.TOC,new_record));  //将该条记录存储到data中
            //_data[new_record.PRN][new_record.TOC] = new_record;
        }
    }
    _header.nrec_len = line_count;

    log << "数据记录读取完毕!   共读取" <<Rec_count<<"条记录。"<< endl;
    log << endl;
    log.close();
    File.close();
      return true;
}

/**************定义获取导航电文数据记录的函数*********/
/*
函数名：Data()
输入：所需要的导航电文数据记录的序号
输出：导航电文数据记录
*/
//const NavRecord& NavFile::Data(int num)
//{
//	if (num >= 0 && num < _data.size())
//	{
//		return NavRecord();
//	}
//
//	else return _data[num];
//}

/**************定义获取导航电文数据记录数的函数*********/
/*
函数名：getDataNum()
输入：
输出：导航电文数据记录数
*/
int NavFile::getDataNum() const
{
	return _data.size();
}


/************定义获取导航电文文件名函数**************/
/*
函数名：getName()
输入：
输出：导航电文文件的文件名
*/
string NavFile::FileName() const
{
	return _filename;
}

/***************定义获取导航电文记录的函数***************************/
/*
   函数名：GetRecord()
   输入：t    卫星发射信号的GPS时
        prn  卫星的PRN号
   输出：卫星的星历记录
*/
NavRecord NavFile::GetRecord(GPSTime t, string prn) const
{

	GPSTime temp;
	temp = NavFile::TimeChange(t);
    if (_data.find(prn) == _data.end())
        return NavRecord();
	if (_data.at(prn).find(temp) != _data.at(prn).end())
	{
		return _data.at(prn).at(temp);
	}
	else
	{
		return NavRecord();
	}
	
	//return NavRecord();
}

/******************时间转换函数*******************************/
/*
   函数名：TimeChange()
   输入：t  卫星发射信号的GPS时
   输出：对应的卫星星历参考时刻
*/
GPSTime NavFile::TimeChange(GPSTime t)
{
	int n,m;
	GPSTime ans;
	n = t.tow.sn % 7200;
	m = t.tow.sn / 7200;
	if (n >= 3600) m++;
	if (m >= 84)
	{
		ans.wn = t.wn + 1;
		ans.tow.sn = 0;
		ans.tow.tos = 0;
	}
	else
	{
		ans.wn = t.wn;
		ans.tow.sn = m * 7200;
		ans.tow.tos = 0;
	}
	return ans;
}
