/************************************************************************
Beta1.0
   修改时间：2017.7.17
   修改内容：1）完成了NavFile类的定义。
   
Beta1.1
   修改时间：2017.7.18
   修改内容：1）修改了对应的cpp文件。

Beta1.11
   修改时间：2017.7.18
   修改内容：1）修改了变量名。
            2）增加了const。

Beta1.2
   修改时间：
   修改内容：1）记录改为用map的方式存储 方便后续的查找


*************************************************************************************/
#ifndef NAVFILE_H
#define NAVFILE_H
#include <string>
#include <vector>
#include <map>
#include <math.h>
#include <fstream>
#include "TimeTrans.h"
#include <iostream>

using namespace std;

class NavFile
{
public:
	NavFile();
	NavFile(const string &name); //name表示导航电文文件名
	~NavFile();

private:
	string _filename; //导航电文文件名
	map<string ,map<GPSTime,NavRecord>> _data;  //导航电文数据记录
	GPSNavHdr _header;  //导航电文头文件

public:
	bool ReadFile(const string &name); //读取导航电文文件
	//const NavRecord & Data(int num);  //取出一条记录
	int getDataNum() const;  //得到记录条数
	string FileName() const;  //获取导航电文文件名
	NavRecord GetRecord(GPSTime t, string prn) const;  //根据相应的发射时刻以及PRN号获取相应的导航电文的记录
	static GPSTime TimeChange(GPSTime t);    //时间转换函数
    static string num_std(string s)
    {
        std::size_t found = s.find('D');
        if (found == string::npos)
            return s;
        s.replace(found, 1, "E");
        return s;
    }

	const GPSNavHdr &Header() const { return _header; }
};


#endif // !NAVFILE_H
