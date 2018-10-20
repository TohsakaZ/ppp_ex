/*****************************************************************************
Beta1.0
    修改时间：2017.7.18
    修改内容：1）定义了读取GPS观测数据文件函数
	         2）在ObsFile.CPP中完成了读取GPS观测数据文件函数

Beta1.11
	修改时间：2017.7.18
	修改内容：1）修改了变量名。
			 2）增加了一些函数。


*******************************************************************************/
#ifndef OBSFILE_H
#define OBSFILE_H
#include <string>
#include <vector>
#include "P_Struct.h"
#include<math.h>
#include<fstream>
#include"TimeTrans.h"

using namespace std;
class ObsFile
{
public:
	ObsFile();
	ObsFile(const string &name);//name表示观测数据文件名
	~ObsFile();

private:
	string _filename;//观测数据文件名
	vector<ObsRecord> _data;//观测数据记录
	GPSObsHdr _header;//观测数据头文件
    void convert(); //将观测数据文件传入到P、L、D中

public:
	bool ReadFile(const string &name);//读取观测数据文件
	string getName() const;//获取观测数据文件名
	ObsRecord Data(int num) const; //取出一条记录
	int getDataNum() const; //得到记录条数
    int getObstype(string type); //得到观测类型位置
	const  GPSObsHdr &Header() const;
    
};

#endif;//!OBSFILE_H
