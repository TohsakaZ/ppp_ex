/*****************************************************************************
Beta1.0
    修改时间：2017.7.18
    修改内容：定义矩阵运算函数




*******************************************************************************/

#ifndef CMATRIX_H
#define CMATRIX_H


#include <vector>
#include <fstream>

class Array
{
	friend class CMatrix;
public:
	double operator [](int j)const;//列标
	std::vector<double > _num;
	int _column;
};
class CMatrix
{
public:
	CMatrix();
	CMatrix(int x,int y);
	CMatrix(int unit);
	
	

	CMatrix operator +(const CMatrix &other) const;
	CMatrix operator -(const CMatrix &other) const;
	CMatrix operator *(const CMatrix &other) const;
	CMatrix operator *(const double &k) const;
	CMatrix operator ,(const CMatrix &other) const;
	CMatrix operator |(const CMatrix &other) const;
	Array   operator [](int x) ; //行标

	double Row() const{ return _Row; };
	double Column() const { return _Column; };
	double Num(int x, int y) const;
	void Set_number(int x, int y,double num);
	void print();
	
	static CMatrix SwapRow(int x, int y, const CMatrix&A);
	static CMatrix ScaleRow(int x, double coef, const CMatrix&A);
	static CMatrix AddRow(int x, double coef, int y, const CMatrix&A);
	static CMatrix Transpose(const CMatrix &A);
	static CMatrix Inverse (const CMatrix &A);
	static CMatrix Del_Change(int x, int y, const CMatrix&A);
	static double  Det(const CMatrix &A);

	~CMatrix();

private:

	int _Row, _Column;

	std::vector<std::vector<double>> Number;
   
	Array _arr;
};


#endif
