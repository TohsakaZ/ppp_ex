/*****************************************************************************
 Beta1.0
 修改时间：2017.7.18
 修改内容：实现矩阵运算的函数
 
 
 
 
 *******************************************************************************/

#include "CMatrix.h"
#include <iostream>
#include <math.h>
#include "P_Struct.h"

CMatrix::CMatrix()
{
    _Row = 0;
    _Column = 0;
}

CMatrix::CMatrix(int x, int y)
{
    _Row = x;
    _Column = y;
    int i, j;
    std::vector<double> test(y,0.0);
    for (i = 0; i < x; i++)
    {
        Number.push_back(test);
    }
    
    
}

CMatrix::CMatrix(int unit)
{
    int i;
    _Row = unit;
    _Column = unit;
    std::vector<double> test(unit, 0.0);
    for (i = 0; i < unit; i++)
        Number.push_back(test);
    
    for (i = 0; i < unit; i++)
    {
        Number[i][i] = 1.0;
    }
}

double CMatrix::Num(int x, int y) const
{
    if (x >= 0 && x < _Row &&y >= 0 && y < _Column)
    {
        return Number[x][y];
    }
    else
        return 0;
}

void CMatrix::Set_number(int x,int y,double num)
{
    if (x < _Row && y < _Column)
    {
        Number[x][y] = num;
    }
}

CMatrix CMatrix::operator+(const CMatrix &other) const
{
    CMatrix sum(_Row, _Column);
    int i, j;
    for (i = 0; i < _Row; i++)
    {
        for (j = 0; j < _Column; j++)
        {
            sum.Set_number(i, j, other.Num(i,j)+Number[i][j]);
        }
    }
    return sum;
}

CMatrix CMatrix::operator-(const CMatrix &other) const
{
    CMatrix ans(_Row, _Column);
    int i, j;
    for (i = 0; i < _Row; i++)
    {
        for (j = 0; j < _Column; j++)
        {
            ans.Set_number(i, j, Number[i][j] - other.Num(i, j));
        }
    }
    return ans;
}

CMatrix CMatrix::operator*(const double &k) const
{
    CMatrix ans(_Row, _Column);
    int i, j;
    for (i = 0; i < _Row; i++)
    {
        for (j = 0; j<_Column; j++)
        {
            ans.Set_number(i, j, k*Number[i][j]);
        }
    }
    return ans;
}


CMatrix CMatrix::operator*(const CMatrix &A) const
{
    CMatrix ans(_Row,A.Column());
    int i, j, m;
    double sum;
    for (i = 0; i < _Row; i++)
    {
        for (j = 0; j < A.Column(); j++)
        {
            sum = 0;
            for (m = 0; m < _Column; m++)
            {
                sum = sum + Number[i][m] * A.Num(m,j);
            }
            ans.Set_number(i,j,sum);
        }
    }
    return ans;
}

CMatrix CMatrix::operator,(const CMatrix &other) const
{
    CMatrix ans(this->Row() + other.Row(), this->Column());
    for (int i = 0; i < this->Row();i++)
        for (int j = 0; j < ans.Column(); j++)
            ans.Set_number(i, j, this->Num(i, j));
    for (int i = this->Row(); i < ans.Row();i++)
        for (int j = 0; j < ans.Column(); j++)
            ans.Set_number(i, j, other.Num(i - this->Row(), j));
    return ans;
}

CMatrix CMatrix::operator | (const CMatrix &other) const
{
    CMatrix ans(this->Row(), this->Column() + other.Column());
    for (int i = 0; i < ans.Row();i++)
    {
        for (int j = 0; j < this->Column(); j++)
        {
            ans.Set_number(i, j, this->Num(i, j));
        }
        for (int k = this->Column(); k < ans.Column(); k++)
        {
            ans.Set_number(i, k, other.Num(i, k - this->Column()));
        }
    }
    return ans;
}
Array  CMatrix:: operator [](int x)
{
    if (x>=0 && x < Row())
    {
        _arr._num = Number[x];
        _arr._column = _arr._num.size();
        return _arr;
    }
    else
    {
        _arr._num = { 0 };
        _arr._column = 0;
        return _arr;
    }
}

double Array::operator[](int j) const
{
    if (j >= 0 && j < _column)
        return (_num[j]);
    else
        return 0;
}

CMatrix CMatrix::Transpose(const CMatrix &A)
{
    CMatrix ans(A.Column(),A.Row());
    int i, j;
    for (i = 0; i < A.Column();i++)
    {
        for (j = 0; j < A.Row(); j++)
        {
            ans.Set_number(i, j, A.Num(j, i));
        }
    }
    return ans;
}

CMatrix CMatrix::Del_Change(int x, int y, const CMatrix &A)
{
    CMatrix ans(A.Row()-1,A.Row()-1);
    int i, j;
    for (i = 0; i < x; i++)
    {
        for (j = 0; j < y; j++)
            ans.Set_number(i, j, A.Num(i, j));
        
        for (j = y + 1; j < A.Row(); j++)
            ans.Set_number(i, j - 1, A.Num(i, j));
    }
    
    for (i = x + 1; i < A.Row(); i++)
    {
        for (j = 0; j < y; j++)
            ans.Set_number(i - 1, j, A.Num(i, j));
        
        for (j = y + 1; j < A.Row(); j++)
            ans.Set_number(i - 1, j - 1, A.Num(i, j));
    }
    return ans;
    
}

double CMatrix::Det(const CMatrix &A)
{
    double ans = 0;
    int i;
    if (A.Row() == 2 && A.Column() == 2)
    {
        return(A.Num(0,0)*A.Num(1,1)-A.Num(0,1)*A.Num(1,0));
    }
    else
    {
        for (i = 0; i < A.Row(); i++)
        {
            if (A.Num(0, i)!=0)
            {
                if (i % 2 == 0)
                {
                    ans = ans + A.Num(0, i)*CMatrix::Det(CMatrix::Del_Change(0, i, A));
                }
                else ans = ans - A.Num(0, i)*CMatrix::Det(CMatrix::Del_Change(0, i, A));
            }
        }
        return ans;
    }
}

CMatrix CMatrix::SwapRow(int x, int y, const CMatrix&A)
{
    CMatrix ans(A.Row(), A.Row());
    ans = A;
    int i;
    for (i = 0; i < A.Row(); i++)
    {
        ans.Set_number(x, i, A.Num(y, i));
        ans.Set_number(y, i, A.Num(x, i));
    }
    
    return ans;
}

CMatrix CMatrix::ScaleRow(int x, double coef, const CMatrix&A)
{
    CMatrix ans(A.Row(), A.Row());
    ans = A;
    int i;
    for (i = 0; i < A.Row(); i++)
    {
        ans.Set_number(x, i, A.Num(x, i)*coef);
    }
    return ans;
}

CMatrix CMatrix::AddRow(int x, double coef, int y, const CMatrix&A)
{
    if (coef == 0) return A;
    CMatrix ans(A.Row(), A.Row());
    ans = A;
    int i;
    for (i = 0; i < A.Row(); i++)
    {
        ans.Set_number(y,i,A.Num(x,i)*coef+A.Num(y,i));
    }
    return ans;
}


CMatrix CMatrix::Inverse(const CMatrix &A)
{
    /*CMatrix ans(A.Row(),A.Row()); //伴随矩阵法 效率太低
     int i, j;
     for (i = 0; i < A.Row(); i++)
     {
     for (j = 0; j < A.Row(); j++)
     {
     if ((i + j) % 2 == 0)
     ans.Add_number(i, j, CMatrix::Det(CMatrix::Del_Change(i, j, A)));
     if ((i + j) % 2 != 0)
     ans.Add_number(i, j, -CMatrix::Det(CMatrix::Del_Change(i, j, A)));
     }
     }
     ans = CMatrix::Transpose(ans);
     ans = ans *(1.0 / CMatrix::Det(A));*/
    
    /*return ans;*/
    
    //全选主元法
    CMatrix ans(A.Row());
    CMatrix B;
    int i, j, k;
    double d;
    B = A;
    for (j = 0; j < A.Row(); j++)
    {
        //判断第j行第j列是否为零，若为零则调换至不为零的一行
        if (B.Num(j, j) == 0)
        {
            for (i = j; i < A.Row(); i++)
            {
                if (B.Num(i, j) != 0)
                {
                    B = CMatrix::SwapRow(i, j, B);
                    ans = CMatrix::SwapRow(i, j, ans);
                    break;
                }
            }
        }
        //将第j行标准化
        d = 1.0 / B.Num(j, j);
        B = CMatrix::ScaleRow(j, d, B);
        ans = CMatrix::ScaleRow(j, d, ans);

        //将第j列其他的元素置0
        for (k = 0; k < A.Row(); k++)
        {
            if (k != j)
            {
                d = -B.Num(k, j);
                B = CMatrix::AddRow(j, d, k, B);
                ans = CMatrix::AddRow(j, d, k, ans);
            }
        }
    }
    
//    int n = int(A.Row());
//    CMatrix U(n,n);
//    for  (int i=0;i<n;i++)
//    {
//        U.Set_number(i, i, A.Num(i, i));
//        for (int k=0;k<i;k++)
//        {
//            U.Set_number(i, i, U[i][i]-SQR(U[k][i]));
//        }
//
//        U.Set_number(i, i, sqrt(U[i][i]));
//        for (int j =i+1;j<n;j++)
//        {
//            U.Set_number(i, j, A.Num(i, j));
//            for (int k =0;k<i;k++)
//            {
//                U.Set_number(i,j, U[i][j]-U[k][i]*U[k][j]);
//            }
//            U.Set_number(i,j, U[i][j]/U[i][i]);
//        }
//    }
//    int temp;
//    CMatrix E(n);
//    CMatrix mu;
//    mu  = U|E;
//    for (int i=n-1;i>-1;i--)
//    {
//        temp = mu[i][i];
//        for (int k =0;k<int(mu.Column());k++)
//            mu.Set_number(i, k, mu[i][k]/temp);
//
//        for (int j= i-1;j>-1;j--)
//        {
//            temp = mu[j][i];
//            for (int l= 0;l<int(mu.Column());l++)
//                mu.Set_number(j, l, mu[j][l]-mu[i][l]*temp);
//        }
//    }
//
//    CMatrix ut(n,n);
//    for (int i=0;i<n;i++)
//        for (int j=0;j<n;j++)
//            ut.Set_number(i,j, mu[i][n+j]);
//
//
//    CMatrix ans;
//    ans = ut * CMatrix::Transpose(ut);
    
    
    //lu分解法
    return ans;
}

void CMatrix::print()
{
    
    
//    std::ofstream T;
//    T.open("DataOut.txt");
    int i, j;
    for (i = 0; i < _Row; i++)
    {
        for (j = 0; j < _Column; j++)
        {
            std::cout << Number[i][j]<<" ";
        }
        std::cout << std::endl;
    }
    std::cout <<std::endl;
//    T << std::endl;
//    T.close();
}

CMatrix::~CMatrix()
{
    
}
