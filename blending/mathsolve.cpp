#include "mathsolve.h"
#include<iostream>
using namespace std;
# include<QDebug>

MathSolve::MathSolve()
{

}


//计算行列式
double MathSolve::getA(double arcs[3][3],int n)
{
    if(n==1)
    {
        return arcs[0][0];
    }
    double ans = 0;
    double temp[3][3]={0.0};
    int i,j,k;
    for(i=0;i<n;i++)
    {
        for(j=0;j<n-1;j++)
        {
            for(k=0;k<n-1;k++)
            {
                temp[j][k] = arcs[j+1][(k>=i)?k+1:k];

            }
        }
        double t = getA(temp,n-1);
        if(i%2==0)
        {
            ans += arcs[0][i]*t;
        }
        else
        {
            ans -=  arcs[0][i]*t;
        }
    }
    return ans;
}

// 计算伴随矩阵
void  MathSolve::getAStart(double arcs[3][3],int n,double ans[3][3])
{
    if(n==1)
    {
        ans[0][0] = 1;
        return;
    }
    int i,j,k,t;
    double temp[3][3];
    for(i=0;i<n;i++)
    {
        for(j=0;j<n;j++)
        {
            for(k=0;k<n-1;k++)
            {
                for(t=0;t<n-1;t++)
                {
                    temp[k][t] = arcs[k>=i?k+1:k][t>=j?t+1:t];
                }
            }


            ans[j][i]  =  getA(temp,n-1);  //此处顺便进行了转置
            if((i+j)%2 == 1)
            {
                ans[j][i] = - ans[j][i];
            }
        }
    }
}


bool MathSolve::GetMatrixInverse(double src[3][3],int n,double des[3][3])
{
    double flag=getA(src,n);
    double t[3][3];
    if(0==flag)
    {
        qDebug()<< "原矩阵行列式为0，无法求逆。请重新运行";
        return false;//如果算出矩阵的行列式为0，则不往下进行
    }
    else
    {
        getAStart(src,n,t);
        for(int i=0;i<n;i++)
        {
            for(int j=0;j<n;j++)
            {
                des[i][j]=t[i][j]/flag;
            }

        }
    }

    return true;
}

// 参数  起始点和目标点  要一一对应
double * MathSolve::solve(double begin_points[3][2], double end_points[3][2]){
    //重组矩阵
    double T[3][3];
    for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++){
            if(j != 2)  T[i][j] = begin_points[i][j];
            else        T[i][j] = 1;
        }

    // 计算矩阵T的逆
    double T_inv[3][3];
    bool flag = GetMatrixInverse(T, 3, T_inv);
    if(!flag)  return NULL;

    // 计算最终结果, 按照a11 a12 a21 a22 a31 a32的顺序
    double *res = (double *) malloc(sizeof(double) * 6);
    for(int i = 0; i < 3; i++){
        for(int j = 0; j < 2; j++){
            int index = 2 * i + j;
            res[index] = 0;
            for(int k = 0; k < 3; k++){
                res[index] += T_inv[i][k] * end_points[k][j];
            }
        }
    }
    return res;
}

 double MathSolve::cal2X2Det(double a[2][2])
 {
     double result;
     result = a[0][0]* a[1][1]- a[0][1]* a[1][0];
     qDebug()<<"Det"<<result;
     return result;

 }

  int  MathSolve::calSign(double b)
  {
      if(b<0) return -1;
       if(b==0) return 0;
       if(b>0) return 1;
  }



