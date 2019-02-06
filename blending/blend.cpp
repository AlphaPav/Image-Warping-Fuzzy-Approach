#include "blend.h"
#include <iostream>
#include <qDebug>
#include <algorithm>
#include <math.h>
# define MYPI 3.1415926
Blend::Blend(Polygon * _Pstart, Polygon * _Pend, int density)
{
     // density = 1;
    Pstart= _Pstart;
    Pend= _Pend;
    mDensity=density;

    calFuzzyMinGraph();
    calAffineKeyPoints();
    calAffineMatrix();
    calPolygonLocalPos();

    for(int i=0;i<= density;i++)
    {
        blendPolygon(i*1.0/density);
    }
}

Blend::~Blend()
{
    free(Pstart);
    free(Pend);


}

float Blend::SimT(Triangle *T0, Triangle* T1)
{
    float w1= 0.5;
    float w2= 1-w1;
    float sim=0;
    sim += w2* (1- (fabs(T0->a1 - T1->a1))/360);
    sim += w1* (1 - (fabs(T0->e1*T1->e2- T0->e2* T1->e1))/(T0->e1*T1->e2+T0->e2* T1->e1));
    return sim;
}

float Blend::SimP(Polygon* P0, Polygon* P1,QVector<int> MAP)
{
    int n= P0->n;
    float sum=0;
    int j;
    for(int i =0; i<=n;i++)
    {
        j= MAP.at(i);//对应目标多边形的顶点下标
        sum+=  SimT(P0->mTris.at(i),P1->mTris.at(j));//计算对应三角形相似程度
    }
    sum=sum/(n+1);
    return sum;
}
void Blend::calFuzzyMinGraph()
{
    float _INFINITY= 100000;
    int m= Pstart->n;
    int n= Pend->n;
    qDebug()<<m<<n;
    float** minGraph=(float**)malloc((sizeof (float *)) *(m+2) );
    float** Graph= (float**)malloc((sizeof (float *)) *(m+2) );
    Path** PathGraph= (Path**)malloc((sizeof (Path *)) *(m+2) );
    float** ValueGraph= (float**)malloc((sizeof (float *)) *(m+2) );
    QVector<int*> MapList;
    QVector<float> SumList;
    int i;
    for( i= 0; i< m+2;i++)
    {
        minGraph[i]= (float*) malloc((sizeof(float)*(2*n+3)));
        PathGraph[i]=(Path*) malloc((sizeof(Path)*(n + 2)));
        ValueGraph[i]=(float*) malloc((sizeof(float)*(n + 2)));
        Graph[i]= (float*) malloc((sizeof(float)*(n+2)));
    }
    int j;
    for( i= 0; i< m+2;i++)
    {
        for(j= 0; j<2*n+2;j++)
        {
            if(i== m+1)
            {
                minGraph[i][j]=minGraph[0][j];
            }else
            {   //余图
                minGraph[i][j]= 1 -SimT(Pstart->mTris.at(i), Pend->mTris.at(j%(n+1)));
            }
        }
        j = 2*n+2;
        if(i== m+1)
        {
            minGraph[i][j]=minGraph[0][j];
        }else{
            minGraph[i][j]= minGraph[i][0];
        }
    }
//    for( i= 0; i< m+2;i++)
//    {
//        for(j= 0; j<=2*n+2;j++)
//        {
//            qDebug("%f ",minGraph[i][j]);
//        }
//        qDebug("\n" );
//    }
//    qDebug("********minGraph Done*******" );

     int k;
     float min_sim=10000;
     int min_map_index;
     for(k=0;k<=n;k++) //起点
     {
        int* tempMap = (int*)malloc(sizeof(int)* (m+1));
        for (i=0;i<=m;i++)
        {
            for(j=0;j<=n;j++)
            {
                //初始化
                Graph[i][j]= minGraph[i][j+k];

                ValueGraph[i][j]= _INFINITY;
                PathGraph[i][j].lastCol=-1; PathGraph[i][j].lastRow=-1;
            }
        }

        for (i = 0; i <= m; i++)
        {
            int bound = i;
            if (bound > n) bound = n;
            j = 0;
            if (i == 0)
            {
                ValueGraph[i][ j] = Graph[i][ j];
                PathGraph[i][j].lastCol = -1; PathGraph[i][ j].lastRow = -1;
                continue;
            }
            else
            {   //j=0  只能上方传过来
                ValueGraph[i][j] = Graph[i][j] + ValueGraph[i - 1][ j];
                PathGraph[i][j].lastCol = j; PathGraph[i][j].lastRow = i - 1;
            }

            if (bound >= 1) //bound=0 只能从上方传过来，已经处理过了
            {
                for (j = 1; j <= bound; j++)
                {
                    if ((j == bound) &&(i<=n))  //只能接收左上方传来的
                    {
                        ValueGraph[i][ j] = ValueGraph[i - 1][ j - 1] + Graph[i][ j]; //左边斜上方过来的
                        PathGraph[i][ j].lastCol = j - 1; PathGraph[i][ j].lastRow = i - 1;
                    }
                    else
                    {
                       float temp1 = ValueGraph[i - 1][ j] + Graph[i][ j];   //上面的
                        float temp2 = ValueGraph[i - 1][ j - 1] + Graph[i][ j]; //左边斜上方过来的
                        if (temp1 < temp2)
                        {
                            ValueGraph[i][j] = temp1;
                            PathGraph[i][ j].lastCol = j; PathGraph[i][ j].lastRow = i - 1;
                        }
                        else
                        {
                            ValueGraph[i][j] = temp2;
                            PathGraph[i][ j].lastCol = j - 1; PathGraph[i][ j].lastRow = i - 1;
                        }
                    }
                }
            }
        }
        qDebug("tempPathSumValue: %f \n", ValueGraph[m][n]);
        if(ValueGraph[m][ n]< min_sim)
        {
            min_sim= ValueGraph[m][n];
            min_map_index= k;//源多边形的起点
        }
        tempMap[m] =(n+k)%(n+1);
        int lastCol = PathGraph[m][ n].lastCol;
        int lastRow = PathGraph[m][ n].lastRow;//应该是m-1
        for (i = m-1; i >= 0; i--)
        {   //lastRow 此时应该是i
            tempMap[i] = (lastCol+k)%(n+1);
            lastRow = PathGraph[i][lastCol].lastRow; //应该是i-1
            lastCol= PathGraph[i][ lastCol].lastCol; //更新位置
        }
        MapList.append(tempMap);
        SumList.append(ValueGraph[m][ n]);
    }
    qDebug()<<"minPathSumValue"<< SumList.at(min_map_index) << min_sim;
    qDebug()<<"map_index"<<min_map_index;

    Map= MapList.at(min_map_index);
    for(i=0;i<=m;i++  )
    {
        qDebug()<<Map[i];
    }
}


void Blend::BubbleSort(float  *p, int length, int * ind_diff)
{
    for (int m = 0; m < length; m++)
    {
        ind_diff[m] = m;
    }
    for (int i = 0; i < length; i++)
    {
        for (int j = 0; j < length- i - 1; j++)
        {
            if (p[j] < p[j + 1])
            {
                float temp = p[j];
                p[j] = p[j + 1];
                p[j + 1] = temp;

                int ind_temp = ind_diff[j];
                ind_diff[j] = ind_diff[j + 1];
                ind_diff[j + 1] = ind_temp;
            }
        }
    }
}


void Blend::calAffineKeyPoints()
{
    float  * SmoothArr= (float *) malloc(sizeof(float)*Pstart->mNum);
    int * SmoothIndex= (int *) malloc(sizeof(int)*Pstart->mNum);
    for(int i=0; i<=Pstart->n;i++)
    {
        SmoothIndex[i]= i;
        SmoothArr[i]= calAffineSmooth(i);
    }
    //冒泡排序，取前三个

    BubbleSort(SmoothArr,Pstart->mNum,SmoothIndex );
    for (int i=0; i<3;i++)
    {
        AffineIndex[i]= SmoothIndex[i];
        qDebug()<< SmoothIndex[i]<< SmoothArr[i];
    }

}
float Blend::calAffineSim(int i, int j)
{
    float S=0;
    float w1=0.5,w2=0.5;
    float temp;
    Triangle* Ti =Pstart->mTris.at(i);
    Triangle* Tj =Pend->mTris.at(j);
    temp=(abs(Ti->e0 - Tj->e0)+ abs(Ti->e1 - Tj->e1) +abs(Ti->e0 - Tj->e1));
    temp=temp/(abs(Ti->e0 + Tj->e0)+ abs(Ti->e1 + Tj->e1) +abs(Ti->e0 + Tj->e1));
    temp= 1- temp;
    S+= w1*temp;
    temp=(abs(Ti->a0 - Tj->a0)+ abs(Ti->a1 - Tj->a1) +abs(Ti->a0 - Tj->a1));
    temp=1-temp/180;
    S+= w2*temp;
    return S;
}


float Blend::calAffineR(int i, int j)
{
    Triangle* Ti =Pstart->mTris.at(i);
    Triangle* Tj =Pend->mTris.at(j);
    QPointF V1= Ti->calMiddleVector();
    QPointF V2= Tj->calMiddleVector();
    QPointF Vori; Vori.setX(0);Vori.setY(0);
    float R= Ti->calCosAngle(V1, Vori, V2);
    return R;
}


float Blend::calAffineA(int i, int j)
{
    float A=0;
    Triangle* Ti =Pstart->mTris.at(i);
    Triangle* Tj =Pend->mTris.at(j);
    A= Ti->calArea()+ Tj->calArea();
    A= A * 1.0/(Pstart->mArea + Pend->mArea);
    return A;
}

float Blend::calAffineSmooth(int i)
{
    float a=0.4, b=0.4, c=0.2, d=2;
    float smooth;
    int j= Map[i];
    smooth=calAffineSim(i,j)*a + b*(1-calAffineR(i,j)/180) +c*calAffineA(i,j);
    return smooth;

}

void  Blend::calAffineMatrix()
{
    double MatStart[3][2] ;
    double MatEnd[3][2] ;

     int i,j;
    for(int k =0; k<3;k++)
    {
        i= AffineIndex[k];
        qDebug()<<"AffineIndex"<<AffineIndex[k];
        j = Map[i];
        MatStart[k][0] = Pstart->mPoints.at(i).x();
        MatStart[k][1] = Pstart->mPoints.at(i).y();
        MatEnd[k][0] = Pend->mPoints.at(j).x();
        MatEnd[k][1] = Pend->mPoints.at(j).y();
    }
    MathSolve mMathhelper;
   double *res = mMathhelper.solve(MatStart, MatEnd);
    qDebug()<< "Aff result";
   for(int i = 0; i < 6; i++)
       qDebug()<< res[i];//计算最终结果, 按照a11 a12 a21 a22 a31 a32的顺序
    MatA[0][0]= res[0]; //a11
    MatA[0][1]= res[1]; //a12
    MatA[1][0]= res[2]; //a21
    MatA[1][1]= res[3]; //a22
    MatT[0] =res[4]; //a31
    MatT[1]= res[5]; //a32
    calMatrixBC();


}

void Blend::calMatrixBC()
{   MathSolve mMathhelper;
    double m[2][2];
    m[0][0]= MatA[1][1];
    m[0][1]= -MatA[1][0];
    m[1][0]=- MatA[0][1];
    m[1][1]=MatA[0][0];
    int w=  mMathhelper.calSign(mMathhelper.cal2X2Det(m));
    qDebug()<<"w" <<w;
    for(int i=0;i<2;i++)
    {
        for(int j= 0;j<2;j++)
        {
            MatB[i][j]= MatA[i][j]+ w * m[i][j];
        }
    }
    double mK_2= MatB[0][0]* MatB[0][0] + MatB[1][0]*  MatB[1][0];
     mK= sqrt(mK_2);

    if(mK<0)
    {
        mK= -mK;
    }
    for(int i=0;i<2;i++)
    {
        for(int j= 0;j<2;j++)
        {
            MatB[i][j] = MatB[i][j] * 1.0/ mK;
        }
    }
       qDebug()<< acos(MatB[0][0])<<asin(- MatB[0][1])<<asin(MatB[1][0])<<acos(MatB[1][1]);
    double mycos= MatB[0][0];
    double mysin= MatB[1][0];
    double myangle= abs( asin(MatB[1][0]));
    if(mycos > 0 && mysin >0)
    {
        myangle= myangle;
    }else if( mysin>0 && mycos <0)
    {
        myangle=MYPI- myangle;
    }else if(mysin<0 && mycos<0)
    {
        myangle= MYPI+myangle;
    }else if(mysin<0 && mycos>0)
    {
        myangle= 2*MYPI - myangle;
    }
    if(myangle > MYPI) myangle-=2*MYPI;

    mAngle= myangle;  // 仿射变换的初始旋转角度


    qDebug()<< cos(myangle) << -sin(myangle) <<  sin(myangle) << cos(myangle);

//    for(int i=0;i<2;i++)
//    {
//        for(int j= 0;j<2;j++)
//        {
//            MatB[i][j] = MatB[i][j] * 1.0  * mK;
//            qDebug()<< "MatB*k"<<MatB[i][j];
//        }
//    }

    double wb= 1* 1.0/ (MatB[0][0]* MatB[1][1]- MatB[0][1]* MatB[1][0]);
    double n[2][2];  //MatB的逆矩阵
    n[0][0]= MatB[1][1]* wb;
    n[0][1]= -MatB[0][1]* wb;
    n[1][0]= -MatB[1][0]* wb;
    n[1][1]=MatB[0][0]* wb;
    MatC[0][0]= n[0][0]* MatA[0][0]+  n[0][1]* MatA[1][0];
    MatC[0][1]= n[0][0]* MatA[0][1]+  n[0][1]* MatA[1][1];
    MatC[1][0]= n[1][0]* MatA[0][0]+  n[1][1]* MatA[1][0];
    MatC[1][1]= n[1][0]* MatA[0][1]+  n[1][1]* MatA[1][1];

//    qDebug() << MatB[0][0] *  MatC[0][0] + MatB[0][1] *  MatC[1][0] << MatA[0][0];
//     qDebug() << MatB[0][0] *  MatC[0][1] + MatB[0][1] *  MatC[1][1] << MatA[0][1];
//      qDebug() << MatB[1][0] *  MatC[0][0] + MatB[1][1] *  MatC[1][0] << MatA[1][0];
//       qDebug() << MatB[1][0] *  MatC[0][1] + MatB[1][1] *  MatC[1][1] << MatA[1][1];
}

void Blend::blendAffineMatix(double t)
{
     qDebug()<<  "time: " <<t;
    int k;
    double tMatA[2][2];
    float c11= MatC[0][0];
    float c12= MatC[0][1];
    float c21= MatC[1][0];
    float c22= MatC[1][1];
    float a31= MatT[0];
    float a32= MatT[1];


    for (k=0;k<3;k++)
    {
        int i= AffineIndex[k];
        double u= Pstart->mPoints.at(i).x();
        double v= Pstart->mPoints.at(i).y();
        qDebug()<<  "start  Pos" <<Pstart->mPoints.at(i);
        qDebug()<<  "End  Pos" << Pend->mPoints.at(Map[i]);
        double new_u=u;
        double new_v=v;
        tMatA[0][0]= 1-t+  cos(t*mAngle)*t*c11-  sin(t*mAngle)*t*c21;
        tMatA[0][1]=   cos(t*mAngle)*t*c12-   sin(t*mAngle)*t*c22;
        tMatA[1][0]=   sin(t*mAngle)*t*c11 +  cos(t*mAngle)*t*c21;
        tMatA[1][1]=   sin(t*mAngle)*t*c12 +  cos(t*mAngle)*t*c22+ 1-t;

        new_u= tMatA[0][0]*u+ tMatA[1][0]*v+ a31*t;
        new_v= tMatA[0][1]*u+ tMatA[1][1]*v +a32*t;

        QPointF _Point(new_u,new_v);
        _Point.setX(new_u);
        _Point.setY(new_v);

        TempAffinePos[k]=_Point;
    }
}

void Blend::calPolygonLocalPos()
{
    QPointF A1=Pstart->mPoints.at( AffineIndex[0]);
    QPointF B1=Pstart->mPoints.at( AffineIndex[1]);
    QPointF C1=Pstart->mPoints.at( AffineIndex[2]);
    QPointF A2=Pend->mPoints.at( Map[AffineIndex[0]]);
    QPointF B2=Pend->mPoints.at( Map[AffineIndex[1]]);
    QPointF C2=Pend->mPoints.at( Map[AffineIndex[2]]);
    int i;
    for (i=0; i<Pstart->mNum;i++)
    {
        QPointF localpoint= solveLocalPos(A1,B1,C1,Pstart->mPoints.at(i));
        Pstart->mLocalPos.append(localpoint);
    }
    for (i=0; i<Pend->mNum;i++)
    {
        QPointF localpoint= solveLocalPos(A2,B2,C2,Pend->mPoints.at(i));
        Pend->mLocalPos.append(localpoint);
    }

}


QPointF Blend::solveLocalPos(QPointF A1, QPointF B1, QPointF C1, QPointF X1)
{

    float u,v;
    float Ax= A1.x(), Ay= A1.y();
    float Bx= B1.x(), By= B1.y();
    float Cx= C1.x(), Cy= C1.y();
    float X=X1.x(), Y= X1.y();

    v= (Ax-Bx) * (Y-By) -(Ay- By)* (X-Bx);
    v/= (Cy-By)*(Ax-Bx)-(Cx- Bx)*(Ay-By);
     u= ((X-Bx)-v*(Cx-Bx))/ (Ax-Bx);
     QPointF point(u,v);
     point.setX(u);
     point.setY(v);
     return point;

}

QPointF  Blend::solveGlobalPos(QPointF A, QPointF B, QPointF C, float u, float v)
{

    QPointF X= B+ u*(A-B) + v*(C-B);
    X.setX(B.x()+u*(A.x()-B.x())+ v*(C.x()-B.x()));
     X.setY(B.y()+u*(A.y()-B.y())+ v*(C.y()-B.y()));
    return X;
}

void Blend::blendPolygon(double t)
{
      blendAffineMatix(t);//获得当前得仿射矩阵
      QPointF A= TempAffinePos[0];
      QPointF B= TempAffinePos[1];
      QPointF C= TempAffinePos[2];

      QVector<QPointF> tempblendPoints;

      float u1,v1,u2,v2,u,v;
      for(int i=0;i<Pstart->mNum;i++)
    {
          QPointF startPoint = Pstart->mLocalPos.at(i);//(u1,v1)
          QPointF endPoint = Pend->mLocalPos.at(Map[i]); //对应的(u2,v2)
          qDebug()<<"map"<<i<<Map[i] <<Pstart->mLocalPos.at(i) <<Pend->mLocalPos.at(Map[i]);
          u1= startPoint.x();
          v1=startPoint.y();
          u2=endPoint.x();
          v2= endPoint.y();
          u= (1-t)*u1+ t*u2;
          v= (1-t)*v1+ t*v2;
          QPoint localx(u,v);
          localx.setX(u); localx.setY(v);
          qDebug()<<"local"<<localx  << u<<v ;

          QPointF globalX= solveGlobalPos(A,B,C,u,v);
          qDebug()<<"start"<< Pstart->mPoints.at(i) <<"globalX"<< globalX<< "end"<< Pend->mPoints.at(Map[i]);

          tempblendPoints.append(globalX);

    }
      AllBlendPoints.append(tempblendPoints);


}




