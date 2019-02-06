#include "polygon.h"
#include <iostream>
#include<QDebug>
#include<math.h>
Polygon::Polygon(QVector<QPointF> _points, int _num)
{
    mNum= _num;
    n= _num-1;
    int i;
    for (i=0; i< _num;i++){
        mPoints.append(_points.at(i));
//        qDebug() <<_points.at(i).x()<< _points.at(i).y();
    }
    for(i=0 ;i<=n;i++)
    {
        if(i==n)
        {
            mLineLens.append(calLineLen(mPoints.at(n), mPoints.at(0)));
        }else
        {
            mLineLens.append(calLineLen(mPoints.at(i), mPoints.at(i+1)));
        }

        if(i==0)
        {
            Triangle * tri= new Triangle(mPoints.at(n),mPoints.at(0),mPoints.at(1));

            mTris.append(tri);
        }else if(i==n)
        {
             Triangle* tri= new Triangle(mPoints.at(n-1),mPoints.at(n),mPoints.at(0));
            mTris.append(tri);
        }else
        {
            Triangle* tri= new Triangle(mPoints.at(i-1),mPoints.at(i),mPoints.at(i+1));
            mTris.append(tri);
        }
    }
    calArea();
}
void Polygon::calArea()
{
    mArea= 0;
    if(mNum < 3) return;
    double s= 0;
    for (int i=0; i< mNum;i++)
    {
        s+= mPoints.at(i).x()* mPoints.at((i+1)%mNum).y()
                - mPoints.at(i).y()* mPoints.at((i+1)%mNum).x();
    }
    s= fabs(s/2.0);
    mArea= s;
}

float Polygon::calLineLen(QPointF p1, QPointF p2)
{
    float result= (p1.x()-p2.x())*(p1.x()-p2.x())+
            (p1.y()-p2.y())*(p1.y()-p2.y());
    result=sqrt(result);
    return result;
}
