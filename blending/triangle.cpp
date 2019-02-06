#include "triangle.h"
#include <math.h>
#include <iostream>
#include<QDebug>
# define MY_PI   3.14159265358979323846
Triangle::Triangle( QPointF _pos0,QPointF _pos1, QPointF _pos2)
{
    pos1=_pos1;
    pos0= _pos0;
    pos2= _pos2;
//    qDebug() <<_pos0 <<_pos1<<_pos2;

    e1=calLineLen(pos1,pos0);
    e2=calLineLen(pos1,pos2);
    e0= calLineLen(pos0,pos2);
//   qDebug() <<"e1"<<e1 <<e2<<e0;
    a1=(e1*e1+e2*e2-e0*e0)/(2*e1*e2);
//       qDebug() <<"a1"<<a1;
    a1= acos (a1)* 180/ MY_PI ;
//      qDebug() <<"a1"<<a1;
    a0=(e0*e0+e1*e1-e2*e2)/(2*e1*e0);
    a0= acos (a0)* 180/ MY_PI ;
    a2=(e0*e0+e2*e2-e1*e1)/(2*e2*e0);
    a2= acos (a2)* 180/ MY_PI;

}
float Triangle::calCosAngle(QPointF p0, QPointF p1, QPointF p2)
{
    float l1=calLineLen(p1,p0);
    float l2=calLineLen(p1,p2);
    float l0= calLineLen(p0,p2);
    float angle=(l1*l1+l2*l2-l0*l0)/(2*l1*l2);
    angle= acos (angle)* 180/ MY_PI ;
    return angle;
}

float Triangle::calLineLen(QPointF p1, QPointF p2)
{
    float result= (p1.x()-p2.x())*(p1.x()-p2.x())+
            (p1.y()-p2.y())*(p1.y()-p2.y());
    result=sqrt(result);
    return result;
}
QPointF Triangle::calMiddleVector()
{
    QPointF middlepoint = (pos0+pos2)/2;
    return middlepoint-pos1;

}

float Triangle::calArea()
{
    float angle = a1/180 * MY_PI;
    float A= 0.5 * sin(angle)* e1* e2 ;
//       qDebug() <<"calArea"<<A;
    return A;
}
