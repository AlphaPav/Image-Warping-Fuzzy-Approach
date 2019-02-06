#ifndef TRIANGLE_H
#define TRIANGLE_H
#include<QPainter>

class Triangle
{
public:
    Triangle(QPointF _pos1, QPointF _pos0, QPointF _pos2);
    float a1; //顶角 0-360之间
    float a0;
    float a2;
    float e1; //相邻边
    float e2;
    float e0;
    QPointF pos1;
    QPointF pos0;
    QPointF pos2;
    float calLineLen(QPointF p1, QPointF p2);
    QPointF calMiddleVector();
    float calCosAngle(QPointF p0, QPointF p1, QPointF p2);
    float calArea();
};

#endif // TRIANGLE_H
