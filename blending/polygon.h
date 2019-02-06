#ifndef POLYGON_H
#define POLYGON_H

#include<QPainter>
# include "triangle.h"
class Polygon
{
public:
    Polygon(QVector<QPointF> _points, int _num);
    int mNum;//实际点的数量
    int n; //课件上-1,故n= mNum-1
    float mArea;

    QVector<QPointF> mPoints;
    QVector<QPointF> mLocalPos;

    QVector<float> mLineLens;
    QVector<Triangle*> mTris; //每个点 point 对应一个 triangle
    float calLineLen(QPointF p1, QPointF p2);
    void calArea();

};



#endif // POLYGON_H
