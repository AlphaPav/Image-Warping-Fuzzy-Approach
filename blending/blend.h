#ifndef BLEND_H
#define BLEND_H
#include<QPainter>
#include "polygon.h"
#include "triangle.h"
#include "mathsolve.h"
struct MatrixPath
{
    int lastRow;
    int lastCol;
};
typedef MatrixPath Path;

struct SmoothNode
{
  int index;
  float value;
};
typedef struct SmoothNode smoothNode;


class Blend
{
public:
    Blend(Polygon * _Pstart, Polygon * _Pend,int density);
    Polygon * Pstart; //PO
    Polygon * Pend;  //P1
    int *  Map; //一种map关系： j in P1 =map(i in P0)
    QVector<QPointF> blendPoints;
    QVector<QVector<QPointF>> AllBlendPoints;
    int mDensity;
    int AffineIndex[3];//仿射变换的3个点
    float   MatA[2][2];
    float MatB[2][2];
    float  MatC[2][2];
    float MatT[2];
    QPointF TempAffinePos[3];

    float mAngle; // Affine的初始旋转角度
    float mK; //算Matrix B 纯粹旋转角度的时候 除的系数
    int M; //P0的点数，M>=N
    int N; //P1的点数
    float SimT(Triangle *T0, Triangle* T1);
    float SimP(Polygon* P0, Polygon* P1, QVector<int> MAP);
    void calFuzzyMinGraph();
    void calAffineKeyPoints();
     void calAffineMatrix();
    float calAffineSim(int i, int j);
    float calAffineR(int i, int j);
    float calAffineA(int i, int j);
    float calAffineSmooth(int i);
    void BubbleSort(float  *p, int length, int * ind_diff);
    void calMatrixBC();
    void blendAffineMatix(double t);
    void blendPolygon(double t);
    void calPolygonLocalPos();
    QPointF solveLocalPos(QPointF A1, QPointF B1, QPointF C1, QPointF X1);
    QPointF solveGlobalPos(QPointF A, QPointF B, QPointF C, float u, float v);
    virtual ~Blend();
private:



};

#endif // BLEND_H
