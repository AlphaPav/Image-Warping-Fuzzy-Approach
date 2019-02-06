#ifndef MATHSOLVE_H
#define MATHSOLVE_H


class MathSolve
{
public:
    MathSolve();
    double getA(double arcs[3][3],int n);
    void  getAStart(double arcs[3][3],int n,double ans[3][3]);
    bool GetMatrixInverse(double src[3][3],int n,double des[3][3]);
    double * solve(double begin_points[3][2], double end_points[3][2]);
    double cal2X2Det(double a[2][2]);
    int calSign(double b);
};

#endif // MATHSOLVE_H
