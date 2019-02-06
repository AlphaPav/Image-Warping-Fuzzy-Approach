#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QPixmap>
#include<QPixmapCache>
#include <QPropertyAnimation>
#include<QPainter>
#include<QMouseEvent>
#include<QTime>
#include<QTimer>
#include <QMainWindow>
# include "blend.h"
namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
     bool canPlay;
     bool isGen;

     int flag;
     bool isPlay;
     int tick;
     QVector<QPointF> mPoints;//鼠标输入的点
     QVector<QPointF> mStartP;//起始图形的控制点
     QVector<QPointF> mEndP;//终止图形的控制点
     QVector<int> mPointNum;
     int tempPointNum;
     int mDrawTimes;
     Blend *mBlend;

     void paintEvent(QPaintEvent *ev);
     void mousePressEvent(QMouseEvent *ev);
     void clearAll();

    ~MainWindow();

private slots:
     void timerBlend();

     void on_End_Draw_clicked();

     void on_Play_clicked();

     void on_Clear_clicked();

private:
      QTimer * timer;
    Ui::MainWindow *ui;
};

#endif // MAINWINDOW_H
