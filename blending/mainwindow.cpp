#include "mainwindow.h"
#include "ui_mainwindow.h"
#include<QDebug>
#include <QFile>
MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    flag=1;
    this->setFixedSize(size().width(),size().height());
    tempPointNum=0;
    mDrawTimes=0;
    canPlay=false;
    isGen=false;

    isPlay=false;
    mBlend=NULL;
    tick=0;
    timer = new QTimer();
    connect(timer,SIGNAL(timeout()),this,SLOT(timerBlend()));
    this->setWindowTitle("Fuzzy 2D Blending");
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::paintEvent(QPaintEvent *ev)
{


    QFile qssfile(":/base.qss");//为UI控件设置样式
    qssfile.open(QFile::ReadOnly);
    QString qss;
    qss = qssfile.readAll();
    this->setStyleSheet(qss);
    QPainter p(this);
    //在界面的左边部分设置背景为白色，表示可接受鼠标点击的绘画区域
    p.setClipRect(5,0,size().width()-140,size().height());
    p.fillRect(QRectF(0, 0, size().width(), size().height()), Qt::white);

    p.setPen(QPen(QColor(200, 0 ,0), 2));
    p.setRenderHint(QPainter::Antialiasing, true);// 反走样

    for (int i=0 ;i<mPoints.size();i++)
    {
        p.drawEllipse(mPoints.at(i),5,5);
    }
    p.setPen(QPen(QColor(100,100,100), 2));



   if(mDrawTimes!=0)
    {
           //qDebug()<<"mmDrawTimes!=0"<<tick;
       for(int i=0;i<mStartP.size()-1;i++)
       {
            p.drawLine(mStartP.at(i), mStartP.at(i+1));//绘制初始图形关键点之间的直线
       }
        if(mStartP.size()>2)  p.drawLine(mStartP.at(mStartP.size()-1), mStartP.at(0));
    }
    if(mDrawTimes==2)
    {

        for(int i=0;i<mEndP.size()-1;i++)
        {
             p.drawLine(mEndP.at(i), mEndP.at(i+1));//绘制终止图形关键点之间的直线
        }
        if(mEndP.size()>2) p.drawLine(mEndP.at(mEndP.size()-1), mEndP.at(0));
    }

    if(isPlay)
    {
         //qDebug()<<"paintEventisPlay"<<tick;
        p.setPen(QPen(QColor(	255,140,0), 2));
        for (int i=0 ;i<mBlend->Pstart->mNum;i++)
        {
          //  qDebug()<<i<< mBlend->blendPoints.at(i);
            p.drawEllipse(mBlend->blendPoints.at(i),5,5);
        }
        p.setPen(QPen(QColor(100,100,100), 2));
        for (int i=0 ;i<mBlend->Pstart->mNum-1;i++)
        {

           // qDebug()<<i<< mBlend->blendPoints.at(i) << mBlend->blendPoints.at(i+1);
            p.drawLine(mBlend->blendPoints.at(i), mBlend->blendPoints.at(i+1));//绘制插值点之间的直线
        }
        if(mBlend->Pstart->mNum>2)
        {
           //  qDebug()<<mBlend->Pstart->mNum-1<< mBlend->blendPoints.at(mBlend->Pstart->mNum-1) << mBlend->blendPoints.at(0);
            p.drawLine(mBlend->blendPoints.at(mBlend->Pstart->mNum-1),mBlend->blendPoints.at(0));

        }
    }

}

void MainWindow::mousePressEvent(QMouseEvent *ev)
{
    if(flag==1)//可接收点击画点的模式
    {   //如果鼠标的点击位置是绘画区域
        if(mDrawTimes>2)
        {
            return;
        }
        if((ev->pos().x() <size().width()-140)&& (ev->pos().x() >0 )
                && (ev->pos().y() < (size().height())) &&  (ev->pos().y() >0))
        {
            mPoints.append(ev->pos());
            tempPointNum++;
            if(mDrawTimes==0)
            {
                mStartP.append(ev->pos());

            }else if(mDrawTimes==1)
            {
                mEndP.append(ev->pos());

            }
            update();
        }
    }
}
void MainWindow::timerBlend()
{
    isPlay=false;
    tick+=1;
    if(tick > mBlend->mDensity )
    {
        timer->stop();
        isPlay=false;
        return;
    }else if(tick>=0)
    {
        isPlay=true;
        mBlend->blendPoints= mBlend->AllBlendPoints.at(tick);
        update();
    }
}

void MainWindow::on_End_Draw_clicked()
{
    if(tempPointNum<=0)
    {
        return;
    }
    canPlay=false;
    mPointNum.append(tempPointNum);//此次绘制的点的个数
    mDrawTimes++;//绘制次数++
    update();
    tempPointNum=0;
    if(mDrawTimes==2)
    {
        qDebug("%d",mDrawTimes);
        qDebug("startpoint num:%d",mPointNum.at(0) );
        qDebug("endpoint num:%d",mPointNum.at(1) );
        if(mPointNum.at(0)<mPointNum.at(1) )
        {
            qDebug("You can't play. start points are less then end points");
            return;
        }
        canPlay=true;

    }else{
       qDebug("You doesn't finish drawing. You can't play.");
    }


}

void MainWindow::on_Play_clicked()
{
    if(canPlay)
    {
        if(isPlay) return;
        flag= 0; //停止绘制，不能接受屏幕上鼠标点击
        if(isGen==false)
        {
            if(mBlend!= NULL) delete mBlend;
            Polygon * P0= new Polygon(mStartP, mStartP.length());
            Polygon * P1= new Polygon(mEndP, mEndP.length());
            mBlend= new Blend(P0,P1,ui->Density->value());
            ui->Density->setReadOnly(true);
            isGen=true;

        }
        tick=-1;
        timer->stop();
        timer->start();
        qDebug()<<"Interval"<<ui->Interval->value();
        timer->setInterval(ui->Interval->value());
    }

}

void MainWindow::on_Clear_clicked()
{
    if(isPlay) return;
    clearAll();
    update();

}

void MainWindow::clearAll(){
    flag=1;
    mDrawTimes=0;
    tempPointNum=0;
    canPlay=false;
    isGen=false;
    isPlay=false;

    tick=0;
    QVector<int> pNullIntVector;
    mPointNum.swap(pNullIntVector);
    QVector<QPointF> pNullPointVector;
    mPoints.swap(pNullPointVector);
    QVector<QPointF> pNullPointVector1;
    mStartP.swap(pNullPointVector1);
    QVector<QPointF> pNullPointVector2;
    mEndP.swap(pNullPointVector2);

    if(mBlend!=NULL)
    {
        delete mBlend;
        mBlend=NULL;

    }

    ui->Density->setReadOnly(false);

}


