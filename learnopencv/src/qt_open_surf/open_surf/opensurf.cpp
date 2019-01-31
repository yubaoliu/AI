#include "opensurf.h"
#include "ui_opensurf.h"
#include <QtGui>
#include <QtCore>
#include "surflib.h"

using namespace std;

OpenSurf::OpenSurf(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::OpenSurf)
{
    open_image_num = 0;
    ui->setupUi(this);
}

OpenSurf::~OpenSurf()
{
    delete ui;
}

void OpenSurf::on_openButton_clicked()
{
    QString img_name = QFileDialog::getOpenFileName(this, "Open Image", "../open_surf",
                                                           tr("Image Files(*.png *.jpeg *.jpg *.bmp)"));
       if(0 == open_image_num)
           ui->textBrowser->clear();
       open_image_num ++;
       if( 1 == open_image_num )
           {
               img1 = cvLoadImage(img_name.toAscii().data());
               img_match1 = cvLoadImage(img_name.toAscii().data());
               cvSaveImage("../open_surf/load_img1.jpg", img1);
               ui->textBrowser->setFixedSize(img1->width, img1->height);
               ui->textBrowser->insertHtml("<img src=../open_surf/load_img1.jpg>");
           }
       else if(2 == open_image_num)
           {
               img2 = cvLoadImage(img_name.toAscii().data());
               img_match2 = cvLoadImage(img_name.toAscii().data());
               cvSaveImage("../open_surf/load_img2.jpg", img2);
               ui->textBrowser->setFixedSize(img1->width+img2->width, std::max(img1->height, img2->height));
               //取消自动换行模式，让2幅图片水平显示
               ui->textBrowser->setWordWrapMode (QTextOption::NoWrap);
               ui->textBrowser->insertHtml("<img src=../open_surf/load_img2.jpg>");
           }
       else if(3 == open_image_num)
           {
               open_image_num = 0;
               ui->textBrowser->clear();
           }
}

void OpenSurf::on_detectButton_clicked()
{
    if( 1 == open_image_num )
        {
            //用surf对特征点进行检测
            surfDetDes(img1, ipts, false, 5, 4, 2, 0.0004f);
            //在图像中将特征点画出来
            drawIpoints(img1, ipts);
            cvSaveImage("../open_surf/detect_img1.jpg", img1);
            ui->textBrowser->clear();
            ui->textBrowser->setFixedSize(img1->width, img1->height);
            ui->textBrowser->insertHtml("<img src=../open_surf/detect_img1.jpg>");
        }
    else if (2 == open_image_num)
        {
            //用surf对特征点进行检测
            surfDetDes(img1, ipts1, false, 5, 4, 2, 0.0004f);
            //在图像中将特征点画出来
            drawIpoints(img1, ipts1);
            cvSaveImage("../open_surf/detect_img1.jpg", img1);
            //用surf对特征点进行检测
            surfDetDes(img2, ipts2, false, 5, 4, 2, 0.0004f);
            //在图像中将特征点画出来
            drawIpoints(img2, ipts2);
            cvSaveImage("../open_surf/detect_img2.jpg", img2);
            ui->textBrowser->clear();
            ui->textBrowser->insertHtml("<img src=../open_surf/detect_img1.jpg>");

            ui->textBrowser->setFixedSize(img1->width+img2->width, std::max(img1->height, img2->height));
            //取消自动换行模式，让2幅图片水平显示
            ui->textBrowser->setWordWrapMode (QTextOption::NoWrap);
            ui->textBrowser->insertHtml("<img src=../open_surf/detect_img2.jpg>");
        }
}

void OpenSurf::on_matchButton_clicked()
{
    if(2 == open_image_num)
        {
            getMatches(ipts1,ipts2,matches);
            for (unsigned int i = 0; i < matches.size(); ++i)
              {
                drawPoint(img_match1,matches[i].first);
                drawPoint(img_match2,matches[i].second);

                const int & w = img1->width;
                const int & h1 = img1->height;
                const int & h2 = img2->height;
                //这里因为我事先已经知道了图片的相对打开后显示的位置，所以在画匹配的直线时加了点常识
                //因此该方法不通用，只是适合本例中给的图片，最好的方法就像Rob Hess的sift算法那样
                //把2张图片合成一张，然后在一张图片上画匹配直线
                cvLine(img_match1,cvPoint(matches[i].first.x,matches[i].first.y),
                       cvPoint(matches[i].second.x+w,matches[i].second.y+std::abs(h1-h2)),
                       cvScalar(255,255,255),1);
                cvLine(img_match2,cvPoint(matches[i].first.x-w,matches[i].first.y-std::abs(h1-h2)),
                       cvPoint(matches[i].second.x,matches[i].second.y),
                       cvScalar(255,255,255),1);
              }
            cvSaveImage("../open_surf/match_img1.jpg", img_match1);
            cvSaveImage("../open_surf/match_img2.jpg", img_match2);
            ui->textBrowser->clear();
            ui->textBrowser->insertHtml("<img src=../open_surf/match_img1.jpg>");

            ui->textBrowser->setFixedSize(img1->width+img2->width, std::max(img1->height, img2->height));
            //取消自动换行模式，让2幅图片水平显示
            ui->textBrowser->setWordWrapMode (QTextOption::NoWrap);
            ui->textBrowser->insertHtml("<img src=../open_surf/match_img2.jpg>");
        }
}


void OpenSurf::on_clusterButton_clicked()
{
    for (int repeat = 0; repeat < 10; ++repeat)
      {

        km.Run(&ipts, 5, true);
        drawPoints(img1, km.clusters);

        for (unsigned int i = 0; i < ipts.size(); ++i)
        {
          cvLine(img1, cvPoint(ipts[i].x,ipts[i].y), cvPoint(km.clusters[ipts[i].clusterIndex].x ,km.clusters[ipts[i].clusterIndex].y),cvScalar(255,255,255));
        }
        cvSaveImage("../open_surf/kmeans_img1.jpg", img1);
        ui->textBrowser->clear();
        ui->textBrowser->setFixedSize(img1->width, img1->height);
        ui->textBrowser->insertHtml("<img src=../open_surf/kmeans_img1.jpg>");
      }
}


void OpenSurf::on_closeButton_clicked()
{
    close();
}


