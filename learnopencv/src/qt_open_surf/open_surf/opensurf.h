#ifndef OPENSURF_H
#define OPENSURF_H

#include <QDialog>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/calib3d/calib3d.hpp>
#include "ipoint.h"
#include "kmeans.h"

using namespace cv;

namespace Ui {
class OpenSurf;
}

class OpenSurf : public QDialog
{
    Q_OBJECT
    
public:
    explicit OpenSurf(QWidget *parent = 0);
    ~OpenSurf();
    
private slots:
    void on_openButton_clicked();

    void on_detectButton_clicked();

    void on_matchButton_clicked();

    void on_closeButton_clicked();

    void on_clusterButton_clicked();

private:
    Ui::OpenSurf *ui;
    IplImage *img1, *img2, *img_match1, *img_match2;
    IpVec ipts, ipts1, ipts2;
    IpPairVec matches;
    Kmeans km;
    int open_image_num;

};

#endif // OPENSURF_H
