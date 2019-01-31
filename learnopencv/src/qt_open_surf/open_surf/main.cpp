#include <QApplication>
#include "opensurf.h"

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    OpenSurf w;
    w.show();
    
    return a.exec();
}
