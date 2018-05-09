#ifndef SEGBYCUT_INNERPROJECTION_H
#define SEGBYCUT_INNERPROJECTION_H

#include "NaiveSCellFilter.h"

class InnerProjection
{

public:
    typedef DGtal::Z2i::SCell SCell;
    typedef DGtal::Z2i::Point Point;
    typedef DGtal::Z2i::Domain Domain;

    typedef DGtal::ImageContainerBySTLVector<Domain, unsigned char> Image2D;

public:

    InnerProjection(KSpace& KImage):KImage(KImage){}

    void operator()(Curve& projection, Curve& toProject, Curve& projectTo)
    {
        NaiveSCellFilter toProjectFilter(KImage,toProject);
        NaiveSCellFilter toProjectFilterCopy = toProjectFilter;

        NaiveSCellFilter projectToFilter(KImage,projectTo);

        toProjectFilterCopy- (toProjectFilter-projectToFilter);
        toProjectFilterCopy.computeBoundary(projection);

        displayProjection(projection,toProject,projectTo);
    }

private:

    void displayProjection(Curve& projection, Curve& toProject, Curve& projectTo)
    {

        Board2D board;
        board << DGtal::CustomStyle(toProject.begin()->className(),
                                    new DGtal::CustomColors( DGtal::Color::Gray, DGtal::Color::Gray ) );
        board <<  toProject;

        board << DGtal::CustomStyle(toProject.begin()->className(),
                                    new DGtal::CustomColors( DGtal::Color::Green, DGtal::Color::Green ) );
        board << projectTo;

        board << DGtal::CustomStyle(toProject.begin()->className(),
                                    new DGtal::CustomColors( DGtal::Color::Red, DGtal::Color::Red ) );
        board << projection;

        board.saveEPS("innerProjection.eps");
    }

private:
    KSpace& KImage;
};

#endif //SEGBYCUT_INNERPROJECTION_H
