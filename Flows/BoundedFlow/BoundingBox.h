#ifndef SEGBYCUT_BOUNDINGBOX_H
#define SEGBYCUT_BOUNDINGBOX_H

#include <DGtal/helpers/StdDefs.h>
#include <DGtal/io/boards/Board2D.h>

class BoundingBox
{
public:
    typedef DGtal::Z2i::SCell SCell;
    typedef DGtal::Z2i::Point Point;
    typedef DGtal::Z2i::Domain Domain;
    typedef DGtal::Z2i::KSpace KSpace;
    typedef DGtal::Z2i::Curve Curve;

public:
    BoundingBox(KSpace& KImage):KImage(KImage){}

    void operator()(Point& lowerLeft, Point& upperRight, Curve& boundingBox,Curve& curveIn)
    {
        upperRight = Point(0,0);
        lowerLeft = Point(100,100);
        for(auto it=curveIn.begin();it!=curveIn.end();++it)
        {
            SCell pixel = KImage.sDirectIncident(*it,KImage.sOrthDir(*it));
            Point tryCoords = KImage.sKCoords(pixel);

            if(tryCoords[0] > upperRight[0]) upperRight[0] = tryCoords[0];
            if(tryCoords[1] > upperRight[1]) upperRight[1] = tryCoords[1];
            if(tryCoords[0] < lowerLeft[0]) lowerLeft[0] = tryCoords[0];
            if(tryCoords[1] < lowerLeft[1]) lowerLeft[1] = tryCoords[1];
        }

        typedef std::pair< SCell,bool > BoundaryPixel; //Pixel itself, it is on a corner

        std::vector<BoundaryPixel> boundaryPixels;
        SCell currentPixel = KImage.sCell(upperRight,true);

        int limitsCheck[4][3] = { {1,1,-1}, {1,0,-1}, {0,1,1}, {0,0,1} }; //Compare with LowerLeft, Use x, increment value

        bool corner = true;
        for(int i=0;i<4;++i)
        {
            Point comparePoint = limitsCheck[i][0]==1?lowerLeft:upperRight;
            int dim = limitsCheck[i][1]==1?0:1; //To use x or y
            int inc = limitsCheck[i][2];


            Point currentPoint = KImage.sKCoords(currentPixel);

            int limValue = comparePoint[dim];
            int coordValue = currentPoint[dim];


            do
            {
                boundaryPixels.push_back( BoundaryPixel(currentPixel,corner) );
                currentPixel = KImage.sGetAdd(currentPixel,dim,inc);
                coordValue+=inc*2;

                corner = false;
            }while(coordValue!=limValue);
            corner = true;

        }

        typedef std::pair<Point,bool> ModelPoint;
        std::vector<SCell> curveSCells;
        ModelPoint models[4] = { ModelPoint( Point(0,1),true ), //Upper Horizontal
                                 ModelPoint( Point(-1,0),true ), //Left Vertical
                                 ModelPoint( Point(0,-1),false ), //Bottom Horizontal
                                 ModelPoint( Point(1,0),false) }; //Right Vertical


        int k=0;
        int previousDim=0; //x

        boundaryPixels.push_back(boundaryPixels[0]);
        for(auto it=boundaryPixels.begin()+1;it!=boundaryPixels.end();++it)
        {
            SCell boundaryPixel = it->first;
            bool corner = it->second;

            if(corner)
            {
                curveSCells.push_back( KImage.sCell( KImage.sKCoords(boundaryPixel) + models[k].first,models[k].second ) );
                k+=1;
                k=k%4;
            }

            curveSCells.push_back( KImage.sCell( KImage.sKCoords(boundaryPixel) + models[k].first,models[k].second) );


        }

        boundingBox.initFromSCellsVector(curveSCells);
        displayBoundingBox(boundingBox,curveIn);
    }

private:
    void displayBoundingBox(Curve& boundingBox, Curve& originalCurve)
    {
        DGtal::Board2D board;
        board << DGtal::CustomStyle( boundingBox.begin()->className(),
                                     new DGtal::CustomColors(DGtal::Color::Green,DGtal::Color::Black) );
        board << boundingBox;
        board << DGtal::CustomStyle( originalCurve.begin()->className(),
                                     new DGtal::CustomColors(DGtal::Color::Red,DGtal::Color::Black) );
        board << originalCurve;
        board.saveEPS("boundingBox.eps");
    }

private:
    KSpace& KImage;
};

#endif //SEGBYCUT_BOUNDINGBOX_H
