#ifndef SEGBYCUT_NAIVECURVEFILTER_H
#define SEGBYCUT_NAIVECURVEFILTER_H

#include <DGtal/helpers/StdDefs.h>

class NaiveSCellFilter
{
public:
    typedef DGtal::Z2i::SCell SCell;
    typedef DGtal::Z2i::Point Point;
    typedef DGtal::Z2i::Domain Domain;
    typedef DGtal::Z2i::KSpace KSpace;
    typedef DGtal::Z2i::Curve Curve;

    typedef Curve::ConstIterator SCellIterator;

    typedef DGtal::ImageContainerBySTLVector<Domain, unsigned char> Image2D;

    typedef enum {NotMarked=0,Marked=1} MarkerStatus;

public:
    NaiveSCellFilter(KSpace& KImage, Curve& baseCurve):KImage(KImage)
    {
        width = KImage.upperBound()[0]+1;
        height = KImage.upperBound()[1]+1;
        size = width*height;

        mask.resize(size);

        for(int i=0;i<size;++i) mask[i] = NotMarked;

        for(auto it=baseCurve.begin();it!=baseCurve.end();++it)
        {
            SCell pixel = KImage.sDirectIncident( *it, KImage.sOrthDir(*it) );
            Point p = KImage.sCoords(pixel);
            int x = p[0];
            int y = p[1];

            int k = y*width + x;
            mask[k] = Marked;
        }

        SCell innerPixel;
        {
            int x, y;
            auto it = baseCurve.begin();
            do {
                SCell curveSCell = *it;
                Dimension dir = KImage.sOrthDir(curveSCell);
                SCell boundaryPixel = KImage.sDirectIncident(curveSCell, dir);

                int sign;
                if (dir == 0) sign = KImage.sSign(curveSCell) ? 1 : -1;
                else sign = KImage.sSign(curveSCell) ? -1 : 1;

                innerPixel = KImage.sGetAdd(boundaryPixel, dir, sign);

                x = KImage.sCoords(innerPixel)[0];
                y = KImage.sCoords(innerPixel)[1];

                ++it;
            } while (mask[y * width + x] == 1);
        }

        fillIn( KImage.sCoords( innerPixel) );

    }

    NaiveSCellFilter& operator=(const NaiveSCellFilter& other)
    {
        KImage = other.KImage;
        height = other.height;
        width = other.width;
        size = other.size;
        mask = other.mask;

        return *this;
    }

    void fillIn(Point innerPixel)
    {
        Point exploreFilter[4] = { Point(1,0),Point(0,1),Point(-1,0),Point(0,-1) };
        std::stack<Point> toExplore;

        Point myInnerPixel( innerPixel[0],innerPixel[1]);

        toExplore.push(myInnerPixel);
        while(!toExplore.empty())
        {
            Point next = toExplore.top();
            toExplore.pop();

            int x = next[0];
            int y = next[1];

            if(x<0 || y <0) continue;
            if(x>=width || y >=height) continue;

            int k = y*width+x;

            if(mask[k]==Marked) continue;
            mask[k] = Marked;

            for(int i=0;i<4;++i)
            {
                toExplore.push( next + exploreFilter[i] );
            }
        }

    }

    NaiveSCellFilter& operator+(const NaiveSCellFilter& other)
    {

        int i=0;
        for(auto it=other.begin();it!=other.end();++it)
        {
            mask[i] = *it;
            ++i;
        }

        return *this;
    }

    NaiveSCellFilter& operator-(const NaiveSCellFilter& other)
    {
        int i=0;
        for(auto it=other.begin();it!=other.end();++it)
        {

            if( *it == Marked )
            {
                mask[i] = NotMarked;
            }

            ++i;
        }

        return *this;
    }


    std::vector<char>::const_iterator begin() const
    {
        return mask.begin();
    }

    std::vector<char>::const_iterator end() const
    {
        return mask.end();
    }

    void computeBoundary(Curve& boundary)
    {
        Domain domain(Point(0,0),Point(width,height));
        Image2D image(domain);

        for(int i=0;i<size;++i)
        {
            if(mask[i]==1)
            {
                int y = i/width;
                int x = i%width;

                image.setValue( Point(x,y), 255 );
            }
        }

        ImageProc::computeBoundaryCurve(image,boundary,100);
    }

    void print()
    {
        for(int i=0;i<size;++i)
        {
            if(i%width==0) printf("\n");
            printf("%d",mask[i]);
        }
    }

private:
    KSpace& KImage;
    int height, width,size;

    std::vector<char> mask;
};

#endif //SEGBYCUT_NAIVECURVEFILTER_H
