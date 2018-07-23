#ifndef SEGBYCUT_PREPROCESSIMAGE_H
#define SEGBYCUT_PREPROCESSIMAGE_H

#include <string>
#include <DGtal/kernel/sets/DigitalSetInserter.h>
#include "../ImageProc/imageProc.h"
#include "../ImageProc/ImageProc.h"

class ImageData
{
public:
    typedef DGtal::Z2i::DigitalSet DigitalSet;
    typedef DGtal::Z2i::Domain Domain;
    typedef DGtal::Z2i::Point Point;

    typedef DGtal::ImageContainerBySTLVector<Domain, unsigned char> Image2D;
    typedef DGtal::DigitalSetInserter<DigitalSet> DigitalSetInserter;

private:
    struct ProcessedDigitalSets
    {
        ProcessedDigitalSets(Domain domain):originalDS(domain),
                                            segmentationDS(domain){}

        DigitalSet originalDS;
        DigitalSet segmentationDS;
    };

public:

    ImageData(std::string imageFilepath):originalImage(GenericReader<SegCut::Image2D>::import(imageFilepath)),
                                         preprocessedImage(originalImage.domain()),
                                         digitalSets(originalImage.domain())
    {
        preprocessingTwo(originalImage);
        preprocessingOne(preprocessedImage);
        //preprocessedImage = originalImage;
    }

private:
    void preprocessingOne(Image2D& image)
    {
        ImageProc::closing(image, image, 1);

        int borderWidth = 10;

        Domain newDomain(image.domain().lowerBound(),
                         image.domain().upperBound() + DGtal::Z2i::Point(2 * borderWidth, 2 * borderWidth)
        );

        SegCut::Image2D paddedImage(newDomain);
        ImageProc::createBorder(paddedImage, image, borderWidth);
        image = paddedImage;
    }

    void preprocessingTwo(Image2D& image)
    {
        ImageProc2::ImageAsDigitalSet(digitalSets.originalDS,image);

        DigitalSet dilated( image.domain() );
        ImageProc2::Dilate(dilated, digitalSets.originalDS,ImageProc2::Dilate::RECT,1);

        DigitalSet dilatedNoHoles( image.domain() );
        ImageProc2::NoHoles(digitalSets.segmentationDS,dilated);

        saveDigitalSetAsImage(preprocessedImage,digitalSets.segmentationDS);
    }

    void saveDigitalSetAsImage(Image2D& image, DigitalSet& ds)
    {
        for(auto it=image.domain().begin();it!=image.domain().end();++it)
        {
            image.setValue(*it,0);
        }

        for(auto it=ds.begin();it!=ds.end();++it)
        {
            image.setValue(*it,255);
        }
    }


public:
    Image2D originalImage;
    Image2D preprocessedImage;

    ProcessedDigitalSets digitalSets;
};

#endif //SEGBYCUT_PREPROCESSIMAGE_H
