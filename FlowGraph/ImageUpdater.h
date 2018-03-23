#ifndef SEGBYCUT_IMAGEUPDATER_H
#define SEGBYCUT_IMAGEUPDATER_H

#include "FlowGraph.h"
#include "FlowGraphQuery.h"

class ImageUpdater
{
public:
    typedef DGtal::Z2i::SCell SCell;
    typedef DGtal::Z2i::Point Point;
    typedef DGtal::Z2i::Domain Domain;

    typedef DGtal::ImageContainerBySTLVector<Domain, unsigned char> Image2D;

public:
    ImageUpdater(KSpace& KImage):KImage(KImage){}

    void operator()(Image2D& imageSolution, FlowGraph& fg, Image2D& baseImage)
    {

        ListDigraph::NodeMap<bool> sourceNodesFilter(fg.graph(),false);
        FlowGraphQuery::sourceComponentNodes(fg,
                                             sourceNodesFilter);


        std::vector<SCell> pixelsInTheGraph;

        FlowGraphQuery::fetchFromNodes(fg,
                                      [ &fg ](ListDigraph::Node& n)->SCell{ return fg.pixel(n);},
                                      sourceNodesFilter,
                                      std::back_inserter(pixelsInTheGraph));


        imageSolution = baseImage;
        std::for_each(pixelsInTheGraph.begin(),
                      pixelsInTheGraph.end(),
                      [this,&imageSolution](SCell& s){ imageSolution.setValue( this->KImage.sCoords(s) ,255); });

        fillHoles(imageSolution,
                  KImage,
                  fg);
    }

private:
    void fillHoles(Image2D& out,
                   KSpace& KImage,
                   FlowGraph& fg)
    {
        ListDigraph::NodeMap<bool> sourceNodesFilter(fg.graph(),true);
        sourceNodesFilter[fg.source()] = false;
        sourceNodesFilter[fg.target()] = false;

        std::vector<SCell> pixelsInTheGraph;
        FlowGraphQuery::fetchFromNodes(fg,
                                       [&fg](ListDigraph::Node& n)->SCell{ return fg.pixel(n);},
                                       sourceNodesFilter,
                                       std::back_inserter(pixelsInTheGraph));

        for(auto it=pixelsInTheGraph.begin();it!=pixelsInTheGraph.end();++it)
        {
            Z2i::SCell pixel = *it;
            Z2i::Point p = KImage.sCoords(pixel);
            Z2i::SCells N = KImage.sNeighborhood(pixel);

            unsigned char v = out(p);
            unsigned char nv;
            bool hole=true;
            for(int j=1;j<5;++j){
                nv = out(KImage.sCoords(N[j]));
                if(v==nv){
                    hole=false;
                    break;
                }
            }
            if(hole){
                out.setValue(p,nv);
            }
        }
    }
private:
    KSpace& KImage;
};

#endif //SEGBYCUT_IMAGEUPDATER_H
