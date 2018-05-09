#ifndef SEGBYCUT_CURVEPROJECTOR_H
#define SEGBYCUT_CURVEPROJECTOR_H

#include <DGtal/helpers/StdDefs.h>

class CurveProjector
{
public:
    typedef DGtal::Z2i::KSpace KSpace;
    typedef DGtal::Z2i::Curve Curve;

    typedef KSpace::SCell SCell;
    typedef std::map<SCell,SCell> Projection;

    typedef Curve::ConstIterator SCellIterator;

    class UnsignedSCellComparison{
    public:
        bool operator()(const SCell& s1, const SCell& s2) const{
            return s1.preCell().coordinates < s2.preCell().coordinates;
        }
    };

public:
    CurveProjector(KSpace& KImage,
                   Curve& toProject,
                   Curve& projectOn):KImage(KImage),
                                     toProject(toProject),
                                     projectOn(projectOn),
                                     toProjectSCells(UnsignedSCellComparison()),
                                     projectOnSCells(UnsignedSCellComparison())
    {
        std::for_each(toProject.begin(),toProject.end(), [this](const SCell& s){this->toProjectSCells.insert(s);});
        std::for_each(projectOn.begin(),projectOn.end(), [this](const SCell& s){this->projectOnSCells.insert(s);});

        notMappedSCell = KImage.sCell( DGtal::PointVector<2,int>(1,1) );    //Notice it is a pixel
    }

    void operator()(Projection& p)
    {
        SCellIterator it = toProject.begin();
        do
        {
            if( projectOnSCells.find(*it)!=projectOnSCells.end() )
            {
                p[*it] = *it;
            }
            else
            {

                int dist = 1;
                bool chooseC1 = false;
                bool chooseC2 = false;
                SCell candidate1;
                SCell candidate2;
                while(dist<8 && !chooseC1 && !chooseC2)
                {
                    DGtal::Dimension  dir = KImage.sOrthDir(*it);

                    if(candidate1.preCell().coordinates[dir]-dist < KImage.lowerBound()[dir] ||
                       candidate1.preCell().coordinates[dir]+dist > KImage.upperBound()[dir] ) break;

                    if(candidate2.preCell().coordinates[dir]-dist < KImage.lowerBound()[dir] ||
                       candidate2.preCell().coordinates[dir]+dist > KImage.upperBound()[dir] ) break;


                    candidate1 = KImage.sGetAdd(*it, dir, dist);
                    candidate2 = KImage.sGetAdd(*it, dir, -dist);


                    if (projectOnSCells.find(candidate1) != projectOnSCells.end()) chooseC1 = true;
                    if (projectOnSCells.find(candidate2) != projectOnSCells.end()) chooseC2 = true;

                    ++dist;
                }

                if( chooseC1 && chooseC2)
                {
                    throw std::runtime_error("Two possible projections.");
                }
                else if (chooseC1)
                {
                    p[*it] = candidate1;
                }
                else if(chooseC2)
                {
                    p[*it] = candidate2;
                }
                else if( !chooseC1 && !chooseC2)
                {
                    p[*it] = notMappedSCell;
                }

            }

            ++it;
        }while(it!=toProject.end());
    }

    bool validProjection(SCell& s)
    {
        return s!=notMappedSCell;
    }

private:
    void displayCurves()
    {
        Board2D board;

        board << CustomStyle( toProject.begin()->className(), new CustomColors( DGtal::Color::Red, DGtal::Color::Black ) );
        board << toProject;
        board << CustomStyle( projectOn.begin()->className(), new CustomColors( DGtal::Color::Green, DGtal::Color::Black ) );
        board << projectOn;
        board << CustomStyle( projectOn.begin()->className(), new CustomColors( DGtal::Color::Purple, DGtal::Color::Black ) );
        board << KImage.sCell( DGtal::PointVector<2,int>(111,90),true);
        board.saveEPS("dilationCurves.eps");
    }

private:
    KSpace& KImage;
    Curve& toProject;
    Curve& projectOn;

    std::set<SCell,UnsignedSCellComparison> toProjectSCells;
    std::set<SCell,UnsignedSCellComparison> projectOnSCells;

    SCell notMappedSCell;
};

#endif //SEGBYCUT_CURVEPROJECTOR_H
