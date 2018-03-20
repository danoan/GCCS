#ifndef SEGBYCUT_IO_H
#define SEGBYCUT_IO_H

#include "boost/filesystem.hpp"

#include <DGtal/images/ImageContainerBySTLVector.h>
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/writers/GenericWriter.h"

#include <string>


namespace IO
{

    using namespace boost::filesystem;
    using namespace DGtal;

    typedef boost::filesystem::path Path;
    typedef std::vector<Path> VectorOfPath;
    typedef DGtal::ImageContainerBySTLVector<DGtal::Z2i::Domain, unsigned char> Image2D;


    void listFiles(VectorOfPath& files,
                   std::string directory);

    std::string saveImage(Image2D& out,
                          std::string outputFolder,
                          std::string suffix);
}



#endif //SEGBYCUT_IO_H
