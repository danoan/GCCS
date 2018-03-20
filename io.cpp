#include "io.h"


namespace IO {
    void listFiles(VectorOfPath &files,
                   std::string directory) {

        path p(directory);
        if (boost::filesystem::exists(p)) {
            if (is_directory(p)) {
                for (auto di = directory_iterator(p); di != directory_iterator(); ++di) {
                    files.push_back(*di);
                }
            } else {
                throw filesystem_error("Input folder is not a directory.",
                                       boost::system::errc::make_error_code(boost::system::errc::invalid_argument));
            }
        }
    }

    std::string saveImage(Image2D& out,
                          std::string outputFolder,
                          std::string suffix)
    {

        std::string imageOutputPath = outputFolder + "/" + suffix + ".pgm";

        boost::filesystem::path p2(imageOutputPath.c_str());
        p2.remove_filename();
        boost::filesystem::create_directories(p2);

        GenericWriter<Image2D>::exportFile(imageOutputPath.c_str(),out);
        return imageOutputPath;
    }
}