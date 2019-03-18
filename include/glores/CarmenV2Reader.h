#pragma once


#include <iostream>
#include <string>
#include <fstream>

#include <glores/geometry.h>

namespace glores {

    struct RobotLaserMess {
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        VectorPoint2d points;
        Eigen::Vector3d wTb;
        Eigen::Vector3d wTl;
        int vertexId;
        double timestamp;
    };

    class CarmenV2Reader {
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        std::vector<RobotLaserMess> scans;

        /**
         * Creates an empty reader of CarMeN (Carnegie Mellon Navigation toolkit) 
         * modified log files.
         */
        CarmenV2Reader();

        /**
         * Reads scan from reader of CarMeN (Carnegie Mellon Navigation toolkit) 
         * modified log files.
         * @param filename the log file to be read.
         */
        CarmenV2Reader(const std::string& filename);

        /**
         * Reads scan from reader of CarMeN (Carnegie Mellon Navigation toolkit) 
         * modified log files.
         * @param filename the log file to be read.
         */
        void readCarmen(const std::string& filename);
    };

} // namespace 