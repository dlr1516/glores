/* 
 * File:   testBBOptimizerNd.cpp
 * Author: dario
 *
 * Created on 12 October 2018, 12:00
 */

#include <iostream>
#include <Eigen/Dense>
#include <random>
#include <glores/PointSetGlobalRegistration.h>
#include <glores/geometry.h>
#include <glores/ParamMap.h>
#include <glores/fileutils.h>
#include <glores/thirdparty/gnuplot-iostream.h>
#include <glores/Profiler.h>

using namespace std;

void saturateRanges(glores::VectorPoint2d& points, double distMax);

void samplePoints(glores::VectorPoint2d& points, double sampleRate);

void printPointSets(std::ostream& out, const Eigen::VectorXd & pose, glores::VectorPoint2d pointSrc, glores::VectorPoint2d pointDst);

int main(int argc, char** argv) {
    std::string filenameCfg;
    std::string filenameSrc;
    std::string filenameDst;
    double scanDistMax, scanSampleRate;
    glores::VectorPoint2d pointsSrc;
    glores::VectorPoint2d pointsDst;
    glores::Point2d psrc, pdst;
    int pointNum, trialNum;
    Eigen::VectorXd poseTrue(3), poseLow(3), poseUpp(3), poseOpt(3), poseCoder(3);
    double funcLow, funcUpp;
    double optXTol, optYTol, optInlierRatio;
    double optRelaxX, optRelaxY, optRelaxT;
    int optDepthMax, optIterationMax, optVerboseLevel;
    bool optRelaxationBoundOn, optCheaBoundQueueOn;
    glores::Transformation2d transfTrue;
    glores::PointSetGlobalRegistration optimizer;
    //glores::FastPointSetGlobalRegistration optimizer;
    glores::ParamMap params;
    bool fileInputOn;

    // Random noise
    //std::random_device noiseRd;
    std::default_random_engine noiseRd{0};
    std::mt19937 noiseGen(0);
    double noiseSigma;

    params.read(argc, argv);
    params.getParam<std::string>("cfg", filenameCfg, std::string(""));
    if (!params.read(filenameCfg)) {
        std::cout << "Cannot open configuration file \"" << filenameCfg << "\": using default values" << std::endl;
    }

    params.read(argc, argv);
    params.getParam<double>("tx", poseTrue(0), double(2.1));
    params.getParam<double>("ty", poseTrue(1), double(-3.0));
    params.getParam<double>("rot", poseTrue(2), double(42.0));
    poseTrue(2) *= M_PI / 180.0;
    params.getParam<double>("txL", poseLow(0), double(1.8));
    params.getParam<double>("txU", poseUpp(0), double(2.4));
    params.getParam<double>("tyL", poseLow(1), double(-3.3));
    params.getParam<double>("tyU", poseUpp(1), double(-2.7));
    params.getParam<double>("rotL", poseLow(2), double(37.0));
    params.getParam<double>("rotU", poseUpp(2), double(56.0)); 
    poseLow(2) *= M_PI / 180.0;
    poseUpp(2) *= M_PI / 180.0;
    params.getParam<double>("optXTol", optXTol, double(0.2));
    params.getParam<double>("optYTol", optYTol, double(0.01));
    params.getParam<double>("optInlierRatio", optInlierRatio, double(1.0));
    params.getParam<double>("optRelaxX", optRelaxX, double(0.05));
    params.getParam<double>("optRelaxY", optRelaxY, double(0.05));
    params.getParam<double>("optRelaxT", optRelaxT, double(1.0));
    optRelaxT *= (M_PI / 180.0);
    params.getParam<int>("optDepthMax", optDepthMax, int(50));
    params.getParam<int>("optIterationMax", optIterationMax, int(-1));
    params.getParam<int>("optVerboseLevel", optVerboseLevel, int(0));
    params.getParam<bool>("optRelaxationBoundOn", optRelaxationBoundOn, bool(true));
    params.getParam<bool>("optCheaBoundQueueOn", optCheaBoundQueueOn, bool(true));
    params.getParam<double>("noiseSigma", noiseSigma, double(0.01));
    params.getParam<bool>("fileIn", fileInputOn, bool(false));
    params.getParam<std::string>("src", filenameSrc, std::string(""));
    params.getParam<std::string>("dst", filenameDst, std::string(""));
    params.getParam<double>("scanDistMax", scanDistMax, double(10.0));
    params.getParam<double>("scanSampleRate", scanSampleRate, double(1.0));
    params.getParam<int>("pointNum", pointNum, int(100));
    params.getParam<int>("trialNum", trialNum, int(1));

    std::cout << "\nParams:\n";
    params.write(std::cout);
    std::cout << std::endl;

    // Inserts few points to be matched
    std::normal_distribution<> noiseDistr(0.0, noiseSigma);
    std::uniform_real_distribution<> randPoint(-5.0, 5.0);
    for (int i = 0; i < pointNum; ++i) {
        psrc << randPoint(noiseGen), randPoint(noiseGen);
        pointsSrc.push_back(psrc);
    }
    glores::poseToIsometry2d(poseTrue, transfTrue);
    for (auto& p : pointsSrc) {
        pdst = transfTrue * p;
        pdst(0) += noiseDistr(noiseGen);
        pdst(1) += noiseDistr(noiseGen);
        pointsDst.push_back(pdst);
    }

    //    std::cout << "Point Src:\n";
    //    for (auto& p : pointsSrc) {
    //        std::cout << "  " << p.transpose() << std::endl;
    //    }
    //    std::cout << std::endl;
    //    std::cout << "Point Dst:\n";
    //    for (auto& p : pointsDst) {
    //        std::cout << "  " << p.transpose() << std::endl;
    //    }
    //    std::cout << std::endl;
    std::cout << "Generated points: pointsSrc " << pointsSrc.size() << ", pointsDst " << pointsDst.size() << std::endl;

    {
        glores::Transformation2d dstTsrc;
        std::vector<std::pair<int, int> > associations;
        double dist, func;

        std::cout << "Testing Error:\n";
        glores::poseToIsometry2d(poseTrue, dstTsrc);
        GLORES_PRINT_VARIABLE(poseTrue.transpose());
        std::cout << "transformation dstTsrc\n" << dstTsrc.matrix() << std::endl;
        std::cout << "found associations:\n";
        glores::associateNN(pointsSrc, pointsDst, dstTsrc, 1e+6, associations);

        func = 0.0;
        for (auto& a : associations) {
            dist = (pointsDst[a.second] - dstTsrc * pointsSrc[a.first]).norm();
            std::cout << "  " << a.first << "  " << a.second << ": dist " << dist << std::endl;
            func += dist * dist;
        }
        std::cout << "objective function " << func << std::endl;
    }


    // Optimizer
    optimizer.setPointSource(pointsSrc);
    optimizer.setPointDestination(pointsDst);
    optimizer.setXTolerance(optXTol);
    optimizer.setYTolerance(optYTol);
    optimizer.enableXTolerance(false);
    optimizer.enableYTolerance(true);
    GLORES_PRINT_VARIABLE(optInlierRatio);
    optimizer.setInlierRatio(optInlierRatio);
    optimizer.setRelaxationThres(optRelaxX, optRelaxY, optRelaxT);
    optimizer.setDepthMax(optDepthMax);
    optimizer.setIterationMax(optIterationMax);
    optimizer.enableRelaxationBound(optRelaxationBoundOn);
    optimizer.enableCheapBoundQueueOn(optCheaBoundQueueOn);
    optimizer.setVerboseLevel(optVerboseLevel);

    GLORES_PRINT_VARIABLE(poseTrue.transpose());
    GLORES_PRINT_VARIABLE(optimizer.evaluateFunction(poseTrue));
    GLORES_PRINT_VARIABLE(poseCoder.transpose());
    GLORES_PRINT_VARIABLE(optimizer.evaluateFunction(poseCoder));

    std::cout << "\n*********\nStarting optimization:\n";
    for (int k = 0; k < trialNum; ++k) {
        glores::ScopedTimer timer("findGlobalMin()");
        optimizer.findGlobalMin(poseLow, poseUpp, poseOpt, funcLow, funcUpp);
    }
    std::cout << "\n*****\nFound in " << poseOpt.transpose() << ": estimate between " << funcLow << " and " << funcUpp << std::endl;
    //std::cout << "  true value " << poseTrue.transpose() << std::endl;
    std::cout << "\nTimes:\n";
    glores::Profiler::getProfiler().printStats(std::cout);
    std::cout << std::endl;

    std::stringstream ss;
    std::string timeStatsFilename;
    ss << "time_stats_cheap_queue_" << (optCheaBoundQueueOn ? "ON_" : "OFF_") << "pointnum_" << pointNum << "_";
    timeStatsFilename = glores::generateStampedString(ss.str(), ".txt");
    std::ofstream timeStatsFile(timeStatsFilename);
    params.write(timeStatsFile, "# ");
    timeStatsFile << pointNum << "\t ";
    glores::Profiler::getProfiler().printStats(timeStatsFile);
    timeStatsFile.close();
    //    std::ofstream boundStatFile(glores::generateStampedString("stats_bound", ".txt"));
    //    params.write(boundStatFile, "# ");
    //    boundStatFile << "#iteration lower_bound upper_bound relaxation_bound_on queue_size\n";
    //    for (auto& bs : optimizer.getBoundStats()) {
    //        boundStatFile << bs.iteration << "\t " << bs.lowerBound << "\t " << bs.upperBound
    //                << "\t " << (bs.relaxationBoundOn ? 1 : 0) << "\t " << bs.queueNum << "\n";
    //    }
    //    boundStatFile.close();
    //
    //    std::ofstream scanFile(glores::generateStampedString("scan_aligned", ".plot"));
    //    printPointSets(scanFile, poseOpt, pointsSrc, pointsDst);
    //    scanFile.close();

    return 0;
}

void saturateRanges(glores::VectorPoint2d& points, double distMax) {
    glores::VectorPoint2d tmp;
    for (int i = 0; i < points.size(); ++i) {
        if (points[i].norm() < distMax) {
            tmp.push_back(points[i]);
        }
    }
    tmp.swap(points);
}

void samplePoints(glores::VectorPoint2d& points, double sampleRate) {
    glores::VectorPoint2d tmp;
    double sampleFrac, sampleInc;

    if (sampleRate <= 1e-6) {
        return;
    }
    sampleInc = 1.0 / sampleRate;
    //GLORES_PRINT_VARIABLE(sampleInc);

    sampleFrac = 0.0;
    for (int i = 0; i < points.size(); ++i) {
        if (sampleFrac >= 1.0) {
            tmp.push_back(points[i]);
            sampleFrac = sampleFrac - floor(sampleFrac);
            //GLORES_PRINT_MSG("sampling in i " << i << " with sampleFrac " << sampleFrac);
        }
        sampleFrac += sampleInc;
        //GLORES_PRINT_VARIABLE(sampleFrac);
    }
    points.swap(tmp);
}

void printPointSets(std::ostream& out, const Eigen::VectorXd & pose, glores::VectorPoint2d pointSrc, glores::VectorPoint2d pointDst) {
    glores::Transformation2d dstTsrc;

    glores::poseToIsometry2d(pose, dstTsrc);

    //out << "set term wxt 0\n";
    out << "set terminal postscript eps size 4.3,4.0 enhanced color font 'Arial,20'\n"
            << "set output 'scans_algined.eps'\n";
    out << "set size ratio -1\n";
    out << "set xlabel \"x\"\n";
    out << "set ylabel \"y\"\n";
    out << "set key outside\n";
    out << "plot "
            << "'-' u 1:2 title \"src\" w p pt 7 ps 1.5, "
            << "'-' u 1:2 title \"dst\" w p pt 4 ps 1.5, "
            << "'-' u 1:2 title \"rot\" w p pt 1 ps 1.1  "
            << "\n";
    for (auto& p : pointSrc) {
        out << p.transpose() << "\n";
    }
    out << "e\n";
    for (auto& p : pointDst) {
        out << p.transpose() << "\n";
    }
    out << "e\n";
    for (auto& p : pointSrc) {
        out << (dstTsrc * p).transpose() << "\n";
    }
    out << "e\n";
}


