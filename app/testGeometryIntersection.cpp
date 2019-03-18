/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   testProjectSE2Interval.cpp
 * Author: dario
 *
 * Created on 22 October 2018, 14:14
 */

#include <iostream>
#include <glores/geometry.h>
#include <glores/ParamMap.h>
#include <glores/thirdparty/gnuplot-iostream.h>

#include <boost/geometry.hpp>

#include "glores/fileutils.h"

#define DIFF_ANGLE_DEG_EXAMPLE(BEG,END) std::cout << "normalize((" << (END) << ") - (" << (BEG) << ")) = " \
          << ( 180.0/M_PI * glores::normalizeAngle2Pi(M_PI/180.0*((END) - (BEG))) ) << std::endl;


using namespace std;

void plotLine(std::ostream& out, double theta, double rho, double l = 1.0);

void plotArcCircle(std::ostream& out, double radius, double angleBeg, double angleEnd, double dtheta = M_PI / 180.0 * 0.2);

void plotCircle(std::ostream& out, double radius, const glores::Point2d& center, double dtheta = M_PI / 180.0 * 0.2);

void plotSegment(std::ostream& out, const glores::Point2d& p1, const glores::Point2d& p2);

int main(int argc, char** argv) {
    // General application parameters
    glores::ParamMap params;
    Eigen::Vector3d poseLow, poseUpp;
    glores::Point2d p, q, ptr;
    int vnum, sampleNum;
    glores::Polygon2d poly;
    glores::VectorPoint2d samples;
    glores::VectorPoint2d samplesRot;
    glores::VectorPoint2d samplesTra;
    glores::VectorPoint2d testPoints;
    double dmin, dmax;
    glores::Point2d pmin, pmax, nearestLine, nearestArc, seg1, seg2;
    double radius, rho, theta, angleIntersBeg, angleIntersEnd, arcBeg, arcEnd;
    bool plot0, plot1, plot2, plot3, plot4, plot5, plot6;
    std::string filenameCfg;
    //Gnuplot gp("gnuplot -persist");
    std::ostream& gp = std::cerr;

    params.read(argc, argv);
    params.getParam<std::string>("cfg", filenameCfg, std::string(""));
    if (!params.read(filenameCfg)) {
        std::cout << "Cannot open configuration file \"" << filenameCfg << "\": using default values" << std::endl;
    }

    // Reads file params from command line
    params.read(argc, argv);
    params.getParam<double>("px", p(0), double(1.0));
    params.getParam<double>("py", p(1), double(1.0));
    params.getParam<double>("qx", q(0), double(2.0));
    params.getParam<double>("qy", q(1), double(2.0));
    params.getParam<double>("txL", poseLow(0), double(0.0));
    params.getParam<double>("txU", poseUpp(0), double(1.0));
    params.getParam<double>("tyL", poseLow(1), double(-0.3));
    params.getParam<double>("tyU", poseUpp(1), double(+0.5));
    params.getParam<double>("rotL", poseLow(2), double(+10.0));
    params.getParam<double>("rotU", poseUpp(2), double(+36.0));
    poseLow(2) *= M_PI / 180.0;
    poseUpp(2) *= M_PI / 180.0;
    params.getParam<int>("vnum", vnum, int(5));
    params.getParam<int>("sampleNum", sampleNum, int(1000));
    params.getParam<double>("radius", radius, double(2.0));
    params.getParam<double>("rho", rho, double(1.0));
    params.getParam<double>("theta", theta, double(65.0));
    theta *= M_PI / 180.0;
    params.getParam<double>("arcBeg", arcBeg, double(95.0));
    params.getParam<double>("arcEnd", arcEnd, double(165.0));
    arcBeg *= M_PI / 180.0;
    arcEnd *= M_PI / 180.0;
    params.getParam<double>("seg1x", seg1(0), double(2.5));
    params.getParam<double>("seg1y", seg1(1), double(0.5));
    params.getParam<double>("seg2x", seg2(0), double(-0.5));
    params.getParam<double>("seg2y", seg2(1), double(2.8));
    params.getParam<bool>("plot0", plot0, bool(false));
    params.getParam<bool>("plot1", plot1, bool(true));
    params.getParam<bool>("plot2", plot2, bool(false));
    params.getParam<bool>("plot3", plot3, bool(false));
    params.getParam<bool>("plot4", plot4, bool(false));
    params.getParam<bool>("plot5", plot5, bool(false));
    params.getParam<bool>("plot6", plot6, bool(false));

    std::cout << "\nParams:\n";
    params.write(std::cout);
    std::cout << std::endl;


    if (plot0) {
        std::cout << "\n-----------------\n";
        glores::projectSE2Interval(p, poseLow, poseUpp, poly, vnum);
        std::cout << "\nBounding polygon for transformation of p " << p.transpose() << ": "
                << "\n  lower " << poseLow.transpose()
                << "\n  upper " << poseUpp.transpose()
                << "\n";
        //    for (auto pit = boost::geometry::begin_points(poly); pit != boost::geometry::end_points(poly); ++pit) {
        //        std::cout << "  " << pit->transpose() << "\n";
        //    }
        for (auto pit = boost::begin(boost::geometry::exterior_ring(poly)); pit != boost::end(boost::geometry::exterior_ring(poly)); ++pit) {
            //        double x = boost::geometry::get<0>(*pit);
            //        double y = boost::geometry::get<1>(*pit);
            //        std::cout << "  " << x << " " << y << "\n";
            std::cout << "  " << pit->transpose() << "\n";
        }

        // Computes the random transformation of point p over the given transformation
        // interval [poseLow, poseUpp]. 
        samples.clear();
        for (int i = 0; i < sampleNum; ++i) {
            double tx = poseLow(0) + (poseUpp(0) - poseLow(0)) * rand() / RAND_MAX;
            double ty = poseLow(1) + (poseUpp(1) - poseLow(1)) * rand() / RAND_MAX;
            double theta = poseLow(2) + glores::normalizeAngle2Pi(poseUpp(2) - poseLow(2)) * rand() / RAND_MAX;
            if (tx < poseLow(0) || tx > poseUpp(0) || ty < poseLow(1) || ty > poseUpp(1)) {
                std::cout << " translation out of bound: tx " << tx << ", ty " << ty << std::endl;
            }
            glores::Transformation2d tr = glores::Transformation2d::Identity();
            tr.prerotate(Eigen::Rotation2Dd(theta));
            tr.pretranslate(Eigen::Vector2d(tx, ty));
            //std::cout << "transformation\n" << tr.matrix() << std::endl;
            samples.push_back(tr * p);
        }


        gp << "set term wxt 0\n";
        gp << "set size ratio -1\n";
        gp << "set xlabel \"x\"\n";
        gp << "set ylabel \"y\"\n";
        gp << "set key outside\n";
        gp << "plot "
                << "'-' u 1:2 title \"samples\" w p pt 7 ps 0.5"
                << ", '-' u 1:2 title \"bound\" w l lt 1"
                << "\n";
        for (auto& sp : samples) {
            gp << sp(0) << " " << sp(1) << "\n";
        }
        gp << "e\n";
        for (auto pit = boost::begin(boost::geometry::exterior_ring(poly)); pit != boost::end(boost::geometry::exterior_ring(poly)); ++pit) {
            double cx = boost::geometry::get<0>(*pit);
            double cy = boost::geometry::get<1>(*pit);
            gp << "  " << cx << " " << cy << "\n";
        }
        gp << "e\n";
    }


    DIFF_ANGLE_DEG_EXAMPLE(-10.0, 160);
    DIFF_ANGLE_DEG_EXAMPLE(160.0, -10.0);
    DIFF_ANGLE_DEG_EXAMPLE(280.0, -15.0);
    DIFF_ANGLE_DEG_EXAMPLE(-15.0, 280.0);

    if (plot1) {

        // Computes the random transformation of point p over the given transformation
        // interval [poseLow, poseUpp]. 
        std::cout << "\n-----------------\n ARC-RECT\n";
        for (int i = 0; i < sampleNum; ++i) {
            double tx = poseLow(0) + (poseUpp(0) - poseLow(0)) * rand() / RAND_MAX;
            double ty = poseLow(1) + (poseUpp(1) - poseLow(1)) * rand() / RAND_MAX;
            double theta = poseLow(2) + glores::normalizeAngle2Pi(poseUpp(2) - poseLow(2)) * rand() / RAND_MAX;
            if (tx < poseLow(0) || tx > poseUpp(0) || ty < poseLow(1) || ty > poseUpp(1)) {
                std::cout << " translation out of bound: tx " << tx << ", ty " << ty << std::endl;
            }
            glores::Transformation2d rot = glores::Transformation2d::Identity();
            rot.prerotate(Eigen::Rotation2Dd(theta));
            samplesRot.push_back(rot * p);
            glores::Transformation2d tra = glores::Transformation2d::Identity();
            tra.pretranslate(Eigen::Vector2d(-tx, -ty));
            samplesTra.push_back(tra * q);
        }
        glores::computeDistanceMinMax2(p, q, poseLow, poseUpp, dmin, dmax, pmin, pmax, testPoints);
        std::cout << "\ndistance min " << dmin << " max " << dmax << "\n"
                << "  pmin " << pmin.transpose() << "\n"
                << "  pmax " << pmax.transpose() << "\n"
                << "  lower " << poseLow.transpose() << "\n"
                << "  upper " << poseUpp.transpose() << "\n"
                << std::endl;
        std::cout << "found " << testPoints.size() << " test points for distances:\n";
        for (auto& tp : testPoints) {
            std::cout << "  " << tp.transpose() << "\n";
        }
        double dmin3, dmax3;
        glores::computeDistanceMinMax3(p, q, poseLow, poseUpp, dmin3, dmax3);
        GLORES_PRINT_MSG("computeDistanceMinMax3(): dmin3 " << dmin3 << " vs " << dmin << ", "
                         "dmax3 " << dmax3 << " vs " << dmax);
        
        double bbXL = q(0) - poseUpp(0);
        double bbXU = q(0) - poseLow(0);
        double bbYL = q(1) - poseUpp(1);
        double bbYU = q(1) - poseLow(1);

        double dist = (poseUpp.head<2>() - poseLow.head<2>()).norm() + q.norm();

        double radiusPlot = p.norm();
        double betaPlot = atan2(p.y(), p.x());
        double arcBegPlot = betaPlot + poseLow(2);
        double arcEndPlot = betaPlot + poseUpp(2);

        glores::Point2d v1, v2, v3, v4;
        v1 << bbXL, bbYL;
        v2 << bbXU, bbYL;
        v3 << bbXU, bbYU;
        v4 << bbXL, bbYU;

        gp << "set term wxt 1\n";
        gp << "set size ratio -1\n";
        gp << "set title \"" << filenameCfg << " dmin " << std::fixed << std::setprecision(4) << dmin << " dmax " << std::setprecision(4) << dmax << "\"\n";
        gp << "set xlabel \"x\"\n";
        gp << "set ylabel \"y\"\n";
        gp << "set key outside\n";
        gp << "plot "
//                << "  '-' u 1:2 title \"samplesTra\" w p pt 7 ps 0.1"
//                << ", '-' u 1:2 title \"samplesRot\" w p pt 7 ps 0.1"
                << " '-' u 1:2 title \"lines\" w l lt 3 lw 1.2"
                << ", '-' u 1:2 title \"arc\" w l lt 3 lw 1.2"
                << ", '-' u 1:2 title \"q\" w p pt 7 ps 1.0"
//                << ", '-' u 1:2 title \"pmin\" w p pt 7 ps 1.0"
//                << ", '-' u 1:2 title \"pmax\" w p pt 7 ps 1.0"
                //<< ", '-' u 1:2 title \"test points\" w p pt 2 ps 1.0"
                << "\n";
//        for (auto& sp : samplesRot) {
//            gp << sp(0) << " " << sp(1) << "\n";
//        }
//        gp << "e\n";
//        for (auto& sp : samplesTra) {
//            gp << sp(0) << " " << sp(1) << "\n";
//        }
//        gp << "e\n";
        //        plotLine(gp, 0.0, bbXL, dist);
        //        plotLine(gp, M_PI, -bbXU, dist);
        //        plotLine(gp, 0.5 * M_PI, bbYL, dist);
        //        plotLine(gp, 1.5 * M_PI, -bbYU, dist);
        plotSegment(gp, v1, v2);
        plotSegment(gp, v2, v3);
        plotSegment(gp, v3, v4);
        plotSegment(gp, v4, v1);
        gp << "e\n";
        plotArcCircle(gp, radiusPlot, arcBegPlot, arcEndPlot);
        gp << "e\n";
        gp << q(0) << " " << q(1) << "\n";
        gp << "e\n";
//        gp << pmin.transpose() << "\ne\n" << pmax.transpose() << "\ne\n";
//        for (auto& tp : testPoints) {
//            gp << tp.transpose() << "\n";
//        }
//        gp << "e\n";
        //    for (auto pit = boost::begin(boost::geometry::exterior_ring(poly)); pit != boost::end(boost::geometry::exterior_ring(poly)); ++pit) {
        //        double cx = boost::geometry::get<0>(*pit);
        //        double cy = boost::geometry::get<1>(*pit);
        //        gp << "  " << cx << " " << cy << "\n";
        //    }
        //    gp << "e\n";
    }

    if (plot2) {

        std::cout << "\n-----------------\n";
        bool hplaneCircleInters = glores::intersectCircleHalfPlane(theta, rho, radius, angleIntersBeg, angleIntersEnd, nearestLine, nearestArc);
        std::cout << "intersection half-plane: x * " << cos(theta) << " + y * " << sin(theta) << " >= " << rho
                << " and circle with radius " << radius << ": intersect? " << hplaneCircleInters << "\n";
        if (hplaneCircleInters) {
            std::cout << "  angle interval: beg " << (180.0 / M_PI * angleIntersBeg)
                    << " end " << (180.0 / M_PI * angleIntersEnd) << " [deg]" << std::endl;
        }


        gp << "set term wxt 2\n";
        gp << "set size ratio -1\n";
        gp << "set xlabel \"x\"\n";
        gp << "set ylabel \"y\"\n";
        gp << "set key outside\n";
        gp << "plot "
                << "'-' u 1:2 title \"line\" w l, "
                << "'-' u 1:2 title \"circle\" w l, ";
        if (hplaneCircleInters) {
            gp << "'-' u 1:2 title \"inters\" w l";
        } else {
            gp << "'-' u 1:2 title \"nearest\" w p pt 7 ps 0.8";
        }
        gp << "\n";
        plotLine(gp, theta, rho, 2.0);
        gp << "e\n";
        plotArcCircle(gp, radius, 0.0, 2.0 * M_PI - 1e-6);
        gp << "e\n";
        if (hplaneCircleInters) {
            plotArcCircle(gp, radius, angleIntersBeg, angleIntersEnd);
        } else {
            gp << nearestLine(0) << " " << nearestLine(1) << "\n"
                    << nearestArc(0) << " " << nearestArc(1) << "\n";
        }
        gp << "e\n";
    }


    if (plot3) {
        std::cout << "\n-----------------\n";
        std::vector<glores::AngleInterval> intersectionArcs;
        bool hplaneCircularArcInters = glores::intersectCircularArcHalfPlane(theta, rho, radius, arcBeg, arcEnd, nearestLine, nearestArc, intersectionArcs);
        std::cout << "intersection? " << hplaneCircularArcInters << "\n"
                << "  half-plane: x * " << cos(theta) << " + y * " << sin(theta) << " - (" << rho << ") >= 0.0\n"
                << "  circular arc: radius " << radius << ", angle interval [" << (180.0 / M_PI * arcBeg) << " : " << (180.0 / M_PI * arcEnd) << "] deg\n";
        if (hplaneCircularArcInters) {
            std::cout << "  intersection occurs on arcs:\n";
            for (auto& arc : intersectionArcs) {
                std::cout << "    [" << (180.0 / M_PI * arc.first) << " : " << (180.0 / M_PI * arc.second) << "] deg\n";
            }
        } else {
            std::cout << "  NO intersection: nearestLine " << nearestLine.transpose() << ", nearestArc " << nearestArc.transpose() << std::endl;
        }


        gp << "set term wxt 3\n";
        gp << "set size ratio -1\n";
        gp << "set xlabel \"x\"\n";
        gp << "set ylabel \"y\"\n";
        gp << "set key outside\n";
        gp << "plot "
                << "'-' u 1:2 title \"line\" w l, "
                << "'-' u 1:2 title \"arc\" w l lw 1.0, ";
        if (hplaneCircularArcInters) {
            gp << "'-' u 1:2 title \"inters\" w l lw 0.3";
        } else {
            gp << "'-' u 1:2 title \"nearest\" w p pt 7 ps 0.8";
        }
        gp << "\n";
        plotLine(gp, theta, rho, (M_PI / 180.0 * 70.0));
        gp << "e\n";
        plotArcCircle(gp, radius, arcBeg, arcEnd);
        gp << "e\n";
        if (hplaneCircularArcInters) {
            for (auto& arc : intersectionArcs) {
                plotArcCircle(gp, radius, arc.first, arc.second);
            }
        } else {
            gp << nearestLine(0) << " " << nearestLine(1) << "\n"
                    << nearestArc(0) << " " << nearestArc(1) << "\n";
            gp << "e\n";
        }
    }

    if (plot4) {
        std::cout << "\n-----------------\n";
        std::vector<glores::AngleInterval> intersectionArcs;
        std::vector<glores::DistancePointPair> candidatePoints;
        bool hplaneCircularArcInters = glores::intersectCircularArcHalfPlaneSegment(seg1, seg2, radius, arcBeg, arcEnd, nearestLine, nearestArc, intersectionArcs, candidatePoints);
        std::cout << "intersection? " << hplaneCircularArcInters << "\n"
                << "  segment [" << seg1.transpose() << "] - [" << seg2.transpose() << "]\n"
                << "  circular arc: radius " << radius << ", angle interval [" << (180.0 / M_PI * arcBeg) << " : " << (180.0 / M_PI * arcEnd) << "] deg\n";
        if (hplaneCircularArcInters) {
            std::cout << "  intersection occurs on arcs:\n";
            for (auto& arc : intersectionArcs) {
                std::cout << "    [" << (180.0 / M_PI * arc.first) << " : " << (180.0 / M_PI * arc.second) << "] deg\n";
            }
        } else {
            std::cout << "  NO intersection: nearestLine " << nearestLine.transpose() << ", nearestArc " << nearestArc.transpose() << std::endl;
        }


        gp << "set term wxt 4\n";
        gp << "set size ratio -1\n";
        gp << "set xlabel \"x\"\n";
        gp << "set ylabel \"y\"\n";
        gp << "set key outside\n";
        gp << "plot "
                << "'-' u 1:2 title \"segment\" w l, "
                << "'-' u 1:2 title \"seg1\" w p pt 2 ps 1.6, "
                << "'-' u 1:2 title \"seg2\" w p pt 3 ps 1.6, "
                << "'-' u 1:2 title \"arc\" w l lw 1.0, ";
        if (hplaneCircularArcInters) {
            gp << "'-' u 1:2 title \"inters\" w l lw 0.3";
        } else {
            gp << "'-' u 1:2 title \"candidates\" w l lw 0.8, "
                    << "'-' u 1:2 title \"nearest\" w p pt 7 ps 0.8";
        }
        gp << "\n";
        //        gp << seg1.transpose() << "\n"
        //           << seg2.transpose() << "\n";
        plotSegment(gp, seg1, seg2);
        gp << "e\n";
        gp << seg1.transpose() << "\n";
        gp << "e\n";
        gp << seg2.transpose() << "\n";
        gp << "e\n";
        plotArcCircle(gp, radius, arcBeg, arcEnd);
        gp << "e\n";
        if (hplaneCircularArcInters) {
            for (auto& arc : intersectionArcs) {
                plotArcCircle(gp, radius, arc.first, arc.second);
            }
        } else {
            for (auto& c : candidatePoints) {
                gp << std::get<1>(c) (0) << " " << std::get<1>(c) (1) << "\n"
                        << std::get<2>(c) (0) << " " << std::get<2>(c) (1) << "\n\n";
            }
            gp << "e\n";
            gp << nearestLine(0) << " " << nearestLine(1) << "\n"
                    << nearestArc(0) << " " << nearestArc(1) << "\n";
            gp << "e\n";
        }
    }



    if (plot5) {
        std::cout << "\n-----------------\n";
        glores::Point2d farthest;
        double distMax;

        distMax = glores::computeMaxDistancePointArc(p, radius, arcBeg, arcEnd, farthest);

        std::cout << "max distance point-arc " << distMax << "\n";

        gp << "set term wxt 5\n";
        gp << "set size ratio -1\n";
        gp << "set title \"max dist arc(radius " << std::setprecision(2) << std::fixed << radius << ", "
                << (180.0 / M_PI * glores::normalizeAngle2Pi(arcBeg)) << ":" << (180.0 / M_PI * glores::normalizeAngle2Pi(arcEnd))
                << ") - point [" << p.transpose() << "]\"\n";
        gp << "set xlabel \"x\"\n";
        gp << "set ylabel \"y\"\n";
        gp << "set key outside\n";
        gp << "plot "
                << "'-' u 1:2 title \"arc\" w l, "
                << "'-' u 1:2 title \"distance\" w l lt 2 lw 0.8,"
                << "'-' u 1:2 title \"p\" w p pt 7 ps 1.6, "
                << "'-' u 1:2 title \"farthest\" w p pt 7 ps 1.3"
                << "\n";
        plotArcCircle(gp, radius, arcBeg, arcEnd);
        gp << "e\n";
        plotCircle(gp, distMax, p);
        gp << "e\n";
        gp << p.transpose() << "\n";
        gp << "e\n";
        gp << farthest.transpose() << "\n";
        gp << "e\n";
    }

    if (plot6) {
        std::cout << "\n-----------------\n";
        glores::VectorPoint2d vertices;
        glores::projectArcBound(arcBeg, arcEnd, vertices, vnum);
        
        std::cout << "arc [" << (180.0 / M_PI * glores::normalizeAngle2Pi(arcBeg)) << ":" << (180.0 / M_PI * glores::normalizeAngle2Pi(arcEnd)) << "] deg"
                << " bounded by polygon:\n";
        for (auto& v : vertices) {
            std::cout << "  " << v(0) << " " << v(1) << "\n";
        }
        std::cout << "\n";

        gp << "set term wxt 5\n";
        gp << "set size ratio -1\n";
        gp << "set title \"bounding arc ("
                << (180.0 / M_PI * glores::normalizeAngle2Pi(arcBeg)) << ":" << (180.0 / M_PI * glores::normalizeAngle2Pi(arcEnd))
                << ")\"\n";
        gp << "set xlabel \"x\"\n";
        gp << "set ylabel \"y\"\n";
        gp << "set key outside\n";
        gp << "plot "
                << "'-' u 1:2 title \"arc\" w l, "
                << "'-' u 1:2 title \"bounding\" w l"
                << "\n";
        plotArcCircle(gp, 1.0, arcBeg, arcEnd);
        gp << "e\n";
        for (int i = 0; i < vertices.size(); ++i) {
            gp << "  " << vertices[i](0) << " " << vertices[i](1) << "\n";
        }
        gp << "  " << vertices[0](0) << " " << vertices[0](1) << "\n";
        gp << "e\n";
    }

    return 0;
}

void plotLine(std::ostream& out, double theta, double rho, double l) {
    double ct, st, x1, y1, x2, y2, x0, y0, xt, yt, dir;
    ct = cos(theta);
    st = sin(theta);

    dir = 0.05 * fabs(l);
    //    if (fabs(rho) < 1e-3) {
    //        l = 1.0;
    //        dir = 0.1;
    //    } else {
    //        l = rho * tan(delta);
    //        dir = 0.1 * fabs(rho);
    //    }
    x1 = rho * cos(theta) - l * sin(theta);
    y1 = rho * sin(theta) + l * cos(theta);
    x2 = rho * cos(theta) + l * sin(theta);
    y2 = rho * sin(theta) - l * cos(theta);
    x0 = rho * cos(theta);
    y0 = rho * sin(theta);
    xt = x0 + dir * cos(theta);
    yt = y0 + dir * sin(theta);


    //    GLORES_PRINT_VARIABLE(rho);
    //    GLORES_PRINT_VARIABLE(180.0 / M_PI * delta);
    //    GLORES_PRINT_VARIABLE(rho * tan(delta));

    out << x1 << " " << y1 << "\n" << x2 << " " << y2 << "\n\n";
    out << x0 << " " << y0 << "\n" << xt << " " << yt << "\n\n";
}

void plotSegment(std::ostream& out, const glores::Point2d& p1, const glores::Point2d& p2) {
    glores::Point2d mid, ver;
    double theta, rho, len;

    glores::computeLine(p1, p2, theta, rho);
    len = (p2 - p1).norm();
    mid = 0.5 * (p1 + p2);
    ver << mid(0) + 0.08 * len * cos(theta), mid(1) + 0.08 * len * sin(theta);

    out << p1(0) << " " << p1(1) << "\n" << p2(0) << " " << p2(1) << "\n\n";
    out << mid(0) << " " << mid(1) << "\n" << ver(0) << " " << ver(1) << "\n\n";
}

void plotArcCircle(std::ostream& out, double radius, double angleBeg, double angleEnd, double dtheta) {
    double x1, y1, x2, y2;
    int n;

    n = floor(glores::normalizeAngle2Pi(angleEnd - angleBeg) / dtheta);
    x1 = radius * cos(angleBeg);
    y1 = radius * sin(angleBeg);
    for (int i = 1; i <= n; ++i) {
        x2 = radius * cos(angleBeg + i * dtheta);
        y2 = radius * sin(angleBeg + i * dtheta);
        out << x1 << " " << y1 << "\n" << x2 << " " << y2 << "\n\n";
        x1 = x2;
        y1 = y2;
    }
}

void plotCircle(std::ostream& out, double radius, const glores::Point2d& center, double dtheta) {
    double x1, y1, x2, y2;
    int n;

    n = floor(2 * M_PI / dtheta);
    x1 = center(0) + radius;
    y1 = center(1);
    for (int i = 1; i <= n; ++i) {
        x2 = center(0) + radius * cos(dtheta * i);
        y2 = center(1) + radius * sin(dtheta * i);
        out << x1 << " " << y1 << "\n" << x2 << " " << y2 << "\n\n";
        x1 = x2;
        y1 = y2;
    }
}
