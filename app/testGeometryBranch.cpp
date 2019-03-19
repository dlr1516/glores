/**
 * GLORES - Globally Optimal Registration based on Fast Branch and Bound
 * Copyright (C) 2019 Luca Consolini, Mattia Laurini, Marco Locatelli, Dario Lodi Rizzini.
 *
 * GLORES is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * GLORES is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with GLORES.  If not, see <http://www.gnu.org/licenses/>.
 */
#include <iostream>
#include <deque>
#include <Eigen/Dense>

#include <glores/geometry.h>
//#include <glores/PointSetGlobalRegistration.h>
#include <glores/ParamMap.h>
#include <glores/thirdparty/gnuplot-iostream.h>

#include <boost/geometry.hpp>

#include "glores/fileutils.h"

#define DIFF_ANGLE_DEG_EXAMPLE(BEG,END) std::cout << "normalize((" << (END) << ") - (" << (BEG) << ")) = " \
          << ( 180.0/M_PI * glores::normalizeAngle2Pi(M_PI/180.0*((END) - (BEG))) ) << std::endl;


using namespace std;

struct Interval {
    int id;
    int parent;
    int depth;
    Eigen::Vector3d poseLow;
    Eigen::Vector3d poseUpp;
    double dmin;
    double dmax;
};

struct PointSort {
    double dmin;
    double dmax;
    int id;

    PointSort(int id_, double dmin_, double dmax_) : id(id_), dmin(dmin_), dmax(dmax_) {
    }
};

bool operator<(const PointSort& ps1, const PointSort& ps2) {
    return (ps1.dmin < ps2.dmin);
}

void plotLine(std::ostream& out, double theta, double rho, double delta = (M_PI / 180.0 * 45.0));

void plotSegment(std::ostream& out, const glores::Point2d& p1, const glores::Point2d& p2);

void plotArcCircle(std::ostream& out, double radius, double angleBeg, double angleEnd, double dtheta = M_PI / 180.0 * 0.2);

void branchInterval(const Interval& father, Interval& childL, Interval& childR, int dim);

std::ostream& operator<<(std::ostream& out, const Interval& interval);

int main(int argc, char** argv) {
    // General application parameters
    glores::ParamMap params;
    std::string filenameCfg;
    Eigen::Vector3d poseLow, poseUpp;
    glores::Point2d p, q, ptr;
    glores::Point2d pSrc;
    glores::Point2d pDst;
    glores::VectorPoint2d pDsts;

    int vnum, sampleNum;
    glores::VectorPoint2d testPoints;
    double dmin, dmax;
    glores::Point2d pmin, pmax;
    int depthMax;
    bool plot;
    Interval curr, left, right;
    int icurr, v, par, plotLeaf;
    std::set<PointSort> psorted;
    //Gnuplot gp("gnuplot -persist");
    std::ostream& gp = std::cerr;


    params.read(argc, argv);
    params.getParam<std::string>("cfg", filenameCfg, std::string(""));
    if (!params.read(filenameCfg)) {
        std::cout << "Cannot open configuration file \"" << filenameCfg << "\": using default values" << std::endl;
    }

    // Reads file params from command line
    params.read(argc, argv);
    params.getParam<double>("px", pSrc(0), double(1.0));
    params.getParam<double>("py", pSrc(1), double(1.0));
    params.getParam<double>("qx", q(0), double(2.0));
    params.getParam<double>("qy", q(1), double(2.0));
    params.getParam<double>("txL", poseLow(0), double(-2.0));
    params.getParam<double>("txU", poseUpp(0), double(+2.0));
    params.getParam<double>("tyL", poseLow(1), double(-2.0));
    params.getParam<double>("tyU", poseUpp(1), double(+2.0));
    params.getParam<double>("rotL", poseLow(2), double(45.0));
    params.getParam<double>("rotU", poseUpp(2), double(90.0));
    poseLow(2) *= M_PI / 180.0;
    poseUpp(2) *= M_PI / 180.0;
    params.getParam<int>("vnum", vnum, int(10));
    params.getParam<int>("sampleNum", sampleNum, int(1000));
    params.getParam<int>("depthMax", depthMax, int(3));
    params.getParam<bool>("plot", plot, bool(true));
    params.getParam<int>("plotLeaf", plotLeaf, int(8));


    std::cout << "\nParams:\n";
    params.write(std::cout);
    std::cout << std::endl;

    std::vector<Interval> intervals;
    std::vector<int> queue;
    std::vector<int> leaves;
    curr.id = 0;
    curr.parent = -1;
    curr.depth = 0;
    curr.poseLow = poseLow;
    curr.poseUpp = poseUpp;
    intervals.push_back(curr);
    queue.push_back(curr.id);

    while (!queue.empty()) {
        icurr = queue.back();
        queue.pop_back();

        glores::computeDistanceMinMax2(p, q, intervals[icurr].poseLow, intervals[icurr].poseUpp, intervals[icurr].dmin, intervals[icurr].dmax, pmin, pmax, testPoints);

        //        for (int j = 0; j < intervals[icurr].depth; ++j) {
        //            std::cout << "  ";
        //        }
        //        std::cout << intervals[icurr] << "\n";

        if (intervals[icurr].depth <= depthMax) {
            branchInterval(intervals[icurr], left, right, ((intervals[icurr].depth + 0) % 3));
            left.id = intervals.size();
            intervals.push_back(left);
            right.id = intervals.size();
            intervals.push_back(right);
            queue.push_back(left.id);
            queue.push_back(right.id);
        } else {
            leaves.push_back(icurr);
        }
    }

    //    std::cout << "\nTesting father-children bounds:\n";
    for (auto& l : leaves) {
        //        std::cout << "leaf " << l << "\n";
        //        std::cout << "  dmin: ";
        v = l;
        while (0 <= v && v < intervals.size()) {
            //            std::cout << intervals[v].dmin << "[" << intervals[v].id << "]";
            par = intervals[v].parent;
            if (par >= 0) {
                //                std::cout << " >= ";
                if (intervals[v].dmin < intervals[par].dmin) {
                    std::cerr << "\n*** ERROR: " << intervals[v].dmin << "[" << intervals[v].id << "] < " << intervals[par].dmin << "[" << intervals[par].id << "]" << std::endl;
                    exit(-1);
                }
            }
            v = intervals[v].parent;
        }
        std::cout << std::endl;
        //        std::cout << "  dmax: ";
        v = l;
        while (0 <= v && v < intervals.size()) {
            //            std::cout << intervals[v].dmax << "[" << intervals[v].id << "]";
            par = intervals[v].parent;
            if (par >= 0) {
                //                std::cout << " <= ";
                if (intervals[v].dmax > intervals[par].dmax) {
                    std::cerr << "\n*** ERROR: " << intervals[v].dmax << "[" << intervals[v].id << "] > " << intervals[par].dmax << "[" << intervals[par].id << "]" << std::endl;
                    exit(-1);
                }
            }
            v = intervals[v].parent;
        }
        //        std::cout << std::endl;
    }


    std::cout << "***\nTest: on the same pose interval I, given src point and dst_j,\n"
            << "  dmin(src, dst*) < dmin(src, dst_j) => dmax(src, dst*) < dmax(src, dst_j)?\n";

    if (plotLeaf < 0 || plotLeaf >= intervals.size()) {
        poseLow(0) = -2.0;
        poseUpp(0) = +1.0;
        poseLow(1) = 0.5;
        poseUpp(1) = 1.8;
        poseLow(2) = M_PI / 180.0 * (70.0);
        poseUpp(2) = M_PI / 180.0 * (350.0);
    } else {
        poseLow = intervals[plotLeaf].poseLow;
        poseUpp = intervals[plotLeaf].poseUpp;
    }


    //pSrc << 1.0, 1.0;
    double L = 5.0;
    unsigned int seed = time(0);
    //    std::cout << "*** SEED " << seed << std::endl;
    //    srand(1547110537);
    //    //srand(100);
    //    for (int i = 0; i < vnum; ++i) {
    //        pDst(0) = pSrc(0) - L + 2.0 * L * rand() / RAND_MAX;
    //        pDst(1) = pSrc(1) - L + 2.0 * L * rand() / RAND_MAX;
    //        pDsts.push_back(pDst);
    //    }

    pDsts.push_back(glores::Point2d(5.8621, 5.47228));
    pDsts.push_back(glores::Point2d(4.99387, 4.67025));
    pDsts.push_back(glores::Point2d(2.07078, -3.16907));
    pDsts.push_back(glores::Point2d(4.58242, 3.12266));
    pDsts.push_back(glores::Point2d(-1.85663, 1.69991));
    pDsts.push_back(glores::Point2d(0.674517, -3.79928));
    pDsts.push_back(glores::Point2d(-2.60072, 0.284387));
    pDsts.push_back(glores::Point2d(5.25282, -2.86258));
    pDsts.push_back(glores::Point2d(1.38624, 2.17236));
    pDsts.push_back(glores::Point2d(2.12282, 4.13161));

    // Checks the size of the queue 
    for (auto& interv : intervals) {
        // Finds the minimum of dmax
        double dmaxMin = 1e+6;
        int imaxMin = -1;
        for (int j = 0; j < pDsts.size(); ++j) {
            glores::computeDistanceMinMax2(pSrc, pDsts[j], interv.poseLow, interv.poseUpp, dmin, dmax, pmin, pmax, testPoints);
            if (dmax < dmaxMin) {
                dmaxMin = dmax;
                imaxMin = j;
            }
        }
        //
        double dminMin = 1e+6;
        int iminMin = -1;
        int countQueue = 0;
        for (int j = 0; j < pDsts.size(); ++j) {
            glores::computeDistanceMinMax2(pSrc, pDsts[j], interv.poseLow, interv.poseUpp, dmin, dmax, pmin, pmax, testPoints);
            if (dmin < dminMin) {
                dminMin = dmin;
                iminMin = j;
            }
            if (dmin < dmaxMin) {
                countQueue++;
            }
        }
        //
        for (int j = 0; j < interv.depth; ++j) {
            std::cout << "  ";
        }
        std::cout << interv.id << ": "
                << "dmin " << dminMin << " [" << iminMin << "], "
                << "dmax " << dmaxMin << " [" << imaxMin << "], "
                << "items with dmin(i) < dmax " << countQueue << " "
                << " par " << interv.parent << " "
                << "low [" << interv.poseLow.transpose() << "] "
                << "upp [" << interv.poseUpp.transpose() << "] "
                << std::endl;
    }


    if (plot) {
        double radiusPlot = pSrc.norm();
        double betaPlot = atan2(pSrc.y(), pSrc.x());
        double arcBegPlot = betaPlot + poseLow(2);
        double arcEndPlot = betaPlot + poseUpp(2);

        params.write(gp, "# ");
        gp << "#\n"
                << "# poseLow \t" << poseLow.transpose() << "\n"
                << "# poseUpp \t" << poseUpp.transpose() << "\n"
                << "# Points:\n";
        int count = 0;
        for (auto& pTmp : pDsts) {
            glores::computeDistanceMinMax2(pSrc, pTmp, poseLow, poseUpp, dmin, dmax, pmin, pmax, testPoints);
            gp << "# " << count << " " << pTmp.transpose() << "\t dmin " << dmin << "\t dmax " << dmax << "\n";
            count++;
        }

        gp << "set term wxt 1\n";
        gp << "set size ratio -1\n";
        gp << "set xlabel \"x\"\n";
        gp << "set ylabel \"y\"\n";
        gp << "set key outside\n";
        gp << "plot "
                << " '-' u 1:2 title \"src arc\" w l";
        count = 0;
        for (auto& pTmp : pDsts) {
            gp << ", '-' u 1:2 title \"box " << count << "\" w l"
                    << ", '-' u 1:2 title \"pmin/pmax " << count << "\" w p pt 7";
            count++;
        }
        gp << "\n";

        plotArcCircle(gp, radiusPlot, arcBegPlot, arcEndPlot);
        gp << "e\n";
        count = 0;

        std::cout << "Box: poseLow " << poseLow.transpose() << ", poseUpp " << poseUpp.transpose() << std::endl;
        for (auto& pTmp : pDsts) {
            glores::computeDistanceMinMax2(pSrc, pTmp, poseLow, poseUpp, dmin, dmax, pmin, pmax, testPoints);

            std::cout << " src [" << pSrc.transpose() << "] dst " << count << " [" << pTmp.transpose() << "]: dmin " << dmin << " dmax " << dmax << std::endl;

            psorted.insert(PointSort(count, dmin, dmax));

            double bbXL = pTmp(0) - poseUpp(0);
            double bbXU = pTmp(0) - poseLow(0);
            double bbYL = pTmp(1) - poseUpp(1);
            double bbYU = pTmp(1) - poseLow(1);
            glores::Point2d v1, v2, v3, v4;
            v1 << bbXL, bbYL;
            v2 << bbXU, bbYL;
            v3 << bbXU, bbYU;
            v4 << bbXL, bbYU;

            plotSegment(gp, v1, v2);
            plotSegment(gp, v2, v3);
            plotSegment(gp, v3, v4);
            plotSegment(gp, v4, v1);
            gp << "e\n";
            gp << pmin.transpose() << "\n"
                    << pmax.transpose() << "\n"
                    << "e\n";
            count++;
        }

        std::cout << "\nPoint sorted by distance:\n";
        for (auto& ps : psorted) {
            std::cout << "  " << ps.id
                    << " dmin " << std::fixed << std::setprecision(3) << ps.dmin
                    << " dmax " << std::fixed << std::setprecision(3) << ps.dmax << "\n";
        }

        std::cout << "\n\n\\begin{tabular}{|l|";
        for (auto& ps : psorted) {
            std::cout << "c|";
        }
        std::cout << "}\n  \\hline\n";
        std::cout << "  $i$";
        for (auto& ps : psorted) {
            std::cout << " & " << ps.id;
        }
        std::cout << " \\\\\n";
        std::cout << "  $d_{min,i}$";
        for (auto& ps : psorted) {
            std::cout << " & " << ps.dmin;
        }
        std::cout << " \\\\\n";
        std::cout << "  \\hline\n\\end{tabular}\n";


        if (!psorted.empty()) {
            int iMaxMin = -1;
            double dMaxMin = 1e+6;
            for (auto& ps : psorted) {
                if (ps.dmax < dMaxMin) {
                    iMaxMin = ps.id;
                    dMaxMin = ps.dmax;
                }
            }
            std::cout << "\n\n\\begin{tabular}{|l|c|}\n  \\hline\n";
            std::cout << "  $i_{max}$ & " << iMaxMin << "\\\\\n";
            std::cout << "  $d_{max}$ & " << dMaxMin << "\\\\\n";
            std::cout << "  \\hline\n\\end{tabular}\n";
        }


        //                << "  '-' u 1:2 title \"samplesRot\" w p pt 7 ps 0.1"
        //                << ", '-' u 1:2 title \"samplesTra\" w p pt 7 ps 0.1"
        //                << ", '-' u 1:2 title \"pmin\" w p pt 7 ps 1.0"
        //                << ", '-' u 1:2 title \"pmax\" w p pt 7 ps 1.0"
        //                << ", '-' u 1:2 title \"test points\" w p pt 2 ps 1.0"
        //                << ", '-' u 1:2 title \"lines\" w l lt 3 lw 0.8"
        //                << "\n";
        //        for (auto& sp : samplesRot) {
        //            gp << sp(0) << " " << sp(1) << "\n";
        //        }
        //        gp << "e\n";
        //        for (auto& sp : samplesTra) {
        //            gp << sp(0) << " " << sp(1) << "\n";
        //        }
        //        gp << "e\n";
        //        gp << pmin.transpose() << "\ne\n" << pmax.transpose() << "\ne\n";
        //        for (auto& tp : testPoints) {
        //            gp << tp.transpose() << "\n";
        //        }
        //        gp << "e\n"; 
        //        plotLine(gp, 0.0, bbXL, (M_PI / 180.0 * 70.0));
        //        plotLine(gp, M_PI, -bbXU, (M_PI / 180.0 * 70.0));
        //        plotLine(gp, 0.5 * M_PI, bbYL, (M_PI / 180.0 * 70.0));
        //        plotLine(gp, 1.5 * M_PI, -bbYU, (M_PI / 180.0 * 70.0));
        //        gp << "e\n";
    }

    return 0;
}

void plotLine(std::ostream& out, double theta, double rho, double l) {
    double ct, st, x1, y1, x2, y2, x0, y0, xt, yt, dir;
    ct = cos(theta);
    st = sin(theta);

    x1 = rho * cos(theta) - l * sin(theta);
    y1 = rho * sin(theta) + l * cos(theta);
    x2 = rho * cos(theta) + l * sin(theta);
    y2 = rho * sin(theta) - l * cos(theta);
    x0 = rho * cos(theta);
    y0 = rho * sin(theta);
    xt = x0 + dir * cos(theta);
    yt = y0 + dir * sin(theta);

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
    //out << mid(0) << " " << mid(1) << "\n" << ver(0) << " " << ver(1) << "\n\n";
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

void branchInterval(const Interval& father, Interval& childL, Interval& childR, int dim) {
    double middle;
    dim = dim % 3;

    //GLORES_PRINT_VARIABLE(dim);

    childL.poseLow = father.poseLow;
    childL.poseUpp = father.poseUpp;
    childR.poseLow = father.poseLow;
    childR.poseUpp = father.poseUpp;

    middle = 0.5 * (father.poseLow(dim) + father.poseUpp(dim));
    childL.poseUpp(dim) = middle;
    childR.poseLow(dim) = middle;
    childL.depth = father.depth + 1;
    childR.depth = father.depth + 1;
    childL.parent = father.id;
    childR.parent = father.id;
}

std::ostream& operator<<(std::ostream& out, const Interval& interval) {
    out << "id " << interval.id << " ";
    for (int i = 0; i < 3; ++i) {
        if (i > 0) out << " x ";
        out << "[" << interval.poseLow(i) << ":" << interval.poseUpp(i) << "]";
    }
    out << " dmin " << interval.dmin << " dmax " << interval.dmax;
    out << " parent " << interval.parent << " depth " << interval.depth;
    return out;
}
