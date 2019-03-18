/**
 * ARS - Angular Radon Spectrum 
 * Copyright (C) 2017 Dario Lodi Rizzini.
 *
 * ARS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * ARS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with ARS.  If not, see <http://www.gnu.org/licenses/>.
 */
#include <glores/PointSetGlobalRegistration.h>
#include <glores/geometry.h>
#include <iomanip>  // std::setw()
#include <fstream>
#include <glores/thirdparty/gnuplot-iostream.h>

namespace glores {

    // --------------------------------------------------------
    // OPTIMIZATION
    // --------------------------------------------------------

    PointSetGlobalRegistration::PointSetGlobalRegistration()
    : nodes_(), dim_(3), ytol_(0.01), xtolOn_(true), ytolOn_(false), inlierRatio_(1.0),
    depthMax_(-1), verboseLevel_(0), iterationCounter_(0), iterationMax_(-1),
    relaxationBoundOn_(true), cheapBoundQueueOn_(true),
    statDminQueueLenAvg_(0.0), statDminQueueLenNum_(0), relaxationThres_(3) {
        xtol_ = Eigen::VectorXd::Constant(3, (double) 0.001);
        relaxationThres_ << 0.05, 0.05, (M_PI / 180.0);
    }

    PointSetGlobalRegistration::PointSetGlobalRegistration(double xtol, double ytol)
    : nodes_(), dim_(3), ytol_(ytol), xtolOn_(true), ytolOn_(false), inlierRatio_(1.0),
    depthMax_(-1), verboseLevel_(0), iterationCounter_(0), iterationMax_(-1),
    relaxationBoundOn_(true), cheapBoundQueueOn_(true),
    statDminQueueLenAvg_(0.0), statDminQueueLenNum_(0), relaxationThres_(3) {
        xtol_ = Eigen::VectorXd::Constant(3, (double) xtol);
        relaxationThres_ << 0.05, 0.05, M_PI / 180.0;
    }

    PointSetGlobalRegistration::PointSetGlobalRegistration(const Eigen::VectorXd& xtol, double ytol)
    : nodes_(), dim_(3), xtol_(xtol), ytol_(ytol), xtolOn_(true), ytolOn_(false), inlierRatio_(1.0),
    depthMax_(-1), verboseLevel_(0), iterationCounter_(0), iterationMax_(-1),
    relaxationBoundOn_(true), cheapBoundQueueOn_(true),
    statDminQueueLenAvg_(0.0), statDminQueueLenNum_(0), relaxationThres_(3) {
        relaxationThres_ << 0.05, 0.05, M_PI / 180.0;
    }

    PointSetGlobalRegistration::~PointSetGlobalRegistration() {
    }

    bool PointSetGlobalRegistration::stopIteration(const IntervalBound& global) const {
        //GLORES_PRINT_VARIABLE(iterationMax_);
        return (global.yupper - global.ylower < ytol_ * global.yupper) || (iterationMax_ >= 0 && iterationCounter_ > iterationMax_);
    }

    void PointSetGlobalRegistration::findGlobalMin(const Eigen::VectorXd& xmin, const Eigen::VectorXd& xmax, Eigen::VectorXd& x, double& ylower, double& yupper) {
        IndexLowerBoundLess comp(&nodes_);
        LeastIndexLowerBoundFirstQueue queue(comp);

        IndexLowerBoundMore compReverse(&nodes_);
        MaxLowerBoundQueue queueMaxLB(compReverse);
        //LeastIndexLowerBoundFirstQueue leaves(comp);
        MinLowerBoundQueue leaves(comp);
        IntervalBound left, right, optimum;
        int icurr, prunedNodeNum, leafNodeNum;
        double cutRatio, volumeTot, volumeCut;
        Gnuplot gp("gnuplot -persist");
        BoundStats bstat;
        bool optimumLowerFromRelaxation;

        if (!checkMinMax(xmin, xmax)) {
            std::cerr << __FILE__ << "," << __LINE__ << ": invalid interval" << std::endl;
            return;
        }

        // At least one tolerance must be enabled to guarantee a stop condition! 
        // If not, default choice (xtol_ on) is chosen. 
        if (!xtolOn_ && !ytolOn_) {
            ytolOn_ = true;
        }

        // Clear the BB tree
        nodes_.clear();

        // Initializes the queue with the minimum and maximum
        optimum.xmin = xmin;
        optimum.xmax = xmax;
        optimum.idx = 0;
        optimum.depth = 0;
        optimum.parent = -1;
        optimum.brother = -1;
        optimum.isLeaf = true;
        GLORES_PRINT_MSG("**** calling computeAllPointDistanceLowers()");
        computeAllPointDistanceLowers(optimum);
        optimum.relaxationBoundUsed = findLU(optimum, optimum.ylower, optimum.yupper);

        //if (verboseLevel_ > 0) {
        std::cout << "initial best bounds: global " << optimum << "; stopIteration(global) " << stopIteration(optimum) << std::endl;
        //}

        nodes_.push_back(optimum);
        queue.push(optimum.idx);
        leaves.push(optimum.idx);
        queueMaxLB.push(optimum.idx);

        // Counter of iterations: variable usefull to monitor the execution and 
        // for selecting output
        iterationCounter_ = 0;
        prunedNodeNum = 0;
        leafNodeNum = 1;
        cutRatio = 0.0;
        volumeTot = intervalVolume(xmin, xmax);

        // Extracts candidate interval from queue and split them to refine lower/upper value estimation
        while (!queue.empty() && !stopIteration(optimum)) {
            icurr = queue.top();
            queue.pop();

            if (iterationCounter_ % 100 == 0 || queue.size() < 10 || verboseLevel_ > 1) {
                std::cout << "iteration " << iterationCounter_ << " bounds " << optimum.ylower << ":" << optimum.yupper
                        << " queue " << queue.size() << " nodes " << nodes_.size() << " leaves " << leafNodeNum << ", stop? " << stopIteration(optimum)
                        << " prunedNodeNum " << prunedNodeNum << " cutRatio " << cutRatio << std::endl;
                if (verboseLevel_ > 1) {
                    GLORES_PRINT_MSG("node " << icurr << ", LB " << optimum.ylower << " UB " << optimum.yupper << " -> queue.size() " << queue.size() << ", nodes " << nodes_.size());
                }
            }

            //            // For debugging on our toy problem:
            //            if (nodes_[icurr].xmin(0) < 2.1 && 2.1 < nodes_[icurr].xmax(0) &&
            //                    nodes_[icurr].xmin(1) < -3.0 && -3.0 < nodes_[icurr].xmax(1) &&
            //                    nodes_[icurr].xmin(2) < 0.73304 && 0.73304 < nodes_[icurr].xmax(2)) {
            //                std::cout << "\n*****\nInterval containing solution:  " << nodes_[icurr] << "\n" << std::endl;
            //            }

            // It processes the interval nodes_[icurr] only if it can improve the current solution. 
            if (nodes_[icurr].ylower <= optimum.yupper * (1.0 - ytol_)) {

                // Computes UPPER BOUND only every 20 iterations!!!
                //                if (iterationCounter_ % 20 == 0) {
                //                    nodes_[icurr].yupper = computeUpperBoundNesterMead2(nodes_[icurr]);
                //                    GLORES_PRINT_MSG("computed upper bound at iteration " << iterationCounter_
                //                            << " for node " << icurr << ": UB " << nodes_[icurr].yupper);
                //                }
                //nodes_[icurr].yupper = computeUpperBound(nodes_[icurr]);
                // If the upper bound value is below the optimum, the results is refined...
                if (nodes_[icurr].yupper <= optimum.yupper) {
                    //nodes_[icurr].yupper = computeUpperBoundNelderMead(nodes_[icurr]);
                }

                // Optimum is updated if the nodes_[icurr] UB (upper bound) is 
                // less then the current upper bound estimate
                if (nodes_[icurr].yupper <= optimum.yupper) {
                    optimum.xmin = nodes_[icurr].xmin;
                    optimum.xmax = nodes_[icurr].xmax;
                    optimum.yupper = nodes_[icurr].yupper;
                    optimum.depth = nodes_[icurr].depth;
                    optimum.idx = icurr;

                    std::cout << "\n**** Update UPPER BOUND " << optimum.yupper << "\n"
                            << "  " << optimum << "\n";
                    //computeUpperBoundNesterMead2(optimum, true);
                    std::cout << std::endl;

                    // Cut the nodes with lower bound larger than upper bound.
                    // Actually the node is nut cut, but its memory consuming dmin/dmax tables 
                    // are disabled.
                    // Goal: to limit memory occupancy
                    while (!queueMaxLB.empty() && nodes_[queueMaxLB.top()].ylower > optimum.yupper) {
                        int deleteId = queueMaxLB.top();
                        queueMaxLB.pop();
                        //                        GLORES_PRINT_MSG("removing node " << deleteId << " with LB " << nodes_[queueMaxLB.top()].ylower 
                        //                                << " > global UB " << optimum.yupper);

                        // Free the memory of the queues pointDistancesMin() and pointDistancesMax
                        // for all the removed nodes
                        nodes_[deleteId].isLeaf = false;
                        nodes_[deleteId].pointDistancesMax.clear();
                        nodes_[deleteId].pointDistancesMin.clear();

                        prunedNodeNum++;
                        volumeCut = intervalVolume(nodes_[deleteId].xmin, nodes_[deleteId].xmax);
                        cutRatio += volumeCut / volumeTot;
                    }
                }

                GLORES_ASSERT(optimum.ylower < optimum.yupper + 1e-10);

                // Improves the estimate of optimum LB (lower) 
                //                if (nodes_[icurr].ylower >= optimum.ylower) {
                //                    optimum.ylower = nodes_[icurr].ylower;
                //                }
                //updateGlobalLowerBound(icurr, optimum.ylower);

                //                // Checks if leaves queue really contains a minimum
                //                for (auto it = leaves.ordered_begin(); it != leaves.ordered_end(); ++it) {
                //                    if (nodes_[*it].ylower < nodes_[leaves.top()].ylower) {
                //                        GLORES_PRINT_MSG("PQUEUE node " << nodes_[*it] << "\n"
                //                                << "  less than leaves.top() " << nodes_[leaves.top()]);
                //                    }
                //                }
                //
                //                // Alternative computation of global LB: find the minimum of LB of leaves
                while (!leaves.empty() && !nodes_[leaves.top()].isLeaf) {
                    //std::cout << "  removing non-leaf " << leaves.top() << " (lb " << nodes_[leaves.top()].ylower << ") from leaf queue\n";
                    leaves.pop();
                }
                // WARNING: lower bound with relaxation ("not-so-cheap") is not necessarily monotonic!!!!
                //                if (nodes_[leaves.top()].ylower < (1.0 - ytol_) * optimum.ylower) {
                //                    GLORES_PRINT_MSG("new LB in node " << leaves.top() << " " << nodes_[leaves.top()].ylower << " < previous LB " << optimum.ylower
                //                            << "\n  " << nodes_[leaves.top()] << "\n  global LB must increase!!!");
                //                    int par = nodes_[leaves.top()].parent;
                //                    std::cout << "  parent of " << leaves.top() << ": " << nodes_[par] << std::endl;
                //
                //                    std::cout << "\nPrinting leaves queue: top " << leaves.top() << "\n";
                //                    for (auto it = leaves.ordered_begin(); it != leaves.ordered_end(); ++it) {
                //                        std::cout << "  " << *it << ": " << nodes_[*it] << "\n";
                //                    }
                //                    GLORES_ASSERT(0);
                //                }
                if (optimum.ylower < nodes_[leaves.top()].ylower) {
                    optimum.ylower = nodes_[leaves.top()].ylower;
                    optimumLowerFromRelaxation = nodes_[leaves.top()].relaxationBoundUsed;
                }
                //std::cout << "global LB " << nodes_[leaves.top()].ylower << " in node " << leaves.top() << std::endl;

                // BRUTE-FORCE computation of LOWER BOUND
                //                int leafLowerBoundIdx;
                //                optimum.ylower = findBestLeafLowerBound(leafLowerBoundIdx);
                //                std::cout << "global LB " << optimum.ylower << " in " << leafLowerBoundIdx
                //                        << " vs leaves.top() " << leaves.top() << " " << nodes_[leaves.top()].ylower << std::endl;

                if (optimum.ylower >= optimum.yupper + 1e-10) {
                    std::cout << "\n***\nPoint Distances for node\n  " << nodes_[leaves.top()] << "\n";
                    for (int is = 0; is < pointSrc_.size(); ++is) {
                        printPointDistance(std::cout, nodes_[leaves.top()], is);
                    }
                    std::cout << "\n***\nEmptying leaves queue:\n";
                    while (!leaves.empty()) {
                        std::cout << "  " << leaves.top() << " LB " << nodes_[leaves.top()].ylower << ": " << nodes_[leaves.top()] << "\n";
                        leaves.pop();
                    }

                    GLORES_PRINT_VARIABLE(optimum.ylower);
                    GLORES_PRINT_VARIABLE(optimum.yupper);
                    GLORES_ASSERT(optimum.ylower < optimum.yupper + 1e-10);
                }

                // Splits curr into intervals left and right, if the interval does not satisfy stop criteria
                //if (isIntervalSplitActive(nodes_[icurr])) {
                if (true) {
                    nodes_[icurr].isLeaf = false;
                    leafNodeNum--;
                    //GLORES_PRINT_MSG("splitting interval " << icurr << " " << curr);
                    splitInterval(nodes_[icurr], (nodes_[icurr].depth % dim_), left, right);
                    // Inserts left and right into nodes_ 
                    left.idx = nodes_.size();
                    left.brother = left.idx + 1;
                    left.isLeaf = true;
                    if (cheapBoundQueueOn_) {
                        updateAllPointDistanceLowers(nodes_[icurr], left, true);
                    } else {
                        computeAllPointDistanceLowersNoQueue(left);
                    }
                    left.relaxationBoundUsed = findLU(left, left.ylower, left.yupper);
                    // Upperbound disabled in findLU(): too inefficient to compute it for every node!!!
                    if (left.ylower <= (1.0 - ytol_) * optimum.yupper) {
                        nodes_.push_back(left);
                        leaves.push(left.idx);
                        queue.push(left.idx);
                        queueMaxLB.push(left.idx);
                    } else {
                        prunedNodeNum++;
                        volumeCut = intervalVolume(left.xmin, left.xmax);
                        cutRatio += volumeCut / volumeTot;
                        if (verboseLevel_ > 0) {
                            GLORES_PRINT_MSG("PRUNING node " << left << " optimum.yupper " << optimum.yupper << ": cut ratio " << cutRatio);
                        }
                    }

                    right.idx = nodes_.size();
                    right.brother = left.idx;
                    right.isLeaf = true;
                    //computeAllPointDistanceLowers(right);
                    //updateAllPointDistanceLowers(nodes_[icurr], right, false);
                    if (cheapBoundQueueOn_) {
                        updateAllPointDistanceLowers(nodes_[icurr], right, false);
                    } else {
                        computeAllPointDistanceLowersNoQueue(right);
                    }
                    right.relaxationBoundUsed = findLU(right, right.ylower, right.yupper);
                    // Upperbound disabled in findLU(): too inefficient to compute it for every node!!!
                    if (right.ylower < (1.0 - ytol_) * optimum.yupper) {
                        nodes_.push_back(right);
                        leaves.push(right.idx);
                        queue.push(right.idx);
                        queueMaxLB.push(right.idx);
                    } else {
                        prunedNodeNum++;
                        volumeCut = intervalVolume(right.xmin, right.xmax);
                        cutRatio += volumeCut / volumeTot;
                        if (verboseLevel_ > 0) {
                            GLORES_PRINT_MSG("PRUNING node " << right << " optimum.yupper " << optimum.yupper << ": cut ratio " << cutRatio);
                        }
                    }

                    //                    if (nodes_[icurr].ylower > nodes_[left.idx].ylower) {
                    //                        GLORES_PRINT_MSG("nodes_[icurr].ylower " << nodes_[icurr].ylower << " <= nodes_[left.idx].ylower " << nodes_[left.idx].ylower);
                    //                        GLORES_ASSERT(nodes_[icurr].ylower <= nodes_[left.idx].ylower);
                    //                    }
                    //                    if (nodes_[icurr].ylower > nodes_[right.idx].ylower) {
                    //                        GLORES_PRINT_MSG("nodes_[icurr].ylower " << nodes_[icurr].ylower << " <=  nodes_[right.idx].ylower " << nodes_[right.idx].ylower);
                    //                        GLORES_ASSERT(nodes_[icurr].ylower <= nodes_[right.idx].ylower);
                    //                    }
                    //                    if (nodes_[leaves.top()].ylower > nodes_[left.idx].ylower) {
                    //                        GLORES_PRINT_MSG("nodes_[leaves.top()].ylower " << nodes_[leaves.top()].ylower << " <= nodes_[left.idx].ylower " << nodes_[left.idx].ylower);
                    //                        GLORES_ASSERT(nodes_[leaves.top()].ylower <= nodes_[left.idx].ylower);
                    //                    }
                    //                    if (nodes_[leaves.top()].ylower > nodes_[right.idx].ylower) {
                    //                        GLORES_PRINT_MSG("nodes_[leaves.top()].ylower " << nodes_[leaves.top()].ylower << " <= nodes_[right.idx].ylower " << nodes_[right.idx].ylower);
                    //                        GLORES_ASSERT(nodes_[leaves.top()].ylower <= nodes_[right.idx].ylower);
                    //                    }

                    if (verboseLevel_ > 1) {
                        std::cout << " left  " << left << std::endl;
                        std::cout << " right " << right << std::endl;
                    }
                    leafNodeNum += 2;
                } else if (verboseLevel_ > 0) {
                    GLORES_PRINT_MSG("interval " << nodes_[icurr] << " NO BRANCH");
                }
            } else {
                // Free the memory of the queues pointDistancesMin() and pointDistancesMax
                // for all the removed nodes
                nodes_[icurr].isLeaf = false;
                nodes_[icurr].pointDistancesMax.clear();
                nodes_[icurr].pointDistancesMin.clear();

                prunedNodeNum++;
                volumeCut = intervalVolume(nodes_[icurr].xmin, nodes_[icurr].xmax);
                cutRatio += volumeCut / volumeTot;
                if (verboseLevel_ > 0) {
                    GLORES_PRINT_MSG("PRUNING node " << nodes_[icurr] << " optimum.yupper " << optimum.yupper << ": cut ratio " << cutRatio);
                }
            }

            bstat.iteration = iterationCounter_;
            bstat.lowerBound = optimum.ylower;
            bstat.upperBound = optimum.yupper;
            bstat.relaxationBoundOn = optimumLowerFromRelaxation; //nodes_[icurr].relaxationBoundUsed;
            bstat.queueNum = queue.size();
            boundsStats_.push_back(bstat);

            iterationCounter_++;
        }
        GLORES_PRINT_VARIABLE(stopIteration(optimum));
        GLORES_PRINT_VARIABLE(optimum);
        x = 0.5 * (optimum.xmin + optimum.xmax);
        ylower = optimum.ylower;
        yupper = optimum.yupper;
    }

    bool PointSetGlobalRegistration::findLU(const IntervalBound& ib, double& ylower, double& yupper) {
        Eigen::VectorXd dims;
        double dimMax, ylowerRelax;
        bool useRelaxationLB;
        bool relaxationBoundUsed;

        // Finds the maximum dimension of cell
        dims = ib.xmax - ib.xmin;
        useRelaxationLB = true;
        for (int i = 0; i < dims.rows(); ++i) {
            if (dims(i) > relaxationThres_(i)) {
                dimMax = dims(i);
                useRelaxationLB = false;
                break;
            }
        }

        // Uses different bounds
        ylower = computeLowerBound(ib);
        //if (dimMax < 0.05) {
        relaxationBoundUsed = false;
        if (relaxationBoundOn_ && useRelaxationLB) {
            //GLORES_PRINT_MSG("NOT-SO-CHEAP");
            ylowerRelax = computeLowerBoundRelaxation(ib);
            //GLORES_PRINT_MSG("cheap lower bound " << ylower << ", no-so-cheap " << ylowerRelax);
            //ylower = std::max(ylower, ylowerRelax);
            if (ylowerRelax > ylower) {
                //GLORES_PRINT_MSG("RELAXATION BOUND USED: ylowerRelax " << ylowerRelax << " vs " << ylower);
                ylower = ylowerRelax;
                relaxationBoundUsed = true;
            }
            //            else {
            //                GLORES_PRINT_MSG("RELAXATION BOUND NOT USED");
            //            }
        }

        yupper = computeUpperBound(ib);
        //yupper = computeUpperBoundNesterMead2(ib);
        // Inherits Upper Bound from parent!!!!
        //        if (ib.parent < 0 || ib.parent >= nodes_.size()) {
        //            yupper = computeUpperBoundNesterMead2(ib);
        //        } else {
        //            yupper = nodes_[ib.parent].yupper;
        //        }
        return relaxationBoundUsed;
    }

    double PointSetGlobalRegistration::evaluateFunction(const Eigen::VectorXd& x) const {
        std::vector<double> upperBests;
        int upperBestNum;
        Transformation2d dstTsrc;
        std::vector<std::pair<int, int> > associations;
        double diff, distMax, yupper;

        poseToIsometry2d(x, dstTsrc);
        distMax = 1e+6;
        associateNN(pointSrc_, pointDst_, dstTsrc, distMax, associations);

        //GLORES_PRINT_MSG("associations in " << __FUNCTION__);
        for (auto& a : associations) {
            diff = (pointDst_[a.second] - dstTsrc * pointSrc_[a.first]).norm();
            upperBests.push_back(diff * diff);
            //std::cout << "  " << a.first << " " << a.second << "\n";
        }

        std::sort(upperBests.begin(), upperBests.end());
        yupper = 0.0;
        upperBestNum = ceil(inlierRatio_ * pointSrc_.size());
        for (int i = 0; i < upperBestNum && i < upperBests.size(); ++i) {
            //std::cout << " yupper += " << upperBests[i] << "\n";
            yupper += upperBests[i];
        }
        return yupper;
    }

    double PointSetGlobalRegistration::evaluateFunctionFast(const Eigen::VectorXd& x, const IntervalBound& ib) const {
        std::vector<double> upperBests;
        int upperBestNum;
        //Transformation2d dstTsrc;
        Eigen::Matrix2d R;
        Eigen::Vector2d t;
        Point2d pTransf;
        //std::vector<std::pair<int, int> > associations;
        double dist, distMin, yupper;

        //poseToIsometry2d(x, dstTsrc);
        double cv, sv;
        cv = cos(x(2));
        sv = sin(x(2));
        R << cv, -sv,
                sv, cv;
        t << x(0), x(1);
        //GLORES_PRINT_MSG("x " << x.transpose() << "\n R\n" << R << "\n t " << t.transpose());
        for (int is = 0; is < pointSrc_.size(); ++is) {
            distMin = 1e+6;
            //pTransf = dstTsrc * pointSrc_[is];
            pTransf = R * pointSrc_[is] + t;
            for (auto& pd : ib.pointDistancesMin[is]) {
                dist = (pointDst_[pd.dstId] - pTransf).norm();
                if (dist < distMin) {
                    distMin = dist;
                }
            }
            upperBests.push_back(distMin * distMin);
        }

        std::sort(upperBests.begin(), upperBests.end());
        yupper = 0.0;
        upperBestNum = ceil(inlierRatio_ * pointSrc_.size());
        for (int i = 0; i < upperBestNum && i < upperBests.size(); ++i) {
            yupper += upperBests[i];
        }
        return yupper;
    }

    // ----------------------------------------------------------
    // PRIVATE
    // ----------------------------------------------------------

    bool PointSetGlobalRegistration::isIntervalSplitActive(const IntervalBound& ib) const {
        //std::cout << "ib.xmin.size() " << ib.xmin.size() << ", ib.xmax.size() " << ib.xmax.size() << ", dim_ " << dim_ << std::endl;
        assert(xtol_.size() == dim_);

        if (verboseLevel_ > 1) {
            std::cout << __FUNCTION__ << ": interval to split " << ib << std::endl;
        }

        // Checking tolerance on dimension of domain
        if (ib.xmin.size() != dim_) {
            std::cerr << __FILE__ << "," << __LINE__ << ": invalid xmin size " << ib.xmin.size() << " must be equal to " << dim_ << std::endl;
            return false;
        }
        if (ib.xmax.size() != dim_) {
            std::cerr << __FILE__ << "," << __LINE__ << ": invalid xmax size " << ib.xmax.size() << " must be equal to " << dim_ << std::endl;
            return false;
        }

        if (iterationMax_ >= 0 && iterationCounter_ < iterationMax_) {
            if (verboseLevel_ > 1) {
                GLORES_PRINT_MSG("reached maximum number of iterations " << iterationCounter_ << " == " << iterationMax_ << ": STOP branching!");
            }
            return false;
        }

        // Checking if current interval depth is greater than maximum depth (a limit to 
        // avoid infinite iterations during debug...)
        // If depthMax_ < 0 then there is no bound on depth
        if (depthMax_ >= 0 && ib.depth > depthMax_) {
            if (verboseLevel_ > 1) {
                GLORES_PRINT_MSG("depth " << ib.depth << " reached " << depthMax_ << ": STOP branching!");
            }
            return false;
        }

        // Checking tolerance on dimension of domain
        if (xtolOn_) {
            bool splitExists = false;
            for (int i = 0; i < dim_ && !splitExists; ++i) {
                if (verboseLevel_ > 3) {
                    GLORES_PRINT_MSG("interval on coordinate " << i << "[" << ib.xmin(i) << "," << ib.xmax(i) << "] greater than tolerance " << xtol_(i));
                }
                if (ib.xmax(i) - ib.xmin(i) > xtol_(i)) {
                    splitExists = true;
                }
            }
            if (!splitExists) {
                return false;
            }
        }
        // Checking tolerance on approximation of function value
        if (ytolOn_ && (ib.yupper - ib.ylower < ytol_ * ib.yupper)) {
            if (verboseLevel_ > 3) {
                GLORES_PRINT_VARIABLE(ib.yupper - ib.ylower);
                GLORES_PRINT_VARIABLE(ytol_ * ib.yupper);
            }
            return false;
        }
        return true;
    }

    bool PointSetGlobalRegistration::checkMinMax(const Eigen::VectorXd& xmin, const Eigen::VectorXd& xmax) const {
        //std::cout << "checking xmin [" << xmin.transpose() << "[, xmax [" << xmax.transpose() << "]" << std::endl;

        // Check dimension
        if (xmin.size() != dim_) {
            std::cerr << __FILE__ << "," << __LINE__ << ": invalid xmin size " << xmin.size() << " must be equal to " << dim_ << std::endl;
            return false;
        }
        if (xmax.size() != dim_) {
            std::cerr << __FILE__ << "," << __LINE__ << ": invalid xmax size " << xmax.size() << " must be equal to " << dim_ << std::endl;
            return false;
        }
        // Check that xmin <= xmax on all the components
        for (int i = 0; i < dim_; ++i) {
            if (xmin(i) > xmax(i)) {
                std::cerr << __FILE__ << "," << __LINE__ << ": invalid xmin(" << i << ") = " << xmin(i)
                        << " > " << xmax(i) << " = xmax(" << i << ")" << std::endl;
                return false;
            }
        }
        return true;
    }

    double PointSetGlobalRegistration::intervalVolume(const Eigen::VectorXd& xmin, const Eigen::VectorXd& xmax) const {
        GLORES_ASSERT(checkMinMax(xmin, xmax));
        Eigen::VectorXd diff = xmax - xmin;
        double vol = 1.0;
        for (int i = 0; i < dim_; ++i) {
            vol = vol * diff(i);
        }
        return vol;
    }

    void PointSetGlobalRegistration::splitInterval(const IntervalBound& intervCur, int d, IntervalBound& intervLow, IntervalBound& intervUpp) const {
        intervLow = intervCur;
        intervUpp = intervCur;
        double mid = 0.5 * (intervCur.xmin(d) + intervCur.xmax(d));
        int depthNew = intervCur.depth;
        intervLow.xmax(d) = mid;
        intervUpp.xmin(d) = mid;
        intervLow.depth = depthNew + 1;
        intervUpp.depth = depthNew + 1;
        intervLow.parent = intervCur.idx;
        intervUpp.parent = intervCur.idx;
        // Indices defined outside!
        intervLow.idx = -1;
        intervUpp.idx = -1;
    }

    // --------------------------------------------------------------
    // BOUNDS METHODS
    // --------------------------------------------------------------

    void PointSetGlobalRegistration::computeAllPointDistanceLowers(IntervalBound& ib) {
        PointDistance pdMin;
        double dmax;

        ib.pointDistancesMin.clear();
        ib.pointDistancesMin.resize(pointSrc_.size());
        ib.pointDistancesMax.clear();
        ib.pointDistancesMax.resize(pointSrc_.size());

        //        pdMin.depth = ib.depth;
        //        pdMin.intervalId = ib.idx;
        for (int is = 0; is < pointSrc_.size(); ++is) {
            //GLORES_PRINT_VARIABLE(is);
            //            pdMin.srcId = is;
            // Inserts a fake upper bounf
            //ib.pointDistancesMax[is].srcId = is;
            ib.pointDistancesMax[is].dstId = -1;
            ib.pointDistancesMax[is].distBound = 1e+6;
            for (int id = 0; id < pointDst_.size(); ++id) {
                pdMin.dstId = id;
                if (!computeDistanceMinMax3(pointSrc_[is], pointDst_[id], ib.xmin, ib.xmax, pdMin.distBound, dmax)) {
                    GLORES_PRINT_ERROR("invalid computeDistanceMinMax3()");
                }
                //GLORES_PRINT_MSG(" computeDistanceMinMax3(" << is << "," << id << ") in [" << pdLower.distBound << "," << distUpper << "]");
                ib.pointDistancesMin[is].push_back(pdMin);

                // Updates upper bound
                if (ib.pointDistancesMax[is].distBound > dmax) {
                    ib.pointDistancesMax[is].dstId = id;
                    ib.pointDistancesMax[is].distBound = dmax;
                    //ib.pointDistancesMax[is].intervalId = ib.idx;
                }
            }
            std::sort(ib.pointDistancesMin[is].begin(), ib.pointDistancesMin[is].end(), PointDistanceLess());

            //            if (verboseLevel_ > 0) {
            //                printPointDistance(std::cout, ib, is);
            //            }
        }
    }
    
    void PointSetGlobalRegistration::computeAllPointDistanceLowersNoQueue(IntervalBound& ib) {
        PointDistance pdMin;
        double dmin, dmax;

        ib.pointDistancesMin.clear();
        ib.pointDistancesMin.resize(pointSrc_.size());
        ib.pointDistancesMax.clear();
        ib.pointDistancesMax.resize(pointSrc_.size());

        //        pdMin.depth = ib.depth;
        //        pdMin.intervalId = ib.idx;
        for (int is = 0; is < pointSrc_.size(); ++is) {
            //GLORES_PRINT_VARIABLE(is);
            //            pdMin.srcId = is;
            // Inserts a fake upper bounf
            //ib.pointDistancesMax[is].srcId = is;
            ib.pointDistancesMax[is].dstId = -1;
            ib.pointDistancesMax[is].distBound = 1e+6;
            
            pdMin.dstId = -1;
            pdMin.distBound = 1e+6;
            ib.pointDistancesMin[is].push_back(pdMin);
            for (int id = 0; id < pointDst_.size(); ++id) {
                if (!computeDistanceMinMax3(pointSrc_[is], pointDst_[id], ib.xmin, ib.xmax, dmin, dmax)) {
                    GLORES_PRINT_ERROR("invalid computeDistanceMinMax3()");
                }
                //GLORES_PRINT_MSG(" computeDistanceMinMax3(" << is << "," << id << ") in [" << pdLower.distBound << "," << distUpper << "]");
                //ib.pointDistancesMin[is].push_back(pdMin);
                
                // Updates lower bound
                if (dmin < ib.pointDistancesMin[is][0].distBound) {
                    ib.pointDistancesMin[is][0].dstId = id;
                    ib.pointDistancesMin[is][0].distBound = dmin;
                }

                // Updates upper bound
                if (dmax < ib.pointDistancesMax[is].distBound) {
                    ib.pointDistancesMax[is].dstId = id;
                    ib.pointDistancesMax[is].distBound = dmax;
                }
            }
        }
    }

    void PointSetGlobalRegistration::updateAllPointDistanceLowers(IntervalBound& ibFather, IntervalBound& ibChild, bool copyFather) {
        double dminNew, dmaxNew;
        double dminBest;
        int j, id, dstIdDminBest;
        static const double TOL = 1e-8;
        //bool orderChanged;

        //        if (iterationCounter_ <= 1) {
        //            std::cout << "\n";
        //            GLORES_PRINT_MSG("updating child interval " << ibChild);
        //        }

        statDminQueueLenAvg_ = 0.0;
        for (int is = 0; is < pointSrc_.size(); ++is) {
            // Copies all the distances between pointSrc_[is] and the destination points
            // from father to child list
            if (copyFather) {
                ibChild.pointDistancesMin[is] = ibFather.pointDistancesMin[is];
            } else {
                ibChild.pointDistancesMin[is].swap(ibFather.pointDistancesMin[is]);
            }
            id = ibChild.pointDistancesMin[is][0].dstId;
            GLORES_ASSERT(!ibChild.pointDistancesMin[is].empty());

            // Computes the new min/max of distances on the new 
            // (reduced) ibChild domain for the first distance item 
            if (!computeDistanceMinMax3(pointSrc_[is], pointDst_[id], ibChild.xmin, ibChild.xmax, dminNew, dmaxNew)) {
                GLORES_PRINT_ERROR("invalid computeDistanceMinMax3()");
            }

            // Checks these invariants:
            //    dminNew < dmaxNew
            //    ibChild.pointDistancesMin[is][0].distBound <= dminNew
            //    dminNew < bChild.pointDistancesMax[is].distBound
            GLORES_ASSERT(dminNew < dmaxNew + 1e-10);
            if (ibChild.pointDistancesMin[is][0].distBound > dminNew + 1e-10) { //aggiunta tolleranza
                GLORES_PRINT_MSG("dmin of child " << ibChild.idx << " " << dminNew
                        << " < " << "dmin of father " << ibFather.idx << " " << ibChild.pointDistancesMin[is][0].distBound
                        << ": IMPOSSIBLE!");

                std::cout << "src " << is << " " << pointSrc_[is].transpose() << ", dst " << id << " " << pointDst_[id].transpose() << std::endl;

                double fatherLower, fatherUpper;
                computeDistanceMinMax3(pointSrc_[is], pointDst_[id], ibFather.xmin, ibFather.xmax, fatherLower, fatherUpper);
                std::cout << "  checking on old interval father: " << ibFather << "\n"
                        << "  fatherLower " << fatherLower << ", fatherUpper " << fatherUpper << std::endl;
                double childLower, childUpper;
                computeDistanceMinMax3(pointSrc_[is], pointDst_[id], ibChild.xmin, ibChild.xmax, childLower, childUpper);
                std::cout << "  checking on interval child: " << ibChild << "\n"
                        << "  childLower " << childLower << ", childUpper " << childUpper << std::endl;

                std::cout << "\nfather geometry test params\n";
                printGeometryIntersectionParams(std::cout, ibFather, is, id);
                std::ofstream fileFather("geom_params_father.txt");
                printGeometryIntersectionParams(fileFather, ibFather, is, id);
                fileFather.close();

                std::cout << "\nchild geometry test params\n";
                printGeometryIntersectionParams(std::cout, ibChild, is, id);
                std::ofstream fileChild("geom_params_child.txt");
                printGeometryIntersectionParams(fileChild, ibChild, is, id);
                fileChild.close();

                GLORES_ASSERT(dminNew >= ibChild.pointDistancesMin[is][0].distBound);
            }
            //            if (dminNew >= ibChild.pointDistancesMax[is].distBound) {
            //                GLORES_PRINT_MSG("dmin of child " << ibChild.idx << " " << dminNew
            //                        << " < " << "dmax of father " << ibFather.idx << " " << ibChild.pointDistancesMax[is].distBound);
            //                GLORES_ASSERT(dminNew < ibChild.pointDistancesMax[is].distBound);
            //            }

            // Sets the refined value of distance min on first item: 
            //ibChild.pointDistancesMin[is][0].srcId = is;
            ibChild.pointDistancesMin[is][0].dstId = id;
            ibChild.pointDistancesMin[is][0].distBound = dminNew;
            //ibChild.pointDistancesMin[is][0].depth = ibChild.depth;
            //ibChild.pointDistancesMin[is][0].intervalId = ibChild.idx;
            dminBest = dminNew;
            dstIdDminBest = id;

            // Updates the upper bound if the new distance reducesibChild.pointDistancesMin[is][0].srcId = is;
            //ibChild.pointDistancesMax[is].srcId = is;
            ibChild.pointDistancesMax[is].dstId = ibFather.pointDistancesMax[is].dstId;
            ibChild.pointDistancesMax[is].distBound = ibFather.pointDistancesMax[is].distBound;
            //ibChild.pointDistancesMax[is].depth = ibFather.pointDistancesMax[is].depth;
            //ibChild.pointDistancesMax[is].intervalId = ibFather.pointDistancesMax[is].intervalId;
            if (dmaxNew < ibChild.pointDistancesMax[is].distBound) {
                //ibChild.pointDistancesMax[is].srcId = is;
                ibChild.pointDistancesMax[is].dstId = id;
                ibChild.pointDistancesMax[is].distBound = dmaxNew;
                //ibChild.pointDistancesMax[is].depth = ibChild.depth;
                //ibChild.pointDistancesMax[is].intervalId = ibChild.idx;
            }
            //orderChanged = false;

            double radius = pointSrc_[is].norm();
            double beta = atan2(pointSrc_[is].y(), pointSrc_[is].x());
            double arcBeg = beta + ibChild.xmin(2);
            double arcEnd = beta + ibChild.xmax(2);
            Point2d pArcBeg;
            pArcBeg << radius * cos(arcBeg), radius * sin(arcBeg);
            Point2d pArcEnd;
            pArcEnd << radius * cos(arcEnd), radius * sin(arcEnd);
            // Recomputes only the new distances that needs to be updated            
            for (j = 1; j < ibChild.pointDistancesMin[is].size() && ibChild.pointDistancesMin[is][j].distBound < dminBest; ++j) {
                // Recomputes the lower and upper bounds of j-th: the new ylower
                // is greater or equal to the ylower computed on ibFather
                id = ibChild.pointDistancesMin[is][j].dstId;


                if (!computeDistanceMinMax3(pArcBeg, pArcEnd, beta, radius, pointDst_[id], ibChild.xmin, ibChild.xmax, dminNew, dmaxNew)) {
                    GLORES_PRINT_ERROR("invalid computeDistanceMinMax3()");
                }

                // Check invariants: 
                GLORES_ASSERT(dminNew <= dmaxNew);
                GLORES_ASSERT(dminNew >= ibChild.pointDistancesMin[is][j].distBound);

                // Recomputed value of lower bound on interval ibChild
                //ibChild.pointDistancesMin[is][j].srcId = is;
                ibChild.pointDistancesMin[is][j].dstId = id;
                ibChild.pointDistancesMin[is][j].distBound = dminNew;
                //ibChild.pointDistancesMin[is][j].depth = ibChild.depth;
                //ibChild.pointDistancesMin[is][j].intervalId = ibChild.idx;
                // Updates the minimum of lower distances on ibChild
                if (dminNew < dminBest) {
                    dminBest = dminNew;
                    dstIdDminBest = id;
                }

                if (dmaxNew < ibChild.pointDistancesMax[is].distBound) {
                    //ibChild.pointDistancesMax[is].srcId = is;
                    ibChild.pointDistancesMax[is].dstId = id;
                    ibChild.pointDistancesMax[is].distBound = dmaxNew;
                    //ibChild.pointDistancesMax[is].depth = ibChild.depth;
                    //ibChild.pointDistancesMax[is].intervalId = ibChild.idx;
                }
                statDminQueueLenAvg_ += ibChild.pointDistancesMin[is].size();
            }

            // Checks the consistency of dminBest with some tolerance on floating 
            // point numerical error with very small values of distances
            if (dminBest > ibChild.pointDistancesMax[is].distBound && fabs(dminBest - ibChild.pointDistancesMax[is].distBound) > TOL) {
                GLORES_PRINT_ERROR("iteration " << iterationCounter_ << ": "
                        "it MUST be: dminBest " << dminBest
                        << " < ibChild.pointDistancesMax[is].distBound " << ibChild.pointDistancesMax[is].distBound);

                std::cout << "\ndmax geometry test params\n";
                printGeometryIntersectionParams(std::cout, ibChild, is, ibChild.pointDistancesMax[is].dstId);
                std::ofstream fileDMax("geom_params_dmax.txt");
                printGeometryIntersectionParams(fileDMax, ibChild, is, ibChild.pointDistancesMax[is].dstId);
                fileDMax.close();

                std::cout << "\ndmin geometry test params\n";
                printGeometryIntersectionParams(std::cout, ibChild, is, dstIdDminBest);
                std::ofstream fileDMin("geom_params_dmin.txt");
                printGeometryIntersectionParams(fileDMin, ibChild, is, dstIdDminBest);
                fileDMin.close();

                GLORES_ASSERT(dminBest < ibChild.pointDistancesMax[is].distBound);
            }

            std::sort(ibChild.pointDistancesMin[is].begin(), ibChild.pointDistancesMin[is].end(), PointDistanceLess());

            if (ibChild.pointDistancesMin[is][0].distBound > ibChild.pointDistancesMax[is].distBound && fabs(ibChild.pointDistancesMin[is][0].distBound - ibChild.pointDistancesMax[is].distBound) > TOL) {
                std::cout << "\n---\n";
                GLORES_PRINT_MSG("dmin of interval " << ibChild.idx << " " << ibChild.pointDistancesMin[is][0].distBound
                        << " < "
                        << "dmax " << ibChild.pointDistancesMax[is].distBound);
                //<< "dmax set in interval " << ibChild.pointDistancesMax[is].intervalId << " " << ibChild.pointDistancesMax[is].distBound);

                printPointDistance(std::cout, ibChild, is);
                std::cout << "  *** saving files point_distance_min.txt and point_distance_max.txt for debug" << std::endl;

                int dminDst = ibChild.pointDistancesMin[is][0].dstId;
                std::ofstream fileDMin("point_distance_min.txt");
                printGeometryIntersectionParams(fileDMin, ibChild, is, dminDst);
                fileDMin.close();

                int dmaxDstId = ibChild.pointDistancesMax[is].dstId;
                //int dmaxIntervalId = ibChild.pointDistancesMax[is].intervalId;
                //std::ofstream fileDMax("point_distance_max.txt");
                //printGeometryIntersectionParams(fileDMax, nodes_[dmaxIntervalId], is, dmaxDstId);
                //fileDMax.close();

                GLORES_ASSERT(ibChild.pointDistancesMin[is][0].distBound <= ibChild.pointDistancesMax[is].distBound);
            }

            // Finds if there are lower bound items in  ibChild.pointDistanceLowers[is] 
            // with lower bound greater than ibChild.pointDistanceUppers[is].distBound. 
            // In that cased it removes them from the lower bound queue. 
            ibChild.pointDistancesMax[is].distBound += TOL;
            auto removeIt = std::upper_bound(ibChild.pointDistancesMin[is].begin(), ibChild.pointDistancesMin[is].end(),
                    ibChild.pointDistancesMax[is], PointDistanceLess());
            ibChild.pointDistancesMax[is].distBound -= TOL;

            //            for (auto vit = removeIt; vit != ibChild.pointDistancesMin[is].end(); ++vit) {
            //                if (vit->distBound < ibChild.pointDistancesMax[is].distBound) {
            //                    GLORES_PRINT_ERROR("lower distance(" << vit->srcId << "," << vit->dstId << ") " << vit->distBound
            //                            << " < upper " << ibChild.pointDistancesMax[is].distBound << "!!!");
            //                }
            //            }

            if (removeIt == ibChild.pointDistancesMin[is].begin()) {
                GLORES_PRINT_ERROR("removing all the points in list!!");
                for (auto vit = removeIt; vit != ibChild.pointDistancesMin[is].end(); ++vit) {
                    //                    std::cout << "  dmin(" << vit->srcId << "," << vit->dstId << ") " << vit->distBound
                    //                            << " > dmax " << ibChild.pointDistancesMax[is].distBound << "\n";
                    std::cout << "  dmin(" << vit->dstId << ") " << vit->distBound << "\n";
                }
                GLORES_ASSERT(0);
            }


            int countRemained = removeIt - ibChild.pointDistancesMin[is].begin();
            int countRemoved = ibChild.pointDistancesMin[is].size() - countRemained;
            ibChild.pointDistancesMin[is].erase(removeIt, ibChild.pointDistancesMin[is].end());

            if (verboseLevel_ > 1 && iterationCounter_ <= 1) {

                //                GLORES_PRINT_MSG("remained lower bounds " << countRemained << " removed " << countRemoved);
                GLORES_PRINT_MSG("point src " << is << " bounds:");
                std::cout << "  lower[" << ibChild.pointDistancesMin[is].size() << "]: ";
                for (auto& l : ibChild.pointDistancesMin[is]) {
                    //std::cout << l.distBound << " (" << l.srcId << "," << l.dstId << "; depth " << l.depth << "), ";
                    std::cout << l.distBound << " in " << l.dstId << ", ";
                }
                std::cout << std::endl;
                std::cout << "  upper: " << ibChild.pointDistancesMax[is].dstId << " dist " << ibChild.pointDistancesMax[is].distBound << std::endl;
                //                std::cout << "  upper: " << ibChild.pointDistancesMax[is].distBound
                //                        << " (" << ibChild.pointDistancesMax[is].srcId << "," << ibChild.pointDistancesMax[is].dstId
                //                        << "; depth " << ibChild.pointDistancesMax[is].depth << ")" << std::endl;
            }
        }
        //        if (iterationCounter_ <= 1) {
        //            GLORES_PRINT_MSG("LOWER " << computeLowerBound(ibChild) << ", UPPER " << computeUpperBoundIcp(ibChild));
        //        }
    }

    double PointSetGlobalRegistration::computeLowerBound(const IntervalBound& ib) {
        std::vector<double> lowerBests;
        double ylower;
        int lowerBestNum;

        // Computes the lower bound using the best first inlierRatio_
        for (int is = 0; is < ib.pointDistancesMin.size(); ++is) {
            if (!ib.pointDistancesMin[is].empty()) {
                lowerBests.push_back(ib.pointDistancesMin[is][0].distBound * ib.pointDistancesMin[is][0].distBound);
                //lowerBests.push_back(ib.pointDistancesMin[is][0].distBound);
            } else {
                GLORES_PRINT_ERROR("no resisuals for src point " << is);
                GLORES_ASSERT(0);
            }
        }

        std::sort(lowerBests.begin(), lowerBests.end());
        ylower = 0.0;
        lowerBestNum = ceil(inlierRatio_ * ib.pointDistancesMin.size());
        for (int i = 0; i < lowerBestNum; ++i) {
            ylower += lowerBests[i];
        }
        return ylower;
    }

    double PointSetGlobalRegistration::computeUpperBound(const IntervalBound& ib, bool outputOn) {
        //        std::vector<double> upperBests;
        //        int upperBestNum;
        //        Transformation2d dstTsrc;
        //        std::vector<std::pair<int, int> > associations;
        //        double diff, distMax, yupper;

        //        poseToIsometry2d(0.5 * (ib.xmin + ib.xmax), dstTsrc);
        //        //GLORES_PRINT_VARIABLE(0.5 * (ib.xmin + ib.xmax).transpose());
        //        distMax = 1e+6;
        //        //        for (int i = 0; i < 10; ++i) {
        //        associateNN(pointSrc_, pointDst_, dstTsrc, distMax, associations);
        //        //            computeTransform(pointSrc_, pointDst_, associations, dstTsrc);
        //        // Function value:
        //        //            yupper = 0.0;
        //        //            for (auto& a : associations) {
        //        //                diff = (pointDst_[a.second] - dstTsrc * pointSrc_[a.first]).norm();
        //        //                yupper = diff * diff;
        //        //            }
        //        //            std::cout << "  ICP " << i << ": function value " << yupper << std::endl;
        //        //        }
        //
        //        //GLORES_PRINT_MSG("associations: ");
        //        for (auto& a : associations) {
        //            diff = (pointDst_[a.second] - dstTsrc * pointSrc_[a.first]).norm();
        //            upperBests.push_back(diff * diff);
        //            outputOn && std::cout << "  " << a.first << "  " << a.second << ": " << upperBests.back() << std::endl;
        //        }
        //
        //        std::sort(upperBests.begin(), upperBests.end());
        //        yupper = 0.0;
        //        upperBestNum = round(inlierRatio_ * pointSrc_.size());
        //        for (int i = 0; i < upperBestNum; ++i) {
        //            yupper += upperBests[i];
        //        }
        //        outputOn && std::cout << "  upper bound of " << upperBestNum << ": " << yupper << std::endl;

        return evaluateFunctionFast(0.5 * (ib.xmin + ib.xmax), ib);
    }

    //    double PointSetGlobalRegistration::computeUpperBoundNesterMead(const IntervalBound& ib) {
    //        nm::NelderMeadOptimizer optimizer(3, 0.01);
    //        nm::Vector valNM;
    //        std::map<double, nm::Vector> vertices;
    //        int count = 0;
    //
    //        // Inserts the 8 vertices of box into the vertices list (ordered by the function value in its)
    //        GLORES_PRINT_VARIABLE(ib);
    //        valNM.prepare(3);
    //        valNM[0] = ib.xmin(0);
    //        valNM[1] = ib.xmin(1);
    //        valNM[2] = ib.xmin(2);
    //        GLORES_PRINT_VARIABLE(valNM.dimension());
    //        GLORES_ASSERT(valNM.dimension() == 3);
    //        std::cout << "valNM " << valNM[0] << " " << valNM[1] << " " << valNM[2] << std::endl;
    //        vertices.insert(std::make_pair(evaluateFunctionNM(valNM), valNM));
    //        valNM[0] = ib.xmin(0);
    //        valNM[1] = ib.xmin(1);
    //        valNM[2] = ib.xmax(2);
    //        vertices.insert(std::make_pair(evaluateFunctionNM(valNM), valNM));
    //        valNM[0] = ib.xmin(0);
    //        valNM[1] = ib.xmax(1);
    //        valNM[2] = ib.xmax(2);
    //        vertices.insert(std::make_pair(evaluateFunctionNM(valNM), valNM));
    //        valNM[0] = ib.xmin(0);
    //        valNM[1] = ib.xmax(1);
    //        valNM[2] = ib.xmin(2);
    //        vertices.insert(std::make_pair(evaluateFunctionNM(valNM), valNM));
    //        valNM[0] = ib.xmax(0);
    //        valNM[1] = ib.xmin(1);
    //        valNM[2] = ib.xmin(2);
    //        vertices.insert(std::make_pair(evaluateFunctionNM(valNM), valNM));
    //        valNM[0] = ib.xmax(0);
    //        valNM[1] = ib.xmin(1);
    //        valNM[2] = ib.xmax(2);
    //        vertices.insert(std::make_pair(evaluateFunctionNM(valNM), valNM));
    //        valNM[0] = ib.xmax(0);
    //        valNM[1] = ib.xmax(1);
    //        valNM[2] = ib.xmax(2);
    //        vertices.insert(std::make_pair(evaluateFunctionNM(valNM), valNM));
    //        valNM[0] = ib.xmax(0);
    //        valNM[1] = ib.xmax(1);
    //        valNM[2] = ib.xmin(2);
    //        vertices.insert(std::make_pair(evaluateFunctionNM(valNM), valNM));
    //        count = 0;
    //        for (auto it = vertices.begin(); it != vertices.end() && count < 4; ++it, ++count) {
    //            GLORES_ASSERT(it->second.dimension() == 3);
    //            if (count == 0) {
    //                valNM = it->second;
    //                GLORES_PRINT_VARIABLE(valNM.dimension());
    //                GLORES_ASSERT(valNM.dimension() == 3);
    //                std::cout << "valNM " << valNM[0] << " " << valNM[1] << " " << valNM[2] << std::endl;
    //            } else {
    //                std::cout << "  vertex " << count << " [" << it->second[0] << "," << it->second[1] << "," << it->second[2] << "]\n";
    //                optimizer.insert(it->second);
    //            }
    //        }
    //
    //        GLORES_PRINT_MSG("starting " << __FUNCTION__);
    //        while (!optimizer.done()) {
    //            GLORES_ASSERT(valNM.dimension() == 3);
    //            double score = evaluateFunctionNM(valNM);
    //            valNM = optimizer.step(valNM, (float) score);
    //            GLORES_ASSERT(valNM.dimension() == 3);
    //        }
    //        GLORES_ASSERT(valNM.dimension() == 3);
    //        return evaluateFunctionNM(valNM);
    //    }

    double PointSetGlobalRegistration::computeUpperBoundNelderMead(const IntervalBound& ib) {
        Eigen::VectorXd mid;
        SimplexObjectiveFunction fun(this, &ib);

        auto optimizer = Optimization::Local::build_simplex(// Builder to generate the correct optimizer
                fun, // Used to infer function type
                Optimization::Local::make_and_criteria(Optimization::Local::IterationCriterion(10),
                Optimization::Local::RelativeValueCriterion<double>(0.001)));

        mid = 0.5 * (ib.xmin + ib.xmax);
        optimizer.set_start_point(mid); // Starting parameters
        optimizer.set_delta(1); // Simplex size
        optimizer.optimize(fun); // Optimization start

        //GLORES_PRINT_MSG("best Upper Bound in " << optimizer.get_best_parameters().transpose() << ": UB " << optimizer.get_best_value());
        return optimizer.get_best_value();
    }

    double PointSetGlobalRegistration::computeLowerBoundRelaxation(const IntervalBound& ib) {
        Eigen::Vector2d v;
        VectorPoint2d polygonBoundArc, boxVertices;
        Transformation2d dstTsrc;
        double am, ap, bm, bp, cm, cp, sm, sp, th1, th2, tx0, ty0, ct0, st0, tx, ty, theta, ct, st;
        std::vector<std::pair<int, int> > associations;
        double distMax, distLinAll, distLinAllMin;

        distMax = 1e+6;

        // Computes the bound 
        am = ib.xmin(0);
        ap = ib.xmax(0);
        bm = ib.xmin(1);
        bp = ib.xmax(1);
        th1 = ib.xmin(2);
        th2 = ib.xmax(2);
        //translMid << 0.5 * (am + ap), 0.5 * (bm + bp);
        tx0 = 0.5 * (am + ap);
        ty0 = 0.5 * (bm + bp);

        // Finds the lower and upper of cos and sin (respectively cm, cp, sm, sp)
        // and the middle of cosine/sine (s0, c0)
        if (sin(th1) * sin(th2) < 0) {
            cm = -1.0;
        } else {
            cm = std::min(cos(th1), cos(th2));
        }
        cp = std::max(cos(th1), cos(th2));
        if (th1 < 0.5 * M_PI && th2 > 0.5 * M_PI) {
            sp = 1.0;
            sm = std::min(sin(th1), sin(th2));
        } else if (th1 < 1.5 * M_PI && th2 > 1.5 * M_PI) {
            sp = std::max(sin(th1), sin(th2));
            sm = -1.0;
        } else {
            sp = std::max(sin(th1), sin(th2));
            sm = std::min(sin(th1), sin(th2));
        }
        st0 = 0.5 * (sp + sm);
        ct0 = 0.5 * (cp + cm);

        v << am, bm;
        boxVertices.push_back(v);
        v << ap, bm;
        boxVertices.push_back(v);
        v << ap, bp;
        boxVertices.push_back(v);
        v << am, bp;
        boxVertices.push_back(v);

        // Approximating arc on unitary circle (for [cos(a), sin(a)))
        projectArcBound(th1, th2, polygonBoundArc, 2);

        //GLORES_PRINT_MSG(__FUNCTION__ << "\n  arc points " << polygonBoundArc.size() << " rectangle " << boxVertices.size());

        // Iterating on all the possible bounding [arc, rectangle] vertices
        // computes the local linear approximation
        distLinAllMin = 1e+6;
        for (int ia = 0; ia < polygonBoundArc.size(); ++ia) {
            for (int ir = 0; ir < boxVertices.size(); ++ir) {
                tx = boxVertices[ir](0);
                ty = boxVertices[ir](1);
                ct = polygonBoundArc[ia](0);
                st = polygonBoundArc[ia](1);
                //theta = atan2(st, ct);
                //poseToIsometry2d(tx, ty, theta, dstTsrc);
                //associateNN(pointSrc_, pointDst_, dstTsrc, distMax, associations);
                distLinAll = evaluateAllPointSetLinearizedDistances(tx0, ty0, ct0, st0, tx, ty, ct, st, ib);
                if (distLinAll < distLinAllMin) {
                    distLinAllMin = distLinAll;
                }
            }
        }
        return distLinAllMin;
    }

    void PointSetGlobalRegistration::updateGlobalLowerBound(int idx, double& lowerBoundGlobal) {
        int cur, par, bro;
        double lowerBoundBrotherMin;

        GLORES_ASSERT(0 <= idx && idx < nodes_.size());
        cur = idx;
        if (verboseLevel_ > 0) {
            GLORES_PRINT_MSG("updating global LB " << lowerBoundGlobal << " from node " << idx << " with ylower " << nodes_[idx].ylower);
        }
        while (0 <= cur && cur < nodes_.size()) {
            par = nodes_[cur].parent;
            bro = nodes_[cur].brother;

            // Computes the minimum lower bound between node cur and its node brother bro
            if (0 <= bro && bro < nodes_.size()) {
                lowerBoundBrotherMin = std::min(nodes_[cur].ylower, nodes_[bro].ylower);
                if (verboseLevel_ > 0) {
                    std::cout << "  best between " << cur << " (lb " << nodes_[cur].ylower << ") "
                            << "and its brother " << bro << " (lb " << nodes_[bro].ylower << "): " << lowerBoundBrotherMin
                            << " (depth " << nodes_[cur].depth << ")" << std::endl;
                }
            } else {
                lowerBoundBrotherMin = nodes_[cur].ylower;
                if (verboseLevel_ > 0) {
                    std::cout << "  node " << cur << "(lb " << nodes_[cur].ylower << ") no brother: best LB " << lowerBoundBrotherMin << std::endl;
                }
            }

            // Cases:
            // 1) node cur has no parent node (par == -1);
            // 2) the parent node LB is less than lowerBoundBrotherMin (minimum LB of cur and bro): 
            //    parent LB is IMPROVED;
            // 3) the parent node LB is better than cur and bro estimate of LB
            if (par == -1) {
                lowerBoundGlobal = lowerBoundBrotherMin;
                if (verboseLevel_ > 0) {
                    std::cout << "  root reached: GLOBAL LB " << lowerBoundGlobal << std::endl;
                }
                cur = -1;
            } else if (nodes_[par].ylower < lowerBoundBrotherMin) {
                if (verboseLevel_ > 0) {
                    std::cout << "  parent node " << par << " (lb " << nodes_[par].ylower << ") < LB " << lowerBoundBrotherMin << std::endl;
                }
                nodes_[par].ylower = lowerBoundBrotherMin;
                cur = par;
            } else {
                if (verboseLevel_ > 0) {
                    std::cout << "  cannot improve parent node " << par << " (lb " << nodes_[par].ylower << ")\n";
                }
                cur = -1;
            }
        }
    }

    double PointSetGlobalRegistration::findBestLeafLowerBound(int& leafLower) const {
        // Brute force search:
        double lowerBoundMin = 1e+6;
        for (int i = 0; i < nodes_.size(); ++i) {
            if (nodes_[i].isLeaf && nodes_[i].ylower < lowerBoundMin) {
                GLORES_ASSERT(nodes_[i].idx == i);
                leafLower = nodes_[i].idx;
                lowerBoundMin = nodes_[i].ylower;
            }
        }
        return lowerBoundMin;
    }

    void PointSetGlobalRegistration::printPointDistance(std::ostream& out, const IntervalBound& ib, int is) const {
        out << "point src " << is << " bounds (interval " << ib.idx << " depth " << ib.depth << " iteration " << iterationCounter_ << "):\n";
        out << "  dmin[" << ib.pointDistancesMin[is].size() << "]: ";
        for (auto& l : ib.pointDistancesMin[is]) {
            //std::cout << l.distBound << " (" << l.srcId << "," << l.dstId << "; depth " << l.depth << "), ";
            std::cout << l.distBound << " (" << l.dstId << "), ";
        }
        std::cout << std::endl;
        //        std::cout << "  dmax: " << ib.pointDistancesMax[is].distBound
        //                << " (" << ib.pointDistancesMax[is].srcId << "," << ib.pointDistancesMax[is].dstId
        //                << "; depth " << ib.pointDistancesMax[is].depth << ")" << std::endl;
        std::cout << "  dmax: " << ib.pointDistancesMax[is].distBound
                << " (" << ib.pointDistancesMax[is].dstId << ")" << std::endl;
    }

    void PointSetGlobalRegistration::printGeometryIntersectionParams(std::ostream& out, const IntervalBound& ib, int is, int id) {
        out << "px " << pointSrc_[is](0) << "\n";
        out << "py " << pointSrc_[is](1) << "\n";
        out << "qx " << pointDst_[id](0) << "\n";
        out << "qy " << pointDst_[id](1) << "\n";
        out << "txL " << ib.xmin(0) << "\n";
        out << "txU " << ib.xmax(0) << "\n";
        out << "tyL " << ib.xmin(1) << "\n";
        out << "tyU " << ib.xmax(1) << "\n";
        out << "rotL " << (180.0 / M_PI * ib.xmin(2)) << "\n";
        out << "rotU " << (180.0 / M_PI * ib.xmax(2)) << "\n";
    }

    void PointSetGlobalRegistration::printPointSets(std::ostream& out, const Eigen::VectorXd & pose) const {
        Transformation2d dstTsrc;

        poseToIsometry2d(pose, dstTsrc);

        //        out << "set term wxt 0\n";
        //        out << "set size ratio -1\n";
        //        out << "set xlabel \"x\"\n";
        //        out << "set ylabel \"y\"\n";
        //        out << "set key outside\n";
        //        out << "plot "
        //                << "'-' u 1:2 title \"src\" w p pt 7 ps 1.5, "
        //                << "'-' u 1:2 title \"dst\" w p pt 4 ps 1.5, "
        //                << "'-' u 1:2 title \"rot\" w p pt 1 ps 1.1  "
        //                << "\n";
        //        for (auto& p : pointSrc_) {
        //            out << p.transpose() << "\n";
        //        }
        //        out << "e\n";
        //        for (auto& p : pointDst_) {
        //            out << p.transpose() << "\n";
        //        }
        //        out << "e\n";
        //        for (auto& p : pointSrc_) {
        //            out << (dstTsrc * p).transpose() << "\n";
        //        }
        //        out << "e\n";
    }

    double PointSetGlobalRegistration::evaluateFunctionNM(nm::Vector & x) const {
        Eigen::VectorXd vx(3);
        GLORES_ASSERT(x.dimension() == 3);
        vx << x.at(0), x.at(1), x.at(2);
        return evaluateFunction(vx);
    }

    double PointSetGlobalRegistration::evaluatePairwiseDistanceSquare(double tx, double ty, double ct, double st, const Point2d& p, const Point2d& q) {
        //Transformation2d dstTsrc;
        //Eigen::Matrix2d R;
        //Point2d delta, transl;

        //        double theta;
        //        theta = atan2(st, ct);
        //        poseToIsometry2d(tx, ty, theta, dstTsrc);
        //        delta = (q - dstTsrc * p);

        //        R << ct, -st,
        //                st, ct;
        //        transl << tx, ty;
        //        delta = R * p + transl - q;
        //        return delta.dot(delta);
        double deltaX, deltaY;

        deltaX = p(0) * ct - p(1) * st + tx - q(0);
        deltaY = p(0) * st + p(1) * ct + ty - q(1);
        return (deltaX * deltaX + deltaY * deltaY);
    }

    void PointSetGlobalRegistration::evaluatePairwiseDistanceSquareGradient(double tx, double ty, double ct, double st, const Point2d& p, const Point2d & q, Eigen::Vector4d &grad) {
        //Transformation2d dstTsrc;
        //        Eigen::Matrix2d R;
        //        Point2d delta, transl;
        //        Eigen::MatrixXd deltaJac(2, 4);
        //        Eigen::VectorXd grad;
        //        double theta;
        //
        //        //        theta = atan2(st, ct);
        //        //        poseToIsometry2d(tx, ty, theta, dstTsrc);
        //        //        delta = (q - dstTsrc * p);
        //        R << ct, -st,
        //                st, ct;
        //        transl << tx, ty;
        //        delta = R * p + transl - q;
        //
        //        // Jacobian of delta w.r.t. pose (angle represented by its cos and sin)
        //        // \partial delta / \partial (pose(0), pose(1), cos(pose), sin(pose)                
        //        deltaJac << 1, 0.0, p(0), -p(1),
        //                0.0, 1.0, p(1), p(0);
        //
        //        grad = 2.0 * delta.transpose() * deltaJac;
        //        return grad;

        typedef Eigen::Matrix< double, 1, 4 > VectorRow4d;
        //Eigen::MatrixXd grad(1,4);
        //Eigen::VectorXd grad(4);
        //VectorRow4d grad;

        double deltaX, deltaY;

        deltaX = p(0) * ct - p(1) * st + tx - q(0);
        deltaY = p(0) * st + p(1) * ct + ty - q(1);

        grad << 2.0 * deltaX,
                2.0 * deltaY,
                2.0 * (deltaX * p(0) + deltaY * p(1)),
                2.0 * (-deltaX * p(1) + deltaY * p(0));
        //return grad;

        //return grad.transpose();
    }

    double PointSetGlobalRegistration::evaluateAllPointSetLinearizedDistances(double tx0, double ty0, double ct0, double st0, double tx, double ty, double ct, double st, const IntervalBound& ib) {
        //Eigen::VectorXd df(4), delta(4);
        Eigen::Vector4d df, delta;

        double fun, distLin, distTot;
        std::vector<double> distLinMin(pointSrc_.size(), 1e+6);
        int num;

        delta << tx - tx0, ty - ty0, ct - ct0, st - st0;
        //GLORES_PRINT_VARIABLE();

        //        for (int is = 0; is < pointSrc_.size(); ++is) {
        //            for (int id = 0; id < pointDst_.size(); ++id) {
        //                fun = evaluatePairwiseDistanceSquare(tx0, ty0, ct0, st0, pointSrc_[is], pointDst_[id]);
        //                df = evaluatePairwiseDistanceSquareGradient(tx0, ty0, ct0, st0, pointSrc_[is], pointDst_[id]);
        //                // Computes the value of linearized distance
        //                distLin = fun + df.dot(delta);
        //                if (distLin < distLinMin[is]) {
        //                    distLinMin[is] = distLin;
        //                }
        //            }
        //        }

        for (int is = 0; is < pointSrc_.size(); ++is) {
            for (auto& pd : ib.pointDistancesMin[is]) {
                fun = evaluatePairwiseDistanceSquare(tx0, ty0, ct0, st0, pointSrc_[is], pointDst_[pd.dstId]);
                evaluatePairwiseDistanceSquareGradient(tx0, ty0, ct0, st0, pointSrc_[is], pointDst_[pd.dstId], df);
                distLin = fun + df.dot(delta);
                if (distLin < distLinMin[is]) {
                    distLinMin[is] = distLin;
                }
            }
        }

        std::sort(distLinMin.begin(), distLinMin.end());
        num = ceil(inlierRatio_ * pointSrc_.size());
        distTot = 0.0;
        for (int i = 0; i < num; ++i) {
            distTot += distLinMin[i];
        }
        return distTot;
    }

} // end of namespace

std::ostream& operator<<(std::ostream& out, const glores::PointSetGlobalRegistration::IntervalBound& interv) {
    //    out << "idx " << interv.idx << ", depth " << interv.depth << ": domain [" << interv.xmin.transpose() << "; " << interv.xmax.transpose() << "]: "
    //            << "bounds [" << interv.ylower << ", " << interv.yupper << "]";
    out << "idx " << interv.idx << " depth " << interv.depth << " par " << interv.parent << " leaf " << interv.isLeaf << ": domain ";
    for (int i = 0; i < interv.xmin.rows() && i < interv.xmax.rows(); ++i) {
        if (i > 0) out << " x ";
        out << "[" << interv.xmin(i) << ":" << interv.xmax(i) << "]";
    }
    out << ": bounds [" << interv.ylower << ", " << interv.yupper << "]";
    return out;
}

