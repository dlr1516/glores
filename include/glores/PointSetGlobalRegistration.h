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
#ifndef POINT_SET_GLOBAL_REGISTRATION_H
#define POINT_SET_GLOBAL_REGISTRATION_H

#include <iostream>
#include <Eigen/Dense>
#include <memory>
#include <vector>
#include <queue>
#include <set>
#include <functional>

#include <boost/heap/priority_queue.hpp>
#include <boost/heap/binomial_heap.hpp>

//#include <roboptim/core/function.hh>


#include <glores/fileutils.h>
#include <glores/thirdparty/optimizer.h>
#include <glores/thirdparty/simplex/state.h>
#include <glores/thirdparty/simplex/criteria.h>
#include <glores/thirdparty/simplex/simplex.h>

namespace glores {

    // --------------------------------------------------------
    // BRANCH AND BOUND OPTIMIZER
    // --------------------------------------------------------

    /** 
     * Class OptimizerBB1D is a general interface for optimizing functions on 1D domain 
     * with Branch and Bound approach. 
     */
    class PointSetGlobalRegistration {
    public:

        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        /** 
         * Elements of a queue representing the distance
         */
        class PointDistance {
        public:
            //int srcId; // id of point in source point set
            int dstId; // id of point in destination point set
            double distBound; // lower bound on distance(srcId, dstId) over
            //int depth; // depth of last node who changed PointDistance
            //int intervalId; // intervalId when distBound was set
        };

        struct PointDistanceLess {

            bool operator()(const PointDistance& item1, const PointDistance& item2) const {
                return (item1.distBound < item2.distBound);
            }
        };

        struct PointDistanceMore {

            bool operator()(const PointDistance& item1, const PointDistance& item2) const {
                return (item1.distBound > item2.distBound);
            }
        };
        typedef boost::heap::binomial_heap<PointDistance, boost::heap::compare<PointDistanceMore> > PointDistanceLowerBoundMinQueue;
        typedef std::vector<PointDistance> PointDistanceQueue;

        /** 
         * Struct representing interval on domain, [xmin,xmax], and on the function values
         * given as lower and upper bounds
         */
        struct IntervalBound {
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW
            int idx;
            Eigen::VectorXd xmin;
            Eigen::VectorXd xmax;
            double ylower;
            double yupper;
            int depth; // depth in splitting operation
            int parent;
            int brother;
            bool isLeaf;
            bool relaxationBoundUsed;

            //std::vector<PointDistanceLowerBoundMinQueue> pointDistanceLowers; 
            std::vector<PointDistanceQueue> pointDistancesMin;
            std::vector<PointDistance> pointDistancesMax;
        };
        typedef std::vector<IntervalBound, Eigen::aligned_allocator<IntervalBound> > IntervalBoundVector;

        /**
         * Comparator of intervals to order interval with smaller lower bound before. 
         */
        struct IndexLowerBoundLess {
            IntervalBoundVector* nodes;

            IndexLowerBoundLess(IntervalBoundVector* n) : nodes(n) {
            }

            bool operator()(int i0, int i1) const {
                return (nodes->at(i0).ylower > nodes->at(i1).ylower);
            }
        };

        /**
         * Comparator of intervals to order interval with larger lower bound before. 
         * It is used to quickly search the bounds to be cut, i.e. the nodes
         * with lower bounds greater than global upper bound. 
         */
        struct IndexLowerBoundMore {
            IntervalBoundVector* nodes;

            IndexLowerBoundMore(IntervalBoundVector* n) : nodes(n) {
            }

            bool operator()(int i0, int i1) const {
                return (nodes->at(i0).ylower < nodes->at(i1).ylower);
            }
        };
        typedef std::priority_queue<int, std::vector<int>, IndexLowerBoundLess> LeastIndexLowerBoundFirstQueue;
        typedef boost::heap::binomial_heap<int, boost::heap::compare<IndexLowerBoundLess> > MinLowerBoundQueue;
        typedef boost::heap::binomial_heap<int, boost::heap::compare<IndexLowerBoundMore> > MaxLowerBoundQueue;

        struct SimplexObjectiveFunction {
            typedef double DataType;
            typedef Eigen::VectorXd ParameterType;
            PointSetGlobalRegistration* preg_;
            const IntervalBound* ib_;

            SimplexObjectiveFunction(PointSetGlobalRegistration* preg) : preg_(preg), ib_(0) {
            }

            SimplexObjectiveFunction(PointSetGlobalRegistration* preg, const IntervalBound* ib) : preg_(preg), ib_(ib) {
            }

            double operator()(const ParameterType& x) const {
                if (ib_ == 0) {
                    return preg_->evaluateFunction(x);
                } else {
                    return preg_->evaluateFunctionFast(x, *ib_);
                }
            }
        };

        struct BoundStats {
            int iteration;
            double lowerBound;
            double upperBound;
            bool relaxationBoundOn;
            int queueNum;
        };


        /**
         * Default constructor.
         */
        PointSetGlobalRegistration();

        /**
         * Constructor with tolerances on domain and function value. 
         */
        PointSetGlobalRegistration(double xtol, double ytol);

        /**
         * Constructor with tolerances on domain and function value. 
         */
        PointSetGlobalRegistration(const Eigen::VectorXd& xtol, double ytol);

        /**
         * Default destructor.
         */
        virtual ~PointSetGlobalRegistration();

        /**
         * 
         * @return 
         */
        int getNodeNum() const {
            return nodes_.size();
        }

        /**
         * Sets the source point set.
         */
        void setPointSource(const VectorPoint2d& pointSrc) {
            pointSrc_ = pointSrc;
        }

        /**
         * Sets the target point set.
         */
        void setPointDestination(const VectorPoint2d& pointDst) {
            pointDst_ = pointDst;
        }

        /**
         * Sets the percentage of matches to be used as inlier in the computation of bounds. 
         */
        void setInlierRatio(double ir) {
            if (ir < 0.0 || ir > 1.0) {
                GLORES_PRINT_ERROR("invalid inlier ratio " << ir << ": should be between 0.0 and 1.0: keeping old value " << inlierRatio_);
                return;
            }
            inlierRatio_ = ir;
            GLORES_PRINT_VARIABLE(inlierRatio_);
        }

        /**
         * Sets the tolerance on domain to stop estimation when an accuracy on x is reached. 
         */
        void setXTolerance(double xtol) {
            xtol_ = Eigen::VectorXd::Constant(dim_, xtol);
        }

        /**
         * Sets the tolerance on domain to stop estimation when an accuracy on x is reached. 
         */
        void setXTolerance(const Eigen::VectorXd& xtol) {
            xtol_ = xtol;
        }

        /**
         * Enables halting condition on x tolerance.
         */
        void enableXTolerance(bool xt) {
            xtolOn_ = xt;
        }

        /**
         * Enables halting condition on x tolerance.
         */
        void enableYTolerance(bool yt) {
            ytolOn_ = yt;
        }

        /**
         * Sets the tolerance on domain to stop estimation when an accuracy on y is reached. 
         */
        void setYTolerance(double ytol) {
            ytol_ = ytol;
        }
        
        void setRelaxationThres(double xThres, double yThres, double tThres) {
            relaxationThres_ << xThres, yThres, tThres;
        }

        /**
         * 
         */
        void enableRelaxationBound(bool rb) {
            relaxationBoundOn_ = rb;
        }
        
        void enableCheapBoundQueueOn(bool cbq) {
            cheapBoundQueueOn_ = cbq;
        }

        /**
         * Sets the maximum depth of exploration. 
         * @param dmax
         */
        void setDepthMax(int dmax) {
            depthMax_ = dmax;
        }

        void setRelaxationThres(const Eigen::VectorXd& relaxationThres) {
            if (relaxationThres_.rows() != 3) {
                GLORES_PRINT_ERROR("cannot set the bound for threshold relaxation: dimension " << relaxationThres_.rows() << " != 3");
                return;
            }
            relaxationThres_ = relaxationThres;
        }

        /**
         * 
         * @return 
         */
        bool setVerboseLevel(int vl) {
            verboseLevel_ = vl;
        }

        void setIterationMax(int imax) {
            iterationMax_ = imax;
        }

        /**
         * Returns a reference to the given node of BB tree. 
         */
        const IntervalBound& getInterval(int i) const {
            return nodes_[i];
        }

        /** 
         * Returns the lower and upper bound of the function over the given IntervalBound.
         */
        virtual bool findLU(const IntervalBound& ib, double& ylower, double& yupper);

        /**
         * Returns true when the stop condition is reached (e.g. sufficiently accurate solution. 
         * @param global the current estimate for solution
         */
        virtual bool stopIteration(const IntervalBound& global) const;

        /**
         * Returns a reference to the lower/upper bound statistics. 
         */
        const std::vector<BoundStats>& getBoundStats() const {
            return boundsStats_;
        }

        /**
         * Finds the global minimum point of the function on the given interval [xmin, xmax]. 
         * If multiple minima exists... 
         * @param xmin minimum value of interval domain
         * @param xmax maximum value of interval domain
         * @param x the domain point corresponding to the global optimum
         * @param ylower the lower bound of function on the given interval 
         * @param yupper the upper bound of function of the given interval
         */
        void findGlobalMin(const Eigen::VectorXd& xmin, const Eigen::VectorXd& xmax, Eigen::VectorXd& x, double& ylower, double& yupper);

        /***/
        double evaluateFunction(const Eigen::VectorXd& x) const;

        double evaluateFunctionFast(const Eigen::VectorXd& x, const IntervalBound& ib) const;

        //double evaluateFunction(const Eigen::VectorXd& x) const;


    protected:
        IntervalBoundVector nodes_;
        //std::set<int> leaves_;
        int dim_;
        VectorPoint2d pointSrc_;
        VectorPoint2d pointDst_;
        double inlierRatio_;
        Eigen::VectorXd xtol_;
        double ytol_;
        bool xtolOn_;
        bool ytolOn_;
        Eigen::VectorXd relaxationThres_;
        int depthMax_;
        int verboseLevel_;
        int iterationCounter_;
        int iterationMax_;
        bool relaxationBoundOn_;
        bool cheapBoundQueueOn_;
        double statDminQueueLenAvg_;
        int statDminQueueLenNum_;
        std::vector<BoundStats> boundsStats_;


        /**
         * Checks if the interval satisfies stop criteria. 
         * Stop criteria are based on tolerances, etc.
         */
        bool isIntervalSplitActive(const IntervalBound& ib) const;

        /**
         */
        bool checkMinMax(const Eigen::VectorXd& xmin, const Eigen::VectorXd& xmax) const;

        double intervalVolume(const Eigen::VectorXd& xmin, const Eigen::VectorXd& xmax) const;

        /**
         * Splits the interval at the given dimension d. It keeps the same lower and 
         * upper bounds ymin and ymax of the original input interval
         * @param intervCur input interval to be split
         * @param d spitting dimension
         * @param intervLow lower interval
         * @param intervUpp upper interval  
         */
        void splitInterval(const IntervalBound& intervCur, int d, IntervalBound& intervLow, IntervalBound& intervUpp) const;


        void computeAllPointDistanceLowers(IntervalBound& ib);
        
        void computeAllPointDistanceLowersNoQueue(IntervalBound& ib);

        void updateAllPointDistanceLowers(IntervalBound& ibFather, IntervalBound& ibChild, bool copyFather);

        double computeLowerBound(const IntervalBound& ib);

        double computeUpperBound(const IntervalBound& ib, bool outputOn = false);

        //double computeUpperBoundNesterMead(const IntervalBound& ib);

        double computeUpperBoundNelderMead(const IntervalBound& ib);

        double computeLowerBoundRelaxation(const IntervalBound& ib);

        void updateGlobalLowerBound(int idx, double& lowerBound);

        double findBestLeafLowerBound(int& leafLower) const;

        void printPointDistance(std::ostream& out, const IntervalBound& ib, int is) const;

        void printGeometryIntersectionParams(std::ostream& out, const IntervalBound& ib, int is, int id);

        void printPointSets(std::ostream& out, const Eigen::VectorXd& pose) const;

        double evaluateFunctionNM(nm::Vector& x) const;

        double evaluatePairwiseDistanceSquare(double tx, double ty, double ct, double st, const Point2d& p, const Point2d& q);

        void evaluatePairwiseDistanceSquareGradient(double tx, double ty, double ct, double st, const Point2d& p, const Point2d& q, Eigen::Vector4d &grad);

        double evaluateAllPointSetLinearizedDistances(double tx0, double ty0, double ct0, double st0, double tx, double ty, double ct, double st, const IntervalBound& ib);
    };

} // end of namespace

std::ostream& operator<<(std::ostream& out, const glores::PointSetGlobalRegistration::IntervalBound& interv);

#endif /* OPTIMIZERBBND_H */


