/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
#include <glores/geometry.h>
#include <glores/fileutils.h>
#include <tuple>

namespace glores {

    // --------------------------------------------------------------
    // TYPE CONVERSION
    // --------------------------------------------------------------

    void poseToIsometry2d(double x, double y, double theta, glores::Transformation2d& transf) {
        transf = Transformation2d::Identity();
        transf.prerotate(Eigen::Rotation2Dd(theta));
        transf.pretranslate(Eigen::Vector2d(x, y));
    }

    void poseToIsometry2f(float x, float y, float theta, glores::Transformation2f& transf) {
        transf = Transformation2f::Identity();
        transf.prerotate(Eigen::Rotation2Df(theta));
        transf.pretranslate(Eigen::Vector2f(x, y));
    }

    void poseToIsometry2d(const Eigen::Vector3d& pose, glores::Transformation2d& transf) {
        poseToIsometry2d(pose(0), pose(1), pose(2), transf);
    }

    void poseToIsometry2f(const Eigen::Vector3f& pose, glores::Transformation2f& transf) {
        poseToIsometry2f(pose(0), pose(1), pose(2), transf);
    }

    // --------------------------------------------------------------
    // FUNCTIONS
    // --------------------------------------------------------------

    double normalizeAngle2Pi(double angle) {
        static const double K_2PI = 2.0 * M_PI;
        return (angle - K_2PI * floor(angle / K_2PI));
    }

    bool insideAngleInterval(double x, double y, Point2d pBeg, Point2d pEnd) {
        double scalarE = -pEnd(0) * y + pEnd(1) * x;
        double scalarB = pBeg(0) * y - pBeg(1) * x;
        int reverseCone = 0;
        if (pBeg(0) * pEnd(1) - pEnd(0) * pBeg(1) < 0) {
            reverseCone = 1;
            scalarE = -scalarE;
            scalarB = -scalarB;
        }

        if (!reverseCone) {
            return ((scalarE > 0)&&(scalarB > 0));
        } else {//reverse cone
            return (!((scalarE > 0)&&(scalarB > 0)));
        }
        assert(0);
    }

    bool insideAngleInterval(double angle, double intervBeg, double intervEnd) {
        return (normalizeAngle2Pi(angle - intervBeg) < normalizeAngle2Pi(intervEnd - intervBeg));
    }

    bool insideAngleInterval(double angle, const AngleInterval& interv) {
        return insideAngleInterval(angle, interv.first, interv.second);
    }

    bool insideAngleIntervalEq(double angle, double intervBeg, double intervEnd) {
        return (normalizeAngle2Pi(angle - intervBeg) <= normalizeAngle2Pi(intervEnd - intervBeg));
    }

    int intersectAngleInterval(const AngleInterval& in1, const AngleInterval& in2, std::vector<AngleInterval>& out) {
        out.clear();

        // Intersection of angular intervals [b1, e1] and [b2, e2] where 
        // e1 is reached CCW from b1 and e2 is reached CCW from b2. 
        // Criterion based on this trivial property: for any intersection interval
        // [bi, ei], bi is either b1 or b2.

        // Case 1: b2 lies inside [b1,e1]. Then, b2 is the begin of an angular 
        // intersection, either [b2, e1] or [b2, e2]. 
        if (insideAngleInterval(in2.first, in1)) {
            if (insideAngleInterval(in1.second, in2)) {
                out.push_back(std::make_pair(in2.first, in1.second));
            } else if (insideAngleInterval(in2.second, in1)) {
                out.push_back(in2);
            }
        }

        // Case 2: b1 lies inside [b2,e2]. Then, b1 is the begin of an angular 
        // intersection, either [b1, e1] or [b1, e2]. 
        if (insideAngleInterval(in1.first, in2)) {
            if (insideAngleInterval(in2.second, in1)) {
                out.push_back(std::make_pair(in1.first, in2.second));
            } else if (insideAngleInterval(in1.second, in2)) {
                out.push_back(in1);
            }
        }
        return out.size();
    }

    void projectPointOnLine(const Point2d& p, double theta, double rho, Point2d& q) {
        Point2d u1, u2;
        if (rho < 0.0) {
            theta = theta + M_PI;
            rho = -rho;
        }
        u1 << cos(theta), sin(theta);
        u2 << -sin(theta), cos(theta);
        q = rho * u1 + u2.dot(p) * u2;
    }

    void computeLine(const Point2d& p1, const Point2d& p2, double& theta, double& rho) {
        Point2d d12;
        computeLine(p1, p2, theta, rho, d12);
    }

    void computeLine(const Point2d& p1, const Point2d& p2, double& theta, double& rho, Point2d& d12) {
        d12 = p2 - p1;
        theta = fmod(atan2(d12.x(), -d12.y()) + 2.0 * M_PI, 2 * M_PI);
        rho = p1.x() * cos(theta) + p1.y() * sin(theta);

        //        GLORES_PRINT_VARIABLE(p1.transpose());
        //        GLORES_PRINT_VARIABLE(p2.transpose());
        //        GLORES_PRINT_VARIABLE(d12.transpose());
        //        GLORES_PRINT_VARIABLE(180.0 / M_PI * theta);
        //        GLORES_PRINT_VARIABLE(rho);
        //        GLORES_PRINT_MSG("half-plane: x * " << cos(theta) << " + y * " << sin(theta) << " - (" << rho << ") >= 0.0\n");
    }

    bool projectArcBound(double arcBeg, double arcEnd, VectorPoint2d& polygonBound, int vnum) {
        Point2d v;
        double dtheta, radius, a;

        if (arcBeg < 0.0 || arcBeg > 2.0 * M_PI) {
            std::cerr << __FILE__ << "," << __LINE__ << ": invalid arcBeg " << (180.0 / M_PI * arcBeg) << " [deg] for angle: must be in [0.0, 360.0[" << std::endl;
            return false;
        }
        if (arcEnd < 0.0 || arcEnd > 2.0 * M_PI) {
            std::cerr << __FILE__ << "," << __LINE__ << ": invalid arcEnd " << (180.0 / M_PI * arcEnd) << " [deg] for angle: must be in [0.0, 360.0[" << std::endl;
            return false;
        }
        if (arcBeg > arcEnd) {
            GLORES_PRINT_ERROR("invalid arc bounds [" << (180.0 / M_PI * arcBeg) << "," << (180.0 / M_PI * arcEnd) << "]");
            return false;
        }
        if (vnum < 1) {
            GLORES_PRINT_ERROR("arc vertices number vnum " << vnum << " must be >= 1");
        }

        dtheta = (arcEnd - arcBeg) / (2.0 * vnum);
        radius = 1 / cos(dtheta);

        v << cos(arcBeg), sin(arcBeg);
        polygonBound.push_back(v);
        for (int i = 0; i < vnum; ++i) {
            a = arcBeg + (2 * i + 1) * dtheta;
            v(0) = radius * cos(a);
            v(1) = radius * sin(a);
            polygonBound.push_back(v);
        }
        v << cos(arcEnd), sin(arcEnd);
        polygonBound.push_back(v);

        return true;
    }

    bool projectSE2Interval(const Point2d& p, const Eigen::Vector3d& poseLow, const Eigen::Vector3d& poseUpp, Polygon2d& convexBound, int vnum) {
        double txLow, txUpp, tyLow, tyUpp, thetaLow, thetaUpp;
        double dtheta, beta, radius, a;
        Polygon2d vertices;
        Point2d v, vLL, vLU, vUL, vUU;

        //std::cout << "thetaLow[deg] " << (180.0 / M_PI * poseLow(2)) << ", thetaUpp[deg] " << (180.0 / M_PI * poseUpp(2)) << std::endl;

        // Normalizes angles of poses into interval [0.0, 2*M_PI[
        txLow = poseLow(0);
        tyLow = poseLow(1);
        thetaLow = poseLow(2);
        // Check intervals of rigid motion bounds
        if (poseLow(2) < 0.0 || poseLow(2) > 2.0 * M_PI) {
            std::cerr << __FILE__ << "," << __LINE__ << ": invalid lower bound " << (180.0 / M_PI * poseLow(2)) << " [deg] for angle: must be in [0.0, 360.0[" << std::endl;
            return false;
        }
        txUpp = poseUpp(0);
        tyUpp = poseUpp(1);
        thetaUpp = poseUpp(2);
        if (poseUpp(2) < 0.0 || poseUpp(2) > 2.0 * M_PI) {
            std::cerr << __FILE__ << "," << __LINE__ << ": invalid upper bound " << (180.0 / M_PI * poseUpp(2)) << " [deg] for angle: must be in [0.0, 360.0[" << std::endl;
            return false;
        }
        if (poseLow(0) > poseUpp(0) || poseLow(1) > poseUpp(1) || poseLow(2) > poseUpp(2)) {
            std::cerr << __FILE__ << "," << __LINE__ << ": invalid bounds\n"
                    << "  x [" << poseLow(0) << "," << poseUpp(0) << "]\n"
                    << "  y [" << poseLow(1) << "," << poseUpp(1) << "]\n"
                    << "  theta [" << (180.0 / M_PI * poseLow(2)) << "," << (180.0 / M_PI * poseUpp(2)) << "] deg\n"
                    << std::endl;
            return false;
        }
        //std::cout << "txLow " << txLow << " txUpp " << txUpp << " tyLow "  << tyLow << " tyUpp " << tyUpp << std::endl;
        //std::cout << "thetaLow[deg] " << (180.0 / M_PI * thetaLow) << ", thetaUpp[deg] " << (180.0 / M_PI * thetaUpp) << std::endl;

        // Computes the vertives of polygon enclosing polygon
        //        beta = atan2(p.y(), p.x());
        //        if (beta < 0.0) {
        //            beta += 2.0 * M_PI;
        //        }
        dtheta = (thetaUpp - thetaLow) / (2.0 * vnum);
        radius = 1 / cos(dtheta);
        //std::cout << "dtheta[deg] " << (180.0 / M_PI * dtheta) << std::endl;
        v(0) = p(0) * cos(thetaLow) - p(1) * sin(thetaLow);
        v(1) = p(0) * sin(thetaLow) + p(1) * cos(thetaLow);
        vLL(0) = v(0) + txLow;
        vLL(1) = v(1) + tyLow;
        boost::geometry::append(vertices, vLL);
        vLU(0) = v(0) + txLow;
        vLU(1) = v(1) + tyUpp;
        boost::geometry::append(vertices, vLU);
        vUL(0) = v(0) + txUpp;
        vUL(1) = v(1) + tyLow;
        boost::geometry::append(vertices, vUL);
        vUU(0) = v(0) + txUpp;
        vUU(1) = v(1) + tyUpp;
        boost::geometry::append(vertices, vUU);
        for (int i = 0; i < vnum; ++i) {
            //std::cout << "angle " << i << ": a[deg] " << (180.0 / M_PI * a) << std::endl;
            a = thetaLow + (2 * i + 1) * dtheta;
            v(0) = radius * (p(0) * cos(a) - p(1) * sin(a));
            v(1) = radius * (p(0) * sin(a) + p(1) * cos(a));
            vLL(0) = v(0) + txLow;
            vLL(1) = v(1) + tyLow;
            boost::geometry::append(vertices, vLL);
            vLU(0) = v(0) + txLow;
            vLU(1) = v(1) + tyUpp;
            boost::geometry::append(vertices, vLU);
            vUL(0) = v(0) + txUpp;
            vUL(1) = v(1) + tyLow;
            boost::geometry::append(vertices, vUL);
            vUU(0) = v(0) + txUpp;
            vUU(1) = v(1) + tyUpp;
            boost::geometry::append(vertices, vUU);
        }
        v(0) = p(0) * cos(thetaUpp) - p(1) * sin(thetaUpp);
        v(1) = p(0) * sin(thetaUpp) + p(1) * cos(thetaUpp);
        vLL(0) = v(0) + txLow;
        vLL(1) = v(1) + tyLow;
        boost::geometry::append(vertices, vLL);
        vLU(0) = v(0) + txLow;
        vLU(1) = v(1) + tyUpp;
        boost::geometry::append(vertices, vLU);
        vUL(0) = v(0) + txUpp;
        vUL(1) = v(1) + tyLow;
        boost::geometry::append(vertices, vUL);
        vUU(0) = v(0) + txUpp;
        vUU(1) = v(1) + tyUpp;
        boost::geometry::append(vertices, vUU);

        // Computes the convex hull 
        boost::geometry::convex_hull(vertices, convexBound);

        //        std::cout << __FILE__ << "," << __LINE__ << ": vertices:\n";
        //        for (auto pit = boost::begin(boost::geometry::exterior_ring(vertices)); pit != boost::end(boost::geometry::exterior_ring(vertices)); ++pit) {
        //            //            double x = boost::geometry::get<0>(*pit);
        //            //            double y = boost::geometry::get<1>(*pit);
        //            //            std::cout << "  " << x << " " << y << "\n";
        //            std::cout << "  " << pit->transpose() << "\n";
        //        }
        //        std::cout << __FILE__ << "," << __LINE__ << ": convex hull:\n";
        //        for (auto pit = boost::begin(boost::geometry::exterior_ring(convexBound)); pit != boost::end(boost::geometry::exterior_ring(convexBound)); ++pit) {
        //            std::cout << "  " << pit->transpose() << "\n";
        //        }
    }

    bool computeDistanceMinMax(const Point2d& pSrc, const Point2d& pDst, const Eigen::Vector3d& poseLow, const Eigen::Vector3d& poseUpp, double& dmin, double& dmax) {
        Point2d pMin, pMax;
        VectorPoint2d keypoints;
        pMin << 0.0, 0.0;
        pMax << 0.0, 0.0;
        return computeDistanceMinMax(pSrc, pDst, poseLow, poseUpp, dmin, dmax, pMin, pMax, keypoints);
    }

    bool computeDistanceMinMax(const Point2d& pSrc, const Point2d& pDst, const Eigen::Vector3d& poseLow, const Eigen::Vector3d& poseUpp, double& dmin, double& dmax, Point2d& pMin, Point2d& pMax, VectorPoint2d& keypoints) {
        Point2d v, nearestLine, nearestArc;
        double dist;
        double radius, beta, arcBeg, arcEnd;
        double xl, xu, yl, yu, thetaXL, thetaXU, thetaYL, thetaYU;

        // Negative values of dmin and dmax to signal to the user that this
        // function has been aborted. 
        dmin = -1.0;
        dmax = -1.0;

        // Check intervals of rigid motion bounds
        if (poseLow(2) < 0.0 || poseLow(2) > 2.0 * M_PI - 1e-5) {
            std::cerr << __FILE__ << "," << __LINE__ << ": invalid lower bound " << (180.0 / M_PI * poseLow(2)) << " [deg] for angle: must be in [0.0, 360.0[" << std::endl;
            return false;
        }
        if (poseUpp(2) < 0.0 || poseUpp(2) > 2.0 * M_PI - 1e-5) {
            std::cerr << __FILE__ << "," << __LINE__ << ": invalid upper bound " << (180.0 / M_PI * poseUpp(2)) << " [deg] for angle: must be in [0.0, 360.0[" << std::endl;
            return false;
        }
        if (poseLow(0) > poseUpp(0) || poseLow(1) > poseUpp(1) || poseLow(2) > poseUpp(2)) {
            std::cerr << __FILE__ << "," << __LINE__ << ": invalid bounds\n"
                    << "  x [" << poseLow(0) << "," << poseUpp(0) << "]\n"
                    << "  y [" << poseLow(1) << "," << poseUpp(1) << "]\n"
                    << "  theta [" << (180.0 / M_PI * poseLow(2)) << "," << (180.0 / M_PI * poseUpp(2)) << "] deg\n"
                    << std::endl;
            return false;
        }

        // Distance min and max is given by the distance of arc and four vertices:
        // - the arc is between Rot(poseLow(2))*pSrc and Rot(poseUpp(2))*pSrc
        // - the vertices are pDst -/+ [poseLow(0) : poseUpp(0), poseLow(1) : poseUpp(1)]
        radius = pSrc.norm();
        beta = atan2(pSrc.y(), pSrc.x());
        arcBeg = beta + poseLow(2);
        arcEnd = beta + poseUpp(2);
        xl = pDst(0) - poseUpp(0);
        xu = pDst(0) - poseLow(0);
        yl = pDst(1) - poseUpp(1);
        yu = pDst(1) - poseLow(1);
        //        thetaXL = (xl > 0.0 ? 0.0 : M_PI);
        //        thetaXU = (xu > 0.0 ? M_PI : 0.0);
        //        thetaYL = (yl > 0.0 ? 0.5 * M_PI : 1.5 * M_PI);
        //        thetaYL = (yu > 0.0 ? 1.5 * M_PI : 0.5 * M_PI);
        thetaXL = 0.0;
        thetaXU = M_PI;
        thetaYL = 0.5 * M_PI;
        thetaYU = 1.5 * M_PI;

        //        GLORES_PRINT_VARIABLE(xl);
        //        GLORES_PRINT_VARIABLE(xu);
        //        GLORES_PRINT_VARIABLE(yl);
        //        GLORES_PRINT_VARIABLE(yu);
        //        GLORES_PRINT_VARIABLE(radius);
        //        GLORES_PRINT_VARIABLE(180.0 / M_PI * beta);
        //        GLORES_PRINT_VARIABLE(180.0 / M_PI * arcBeg);
        //        GLORES_PRINT_VARIABLE(180.0 / M_PI * arcEnd);
        //        GLORES_PRINT_VARIABLE(180 / M_PI * thetaXL);
        //        GLORES_PRINT_VARIABLE(180 / M_PI * thetaXU);
        //        GLORES_PRINT_VARIABLE(180 / M_PI * thetaYL);
        //        GLORES_PRINT_VARIABLE(180 / M_PI * thetaYU);

        keypoints.clear();

        // Inserts box vertices in keypoints (points where to evaluate min distance)
        v << xl, yl;
        keypoints.push_back(v);
        v << xu, yl;
        keypoints.push_back(v);
        v << xu, yu;
        keypoints.push_back(v);
        v << xl, yu;
        keypoints.push_back(v);

        // Evaluates intersection with the half-planes defined by box edges
        bool intersectionArcBox = false;
        bool inHalfPlaneXL = glores::intersectCircularArcHalfPlane(thetaXL, xl, radius, arcBeg, arcEnd, nearestLine, nearestArc);
        if (!inHalfPlaneXL && yl <= nearestLine.y() && nearestLine.y() <= yu) {
            keypoints.push_back(nearestLine);
        }
        bool inHalfPlaneXU = glores::intersectCircularArcHalfPlane(thetaXU, -xu, radius, arcBeg, arcEnd, nearestLine, nearestArc);
        if (!inHalfPlaneXU && yl <= nearestLine.y() && nearestLine.y() <= yu) {
            keypoints.push_back(nearestLine);
        }
        bool inHalfPlaneYL = glores::intersectCircularArcHalfPlane(thetaYL, yl, radius, arcBeg, arcEnd, nearestLine, nearestArc);
        if (!inHalfPlaneYL && xl <= nearestLine.x() && nearestLine.x() <= xu) {
            keypoints.push_back(nearestLine);
        }
        bool inHalfPlaneYU = glores::intersectCircularArcHalfPlane(thetaYU, -yu, radius, arcBeg, arcEnd, nearestLine, nearestArc);
        if (!inHalfPlaneYU && xl <= nearestLine.x() && nearestLine.x() <= xu) {
            keypoints.push_back(nearestLine);
        }
        //        GLORES_PRINT_VARIABLE(inHalfPlaneXL);
        //        GLORES_PRINT_VARIABLE(inHalfPlaneXU);
        //        GLORES_PRINT_VARIABLE(inHalfPlaneYL);
        //        GLORES_PRINT_VARIABLE(inHalfPlaneYU);

        // If there is an intersection, i.e. there is a portion of arc in all 
        // the half-planes, then dmin = 0.0
        if (inHalfPlaneXL && inHalfPlaneXU && inHalfPlaneYL && inHalfPlaneYU) {
            dmin = 0.0;
            dmax = 0.0;
            for (auto& kp : keypoints) {
                dist = computeDistancePointArc(kp, radius, arcBeg, arcEnd, pMin);
                if (dist > dmax) {
                    dmax = dist;
                    pMax = kp;
                }
            }
        } else {
            // Checks minimum distance from comparison with keypoints
            dmin = 1e+8;
            dmax = 0.0;
            for (auto& kp : keypoints) {
                dist = computeDistancePointArc(kp, radius, arcBeg, arcEnd);
                if (dist < dmin) {
                    dmin = dist;
                    pMin = kp;
                }
                if (dist > dmax) {
                    dmax = dist;
                    pMax = kp;
                }
            }
        }

        return true;
    }

    bool computeDistanceMinMax2(const Point2d& pSrc, const Point2d& pDst, const Eigen::Vector3d& poseLow, const Eigen::Vector3d& poseUpp, double& dmin, double& dmax) {
        Point2d pMin, pMax;
        VectorPoint2d keypoints;
        pMin << 0.0, 0.0;
        pMax << 0.0, 0.0;
        return computeDistanceMinMax2(pSrc, pDst, poseLow, poseUpp, dmin, dmax, pMin, pMax, keypoints);
    }

    bool computeDistanceMinMax2(const Point2d& pSrc, const Point2d& pDst, const Eigen::Vector3d& poseLow, const Eigen::Vector3d& poseUpp, double& dmin, double& dmax, Point2d& pMin, Point2d& pMax, VectorPoint2d& keypoints) {
        Point2d v1, v2, v3, v4, nearestLine, nearestArc;
        double dist;
        double radius, beta, arcBeg, arcEnd;
        double xl, xu, yl, yu, thetaXL, thetaXU, thetaYL, thetaYU;

        // Negative values of dmin and dmax to signal to the user that this
        // function has been aborted. 
        dmin = -1.0;
        dmax = -1.0;

        // Check intervals of rigid motion bounds
        if (poseLow(2) < 0.0 || poseLow(2) > 2.0 * M_PI - 1e-5) {
            std::cerr << __FILE__ << "," << __LINE__ << ": invalid lower bound " << (180.0 / M_PI * poseLow(2)) << " [deg] for angle: must be in [0.0, 360.0[" << std::endl;
            return false;
        }
        if (poseUpp(2) < 0.0 || poseUpp(2) > 2.0 * M_PI - 1e-5) {
            std::cerr << __FILE__ << "," << __LINE__ << ": invalid upper bound " << (180.0 / M_PI * poseUpp(2)) << " [deg] for angle: must be in [0.0, 360.0[" << std::endl;
            return false;
        }
        if (poseLow(0) > poseUpp(0) || poseLow(1) > poseUpp(1) || poseLow(2) > poseUpp(2)) {
            std::cerr << __FILE__ << "," << __LINE__ << ": invalid bounds\n"
                    << "  x [" << poseLow(0) << "," << poseUpp(0) << "]\n"
                    << "  y [" << poseLow(1) << "," << poseUpp(1) << "]\n"
                    << "  theta [" << (180.0 / M_PI * poseLow(2)) << "," << (180.0 / M_PI * poseUpp(2)) << "] deg\n"
                    << std::endl;
            return false;
        }

        // Distance min and max is given by the distance of arc and four vertices:
        // - the arc is between Rot(poseLow(2))*pSrc and Rot(poseUpp(2))*pSrc
        // - the vertices are pDst -/+ [poseLow(0) : poseUpp(0), poseLow(1) : poseUpp(1)]
        radius = pSrc.norm();
        beta = atan2(pSrc.y(), pSrc.x());
        arcBeg = beta + poseLow(2);
        arcEnd = beta + poseUpp(2);
        xl = pDst(0) - poseUpp(0);
        xu = pDst(0) - poseLow(0);
        yl = pDst(1) - poseUpp(1);
        yu = pDst(1) - poseLow(1);
        //
        //        GLORES_PRINT_VARIABLE(xl);
        //        GLORES_PRINT_VARIABLE(xu);
        //        GLORES_PRINT_VARIABLE(yl);
        //        GLORES_PRINT_VARIABLE(yu);
        //        GLORES_PRINT_VARIABLE(radius);
        //        GLORES_PRINT_VARIABLE(180.0 / M_PI * beta);
        //        GLORES_PRINT_VARIABLE(180.0 / M_PI * arcBeg);
        //        GLORES_PRINT_VARIABLE(180.0 / M_PI * arcEnd);

        // Inserts box vertices in keypoints (points where to evaluate min distance)
        keypoints.clear();
        v1 << xl, yl;
        keypoints.push_back(v1);
        v2 << xu, yl;
        keypoints.push_back(v2);
        v3 << xu, yu;
        keypoints.push_back(v3);
        v4 << xl, yu;
        keypoints.push_back(v4);

        // Evaluates intersection with the half-planes defined by box edges
        bool inHalfPlaneXL = glores::intersectCircularArcHalfPlaneSegment(v1, v2, radius, arcBeg, arcEnd, nearestLine, nearestArc);
        if (!inHalfPlaneXL) { // && yl <= nearestLine.y() && nearestLine.y() <= yu
            keypoints.push_back(nearestLine);
        }
        bool inHalfPlaneYU = glores::intersectCircularArcHalfPlaneSegment(v2, v3, radius, arcBeg, arcEnd, nearestLine, nearestArc);
        if (!inHalfPlaneYU) { // && xl <= nearestLine.x() && nearestLine.x() <= xu
            keypoints.push_back(nearestLine);
        }
        bool inHalfPlaneXU = glores::intersectCircularArcHalfPlaneSegment(v3, v4, radius, arcBeg, arcEnd, nearestLine, nearestArc);
        if (!inHalfPlaneXU) { // && yl <= nearestLine.y() && nearestLine.y() <= yu
            keypoints.push_back(nearestLine);
        }
        bool inHalfPlaneYL = glores::intersectCircularArcHalfPlaneSegment(v4, v1, radius, arcBeg, arcEnd, nearestLine, nearestArc);
        if (!inHalfPlaneYL) { // && xl <= nearestLine.x() && nearestLine.x() <= xu
            keypoints.push_back(nearestLine);
        }

        //        GLORES_PRINT_VARIABLE(inHalfPlaneXL);
        //        GLORES_PRINT_VARIABLE(inHalfPlaneYU);
        //        GLORES_PRINT_VARIABLE(inHalfPlaneXU);
        //        GLORES_PRINT_VARIABLE(inHalfPlaneYL);


        // If there is an intersection, i.e. there is a portion of arc in all 
        // the half-planes, then dmin = 0.0
        dmax = 0.0;
        for (auto& kp : keypoints) {
            dist = computeMaxDistancePointArc(kp, radius, arcBeg, arcEnd);
            if (dist > dmax) {
                dmax = dist;
                pMax = kp;
            }
        }

        if (inHalfPlaneXL && inHalfPlaneXU && inHalfPlaneYL && inHalfPlaneYU) {
            dmin = 0.0;
        } else {
            // Checks minimum distance from comparison with keypoints
            dmin = 1e+8;
            for (auto& kp : keypoints) {
                dist = computeDistancePointArc(kp, radius, arcBeg, arcEnd);
                if (dist < dmin) {
                    dmin = dist;
                    pMin = kp;
                }
            }
        }

        return true;
    }

    bool computeDistanceMinMax3(const Point2d& pArcBeg, const Point2d& pArcEnd, const double beta, double radius, const Point2d& pDst, const Eigen::Vector3d& poseLow, const Eigen::Vector3d& poseUpp, double& dmin, double& dmax) {
        Point2d v1, v2, v3, v4, pArcMid;
        double dist;
        double arcBeg, arcEnd;
        double xl, xu, yl, yu;

        // Negative values of dmin and dmax to signal to the user that this
        // function has been aborted.
        dmin = -1.0;
        dmax = -1.0;

        // Check intervals of rigid motion bounds
        if (poseLow(2) < 0.0 || poseLow(2) > 2.0 * M_PI) {
            std::cerr << __FILE__ << "," << __LINE__ << ": invalid lower bound " << (180.0 / M_PI * poseLow(2)) << " [deg] for angle: must be in [0.0, 360.0[" << std::endl;
            return false;
        }
        if (poseUpp(2) < 0.0 || poseUpp(2) > 2.0 * M_PI) {
            std::cerr << __FILE__ << "," << __LINE__ << ": invalid upper bound " << (180.0 / M_PI * poseUpp(2)) << " [deg] for angle: must be in [0.0, 360.0[" << std::endl;
            return false;
        }
        if (poseLow(0) > poseUpp(0) || poseLow(1) > poseUpp(1) || poseLow(2) > poseUpp(2)) {
            std::cerr << __FILE__ << "," << __LINE__ << ": invalid bounds\n"
                    << "  x [" << poseLow(0) << "," << poseUpp(0) << "]\n"
                    << "  y [" << poseLow(1) << "," << poseUpp(1) << "]\n"
                    << "  theta [" << (180.0 / M_PI * poseLow(2)) << "," << (180.0 / M_PI * poseUpp(2)) << "] deg\n"
                    << std::endl;
            return false;
        }

        // Distance min and max is given by the distance of arc and four vertices:
        // - the arc is between Rot(poseLow(2))*pSrc and Rot(poseUpp(2))*pSrc
        // - the vertices are pDst -/+ [poseLow(0) : poseUpp(0), poseLow(1) : poseUpp(1)]
        xl = pDst(0) - poseUpp(0);
        xu = pDst(0) - poseLow(0);
        yl = pDst(1) - poseUpp(1);
        yu = pDst(1) - poseLow(1);
        v1 << xl, yl;
        v2 << xu, yl;
        v3 << xu, yu;
        v4 << xl, yu;
        arcBeg = beta + poseLow(2);
        arcEnd = beta + poseUpp(2);

        //        GLORES_PRINT_VARIABLE(xl);
        //        GLORES_PRINT_VARIABLE(xu);
        //        GLORES_PRINT_VARIABLE(yl);
        //        GLORES_PRINT_VARIABLE(yu);
        //        GLORES_PRINT_VARIABLE(pArcBeg.transpose());
        //        GLORES_PRINT_VARIABLE(pArcEnd.transpose());

        // Computing dmax
        dist = computeMaxDistancePointArc(v1, radius, pArcBeg, pArcEnd);
        if (dist > dmax) {
            dmax = dist;
        }
        dist = computeMaxDistancePointArc(v2, radius, pArcBeg, pArcEnd);
        if (dist > dmax) {

            dmax = dist;
        }
        dist = computeMaxDistancePointArc(v3, radius, pArcBeg, pArcEnd);
        if (dist > dmax) {
            dmax = dist;
        }
        dist = computeMaxDistancePointArc(v4, radius, pArcBeg, pArcEnd);
        if (dist > dmax) {
            dmax = dist;
        }

        // Computing dmin:
        // Case 1: one arc point is inside the rectangle v1,v2,v3,v4
        // Case 2: test segment arc intersection
        dmin = 1e+6;
        if ((xl <= pArcBeg(0) && pArcBeg(0) <= xu && yl <= pArcBeg(1) && pArcBeg(1) <= yu) ||
                (xl <= pArcEnd(0) && pArcEnd(0) <= xu && yl <= pArcEnd(1) && pArcEnd(1) <= yu) ||
                intersectSegmentArcHor(xl, xu, yl, radius, pArcBeg, pArcEnd) ||
                intersectSegmentArcHor(xl, xu, yu, radius, pArcBeg, pArcEnd) ||
                intersectSegmentArcVer(xl, yl, yu, radius, pArcBeg, pArcEnd) ||
                intersectSegmentArcVer(xu, yl, yu, radius, pArcBeg, pArcEnd)) {
            dmin = 0.0;
        } else {
            // Checks all the distances between arc and rectangle vertices
            dist = computeDistancePointArc(v1, radius, pArcBeg, pArcEnd);
            //            GLORES_PRINT_VARIABLE(computeDistancePointArc(v1, radius, arcBeg, arcEnd));
            if (dist < dmin) {
                dmin = dist;
            }
            dist = computeDistancePointArc(v2, radius, pArcBeg, pArcEnd);
            //            GLORES_PRINT_VARIABLE(computeDistancePointArc(v2, radius, arcBeg, arcEnd));
            if (dist < dmin) {
                dmin = dist;
            }
            dist = computeDistancePointArc(v3, radius, pArcBeg, pArcEnd);
            //            GLORES_PRINT_VARIABLE(computeDistancePointArc(v3, radius, arcBeg, arcEnd));
            if (dist < dmin) {
                dmin = dist;
            }
            dist = computeDistancePointArc(v4, radius, pArcBeg, pArcEnd);
            //            GLORES_PRINT_VARIABLE(computeDistancePointArc(v4, radius, arcBeg, arcEnd));
            if (dist < dmin) {
                dmin = dist;
            }
            // Checks all the distances between each edge (segments) of rectangle and pArcBeg
            //            GLORES_PRINT_VARIABLE(computeDistancePointSegmentHor(pArcBeg, xl, xu, yl));
            dist = computeDistancePointSegmentHor(pArcBeg, xl, xu, yl);
            if (dist < dmin) {
                dmin = dist;
            }
            //            GLORES_PRINT_VARIABLE(computeDistancePointSegmentHor(pArcBeg, xl, xu, yu));
            dist = computeDistancePointSegmentHor(pArcBeg, xl, xu, yu);
            if (dist < dmin) {
                dmin = dist;
            }
            //            GLORES_PRINT_VARIABLE(computeDistancePointSegmentVer(pArcBeg, xl, yl, yu));
            dist = computeDistancePointSegmentVer(pArcBeg, xl, yl, yu);
            if (dist < dmin) {
                dmin = dist;
            }
            //            GLORES_PRINT_VARIABLE(computeDistancePointSegmentVer(pArcBeg, xu, yl, yu));
            dist = computeDistancePointSegmentVer(pArcBeg, xu, yl, yu);
            if (dist < dmin) {
                dmin = dist;
            }
            // Checks all the distances between each edge (segments) of rectangle and pArcEnd
            //            GLORES_PRINT_VARIABLE(computeDistancePointSegmentHor(pArcEnd, xl, xu, yl));
            dist = computeDistancePointSegmentHor(pArcEnd, xl, xu, yl);
            if (dist < dmin) {
                dmin = dist;
            }
            //            GLORES_PRINT_VARIABLE(computeDistancePointSegmentHor(pArcEnd, xl, xu, yu));
            dist = computeDistancePointSegmentHor(pArcEnd, xl, xu, yu);
            if (dist < dmin) {
                dmin = dist;
            }
            //            GLORES_PRINT_VARIABLE(computeDistancePointSegmentVer(pArcEnd, xl, yl, yu));
            dist = computeDistancePointSegmentVer(pArcEnd, xl, yl, yu);
            if (dist < dmin) {
                dmin = dist;
            }
            //            GLORES_PRINT_VARIABLE(computeDistancePointSegmentVer(pArcEnd, xu, yl, yu));
            dist = computeDistancePointSegmentVer(pArcEnd, xu, yl, yu);
            if (dist < dmin) {
                dmin = dist;
            }
            // Checks all the distances between
            //if (insideAngleInterval(0.0, arcBeg, arcEnd) && (yl <= 0.0 && 0 <= yu)) {
            if (insideAngleInterval(1.0, 0.0, pArcBeg, pArcEnd) && (yl <= 0.0 && 0.0 <= yu)) {
                dist = fabs(radius - xl);
                if (dist < dmin) {
                    dmin = dist;
                }
                dist = fabs(radius - xu);
                if (dist < dmin) {
                    dmin = dist;
                }
            }
            //if (insideAngleInterval(M_PI, arcBeg, arcEnd) && (yl <= 0.0 && 0 <= yu)) {
            if (insideAngleInterval(-1.0, 0.0, pArcBeg, pArcEnd) && (yl <= 0.0 && 0.0 <= yu)) {
                dist = fabs(-radius - xl);
                if (dist < dmin) {
                    dmin = dist;
                }
                dist = fabs(-radius - xu);
                if (dist < dmin) {
                    dmin = dist;
                }
            }
            // Checks all the distances between

            //if (insideAngleInterval(0.5 * M_PI, arcBeg, arcEnd) && (xl <= 0.0 && 0 <= xu)) {
            if (insideAngleInterval(0.0, 1.0, pArcBeg, pArcEnd) && (xl <= 0.0 && 0.0 <= xu)) {
                dist = fabs(radius - yl);
                if (dist < dmin) {
                    dmin = dist;
                }
                dist = fabs(radius - yu);
                if (dist < dmin) {
                    dmin = dist;
                }
            }
            //if (insideAngleInterval(1.5 * M_PI, arcBeg, arcEnd) && (xl <= 0.0 && 0 <= xu)) {
            if (insideAngleInterval(0.0, -1.0, pArcBeg, pArcEnd) && (xl <= 0.0 && 0.0 <= xu)) {
                dist = fabs(-radius - yl);
                if (dist < dmin) {
                    dmin = dist;
                }
                dist = fabs(-radius - yu);
                if (dist < dmin) {
                    dmin = dist;
                }
            }
        }
        return true;
    }

    bool computeDistanceMinMax3(const Point2d& pSrc, const Point2d& pDst, const Eigen::Vector3d& poseLow, const Eigen::Vector3d& poseUpp, double& dmin, double& dmax) {
        Point2d v1, v2, v3, v4, pArcBeg, pArcEnd, pArcMid;
        double dist;
        double radius, beta, arcBeg, arcEnd;
        double xl, xu, yl, yu;

        // Negative values of dmin and dmax to signal to the user that this
        // function has been aborted. 
        dmin = -1.0;
        dmax = -1.0;

        // Check intervals of rigid motion bounds
        if (poseLow(2) < 0.0 || poseLow(2) > 2.0 * M_PI) {
            std::cerr << __FILE__ << "," << __LINE__ << ": invalid lower bound " << (180.0 / M_PI * poseLow(2)) << " [deg] for angle: must be in [0.0, 360.0[" << std::endl;
            return false;
        }
        if (poseUpp(2) < 0.0 || poseUpp(2) > 2.0 * M_PI) {
            std::cerr << __FILE__ << "," << __LINE__ << ": invalid upper bound " << (180.0 / M_PI * poseUpp(2)) << " [deg] for angle: must be in [0.0, 360.0[" << std::endl;
            return false;
        }
        if (poseLow(0) > poseUpp(0) || poseLow(1) > poseUpp(1) || poseLow(2) > poseUpp(2)) {
            std::cerr << __FILE__ << "," << __LINE__ << ": invalid bounds\n"
                    << "  x [" << poseLow(0) << "," << poseUpp(0) << "]\n"
                    << "  y [" << poseLow(1) << "," << poseUpp(1) << "]\n"
                    << "  theta [" << (180.0 / M_PI * poseLow(2)) << "," << (180.0 / M_PI * poseUpp(2)) << "] deg\n"
                    << std::endl;
            return false;
        }

        // Distance min and max is given by the distance of arc and four vertices:
        // - the arc is between Rot(poseLow(2))*pSrc and Rot(poseUpp(2))*pSrc
        // - the vertices are pDst -/+ [poseLow(0) : poseUpp(0), poseLow(1) : poseUpp(1)]
        radius = pSrc.norm();
        beta = atan2(pSrc.y(), pSrc.x());
        arcBeg = beta + poseLow(2);
        arcEnd = beta + poseUpp(2);
        pArcBeg << radius * cos(arcBeg), radius * sin(arcBeg);
        pArcEnd << radius * cos(arcEnd), radius * sin(arcEnd);
        xl = pDst(0) - poseUpp(0);
        xu = pDst(0) - poseLow(0);
        yl = pDst(1) - poseUpp(1);
        yu = pDst(1) - poseLow(1);
        v1 << xl, yl;
        v2 << xu, yl;
        v3 << xu, yu;
        v4 << xl, yu;

        //        GLORES_PRINT_VARIABLE(xl);
        //        GLORES_PRINT_VARIABLE(xu);
        //        GLORES_PRINT_VARIABLE(yl);
        //        GLORES_PRINT_VARIABLE(yu);
        //        GLORES_PRINT_VARIABLE(pArcBeg.transpose());
        //        GLORES_PRINT_VARIABLE(pArcEnd.transpose());

        // Computing dmax
        dist = computeMaxDistancePointArc(v1, radius, pArcBeg, pArcEnd);
        if (dist > dmax) {
            dmax = dist;
        }
        dist = computeMaxDistancePointArc(v2, radius, pArcBeg, pArcEnd);
        if (dist > dmax) {

            dmax = dist;
        }
        dist = computeMaxDistancePointArc(v3, radius, pArcBeg, pArcEnd);
        if (dist > dmax) {
            dmax = dist;
        }
        dist = computeMaxDistancePointArc(v4, radius, pArcBeg, pArcEnd);
        if (dist > dmax) {
            dmax = dist;
        }

        // Computing dmin:
        // Case 1: one arc point is inside the rectangle v1,v2,v3,v4
        // Case 2: test segment arc intersection
        dmin = 1e+6;
        if ((xl <= pArcBeg(0) && pArcBeg(0) <= xu && yl <= pArcBeg(1) && pArcBeg(1) <= yu) ||
                (xl <= pArcEnd(0) && pArcEnd(0) <= xu && yl <= pArcEnd(1) && pArcEnd(1) <= yu) ||
                intersectSegmentArcHor(xl, xu, yl, radius, pArcBeg, pArcEnd) ||
                intersectSegmentArcHor(xl, xu, yu, radius, pArcBeg, pArcEnd) ||
                intersectSegmentArcVer(xl, yl, yu, radius, pArcBeg, pArcEnd) ||
                intersectSegmentArcVer(xu, yl, yu, radius, pArcBeg, pArcEnd)) {
            dmin = 0.0;
        } else {
            // Checks all the distances between arc and rectangle vertices
            dist = computeDistancePointArc(v1, radius, pArcBeg, pArcEnd);
            //            GLORES_PRINT_VARIABLE(computeDistancePointArc(v1, radius, arcBeg, arcEnd));
            if (dist < dmin) {
                dmin = dist;
            }
            dist = computeDistancePointArc(v2, radius, pArcBeg, pArcEnd);
            //            GLORES_PRINT_VARIABLE(computeDistancePointArc(v2, radius, arcBeg, arcEnd));
            if (dist < dmin) {
                dmin = dist;
            }
            dist = computeDistancePointArc(v3, radius, pArcBeg, pArcEnd);
            //            GLORES_PRINT_VARIABLE(computeDistancePointArc(v3, radius, arcBeg, arcEnd));
            if (dist < dmin) {
                dmin = dist;
            }
            dist = computeDistancePointArc(v4, radius, pArcBeg, pArcEnd);
            //            GLORES_PRINT_VARIABLE(computeDistancePointArc(v4, radius, arcBeg, arcEnd));
            if (dist < dmin) {
                dmin = dist;
            }
            // Checks all the distances between each edge (segments) of rectangle and pArcBeg
            //            GLORES_PRINT_VARIABLE(computeDistancePointSegmentHor(pArcBeg, xl, xu, yl));
            dist = computeDistancePointSegmentHor(pArcBeg, xl, xu, yl);
            if (dist < dmin) {
                dmin = dist;
            }
            //            GLORES_PRINT_VARIABLE(computeDistancePointSegmentHor(pArcBeg, xl, xu, yu));
            dist = computeDistancePointSegmentHor(pArcBeg, xl, xu, yu);
            if (dist < dmin) {
                dmin = dist;
            }
            //            GLORES_PRINT_VARIABLE(computeDistancePointSegmentVer(pArcBeg, xl, yl, yu));
            dist = computeDistancePointSegmentVer(pArcBeg, xl, yl, yu);
            if (dist < dmin) {
                dmin = dist;
            }
            //            GLORES_PRINT_VARIABLE(computeDistancePointSegmentVer(pArcBeg, xu, yl, yu));
            dist = computeDistancePointSegmentVer(pArcBeg, xu, yl, yu);
            if (dist < dmin) {
                dmin = dist;
            }
            // Checks all the distances between each edge (segments) of rectangle and pArcEnd
            //            GLORES_PRINT_VARIABLE(computeDistancePointSegmentHor(pArcEnd, xl, xu, yl));
            dist = computeDistancePointSegmentHor(pArcEnd, xl, xu, yl);
            if (dist < dmin) {
                dmin = dist;
            }
            //            GLORES_PRINT_VARIABLE(computeDistancePointSegmentHor(pArcEnd, xl, xu, yu));
            dist = computeDistancePointSegmentHor(pArcEnd, xl, xu, yu);
            if (dist < dmin) {
                dmin = dist;
            }
            //            GLORES_PRINT_VARIABLE(computeDistancePointSegmentVer(pArcEnd, xl, yl, yu));
            dist = computeDistancePointSegmentVer(pArcEnd, xl, yl, yu);
            if (dist < dmin) {
                dmin = dist;
            }
            //            GLORES_PRINT_VARIABLE(computeDistancePointSegmentVer(pArcEnd, xu, yl, yu));
            dist = computeDistancePointSegmentVer(pArcEnd, xu, yl, yu);
            if (dist < dmin) {
                dmin = dist;
            }
            // Checks all the distances between
            //if (insideAngleInterval(0.0, arcBeg, arcEnd) || ) && (yl <= 0.0 && 0 <= yu)) {
            if (insideAngleInterval(1.0, 0.0, pArcBeg, pArcEnd) && (yl <= 0.0 && 0 <= yu)) {
                dist = fabs(radius - xl);
                if (dist < dmin) {
                    dmin = dist;
                }
                dist = fabs(radius - xu);
                if (dist < dmin) {
                    dmin = dist;
                }
            }
            //if (insideAngleInterval(M_PI, arcBeg, arcEnd) && (yl <= 0.0 && 0 <= yu)) {
            if (insideAngleInterval(-1.0, 0.0, pArcBeg, pArcEnd) && (yl <= 0.0 && 0 <= yu)) {
                dist = fabs(-radius - xl);
                if (dist < dmin) {
                    dmin = dist;
                }
                dist = fabs(-radius - xu);
                if (dist < dmin) {
                    dmin = dist;
                }
            }
            // Checks all the distances between

            if (insideAngleInterval(0.0, 1.0, pArcBeg, pArcEnd) && (xl <= 0.0 && 0 <= xu)) {
                dist = fabs(radius - yl);
                if (dist < dmin) {
                    dmin = dist;
                }
                dist = fabs(radius - yu);
                if (dist < dmin) {
                    dmin = dist;
                }
            }
            if (insideAngleInterval(0.0, -1.0, pArcBeg, pArcEnd) && (xl <= 0.0 && 0 <= xu)) {
                dist = fabs(-radius - yl);
                if (dist < dmin) {
                    dmin = dist;
                }
                dist = fabs(-radius - yu);
                if (dist < dmin) {
                    dmin = dist;
                }
            }
        }
        return true;
    }

    //    bool computeDistanceMinMax3(const Point2d& pSrc, const Point2d& pDst, const Eigen::Vector3d& poseLow, const Eigen::Vector3d& poseUpp, double& dmin, double& dmax, Point2d& pMin, Point2d& pMax, VectorPoint2d& keypoints) {
    //        double radius, xl, xu, yl, yu;
    //        Point2d v1, v2, v3, v4;
    //        double dist1, dist2, dist3, dist4;
    //        
    //        radius = pSrc.norm();
    //        xl = pDst(0) - poseUpp(0);
    //        xu = pDst(0) - poseLow(0);
    //        yl = pDst(1) - poseUpp(1);
    //        yu = pDst(1) - poseLow(1);
    //        v1 << xl, yl;
    //        v2 << xu, yl;
    //        v3 << xu, yu;
    //        v4 << xl, yu;
    //        dist1 = v1.norm();
    //        dist2 = v2.norm();
    //        dist3 = v3.norm();
    //        dist4 = v4.norm();
    //        
    //        return true;
    //    }

    double computeDistancePointArc(const Point2d& p, double radius, const Point2d pBeg, const Point2d pEnd) {
        Point2d nearest;
        return computeDistancePointArc(p, radius, pBeg, pEnd, nearest);
    }

    double computeDistancePointArc(const Point2d& p, double radius, double angleBeg, double angleEnd) {
        Point2d nearest;
        return computeDistancePointArc(p, radius, angleBeg, angleEnd, nearest);
    }

    double computeDistancePointArc(const Point2d& p, double radius, const Point2d pBeg, const Point2d pEnd, Point2d& nearest) {



        double scalarE = -pEnd(0) * p(1) + pEnd(1) * p(0);
        double scalarB = pBeg(0) * p(1) - pBeg(1) * p(0);
        int reverseCone = 0;
        if (pBeg(0) * pEnd(1) - pEnd(0) * pBeg(1) < 0) {
            reverseCone = 1;
            scalarE = -scalarE;
            scalarB = -scalarB;
        }

        if (!reverseCone) {
            if ((scalarE > 0)&&(scalarB > 0)) {
                nearest = (p / p.norm()) * radius;
                return fabs(p.norm() - radius);
            } else {
                double mB = (p - pBeg).norm();
                double mE = (p - pEnd).norm();
                nearest = ((mB < mE) ? pBeg : pEnd);
                return ((mB < mE) ? mB : mE);
            }
        }
        else {//reverse cone
            if ((scalarE > 0)&&(scalarB > 0)) {
                double mB = (p - pBeg).norm();
                double mE = (p - pEnd).norm();
                nearest = ((mB < mE) ? pBeg : pEnd);
                return ((mB < mE) ? mB : mE);
            } else {
                nearest = (p / p.norm()) * radius;
                return fabs(p.norm() - radius);
            }
        }




        assert(0);
    }

    double computeDistancePointArc(const Point2d& p, double radius, double angleBeg, double angleEnd, Point2d& nearest) {
        Point2d pBeg, pEnd;


        pBeg << radius * cos(angleBeg), radius * sin(angleBeg);
        pEnd << radius * cos(angleEnd), radius * sin(angleEnd);
        double scalarE = -pEnd(0) * p(1) + pEnd(1) * p(0);
        double scalarB = pBeg(0) * p(1) - pBeg(1) * p(0);
        int reverseCone = 0;

        if ((angleEnd - angleBeg) > M_PI) {
            reverseCone = 1;
            scalarE = -scalarE;
            scalarB = -scalarB;
        }

        if (!reverseCone) {
            if ((scalarE > 0)&&(scalarB > 0)) {
                nearest = (p / p.norm()) * radius;
                return fabs(p.norm() - radius);
            }//else if ((scalarE>0)&&(scalarB<0)) {
                //nearest=pBeg;
                //return (p - pBeg).norm();
                //} else if ((scalarE<0)&&(scalarB>0)) {
                //    nearest=pEnd;

                //    return (p - pEnd).norm();
                //}
            else {
                double mB = (p - pBeg).norm();
                double mE = (p - pEnd).norm();
                nearest = ((mB < mE) ? pBeg : pEnd);
                return ((mB < mE) ? mB : mE);
            }
        }
        else {//reverse cone
            if ((scalarE > 0)&&(scalarB > 0)) {
                double mB = (p - pBeg).norm();
                double mE = (p - pEnd).norm();
                nearest = ((mB < mE) ? pBeg : pEnd);
                return ((mB < mE) ? mB : mE);
            } else {
                nearest = (p / p.norm()) * radius;
                return fabs(p.norm() - radius);
            }
        }




        assert(0);
    }

    double computeDistancePointArcOld(const Point2d& p, double radius, double angleBeg, double angleEnd, Point2d& nearest) {
        Point2d pBeg, pEnd;
        double rho, alpha, distBeg, distEnd;

        rho = p.norm();
        alpha = fmod(atan2(p.y(), p.x()) + 2.0 * M_PI, 2.0 * M_PI);
        if (insideAngleInterval(alpha, angleBeg, angleEnd)) {
            nearest << radius * cos(alpha), radius * sin(alpha);
            return fabs(rho - radius);
        } else {
            pBeg << radius * cos(angleBeg), radius * sin(angleBeg);
            pEnd << radius * cos(angleEnd), radius * sin(angleEnd);
            distBeg = (p - pBeg).norm();
            distEnd = (p - pEnd).norm();
            if (distBeg < distEnd) {
                nearest = pBeg;
                return distBeg;
            } else {
                nearest = pEnd;
                return distEnd;
            }
        }
        assert(0); // it must return a value before!
        return 0.0;
    }

    double computeMaxDistancePointArc(const Point2d& p, double radius, double angleBeg, double angleEnd) {
        Point2d farthest;
        return computeMaxDistancePointArc(p, radius, angleBeg, angleEnd, farthest);
    }

    double computeMaxDistancePointArc(const Point2d& p, double radius, double angleBeg, double angleEnd, Point2d& farthest) {
        //Point2d farold;
        //computeDistancePointArcOld(-p, radius, angleBeg, angleEnd, farold);

        computeDistancePointArc(-p, radius, angleBeg, angleEnd, farthest);

        //if ((farthest-farold).norm()>1e-5) {
        //    computeDistancePointArcOld(-p, radius, angleBeg, angleEnd, farold);

        //    computeDistancePointArc(-p, radius, angleBeg, angleEnd, farthest);


        //    double a=1;
        //}
        return (farthest - p).norm();
    }

    double computeMaxDistancePointArc(const Point2d& p, double radius, const Point2d pBeg, const Point2d pEnd) {
        Point2d farthest;
        return computeMaxDistancePointArc(p, radius, pBeg, pEnd, farthest);
    }

    double computeMaxDistancePointArc(const Point2d& p, double radius, const Point2d pBeg, const Point2d pEnd, Point2d& farthest) {
        //Point2d farold;
        //computeDistancePointArcOld(-p, radius, angleBeg, angleEnd, farold);

        computeDistancePointArc(-p, radius, pBeg, pEnd, farthest);

        //if ((farthest-farold).norm()>1e-5) {
        //    computeDistancePointArcOld(-p, radius, angleBeg, angleEnd, farold);

        //    computeDistancePointArc(-p, radius, angleBeg, angleEnd, farthest);


        //    double a=1;
        //}
        return (farthest - p).norm();
    }

    double computeDistancePointSegment(const Point2d& p, const Point2d& seg1, const Point2d& seg2, Point2d& nearest) {
        Point2d dir;
        double segLen, theta, rho, t;

        computeLine(seg1, seg2, theta, rho, dir);
        segLen = dir.norm();
        dir.normalize();
        t = (p - seg1).dot(dir);

        //        GLORES_PRINT_VARIABLE(p.transpose());
        //        GLORES_PRINT_VARIABLE(seg1.transpose());
        //        GLORES_PRINT_VARIABLE(seg2.transpose());
        //        GLORES_PRINT_VARIABLE(t);
        //        GLORES_PRINT_VARIABLE(segLen);

        if (0.0 < t && t < segLen) {
            projectPointOnLine(p, theta, rho, nearest);
            return (p - nearest).norm();
        } else if (t <= 0.0) {
            nearest = seg1;
            return (p - seg1).norm();
        } else {
            nearest = seg2;
            return (p - seg2).norm();
        }
    }

    double computeDistancePointSegmentHor(const Point2d& p, double xl, double xu, double y) {
        if (xl <= p(0) && p(0) <= xu) {
            return fabs(p(1) - y);
        } else if (p(0) < xl) {
            return sqrt((p(0) - xl) * (p(0) - xl) + (p(1) - y) * (p(1) - y));
        } else {
            return sqrt((p(0) - xu) * (p(0) - xu) + (p(1) - y) * (p(1) - y));
        }
    }

    double computeDistancePointSegmentVer(const Point2d& p, double x, double yl, double yu) {
        if (yl <= p(1) && p(1) <= yu) {
            return fabs(p(0) - x);
        } else if (p(1) < yl) {
            return sqrt((p(0) - x) * (p(0) - x) + (p(1) - yl) * (p(1) - yl));
        } else {
            return sqrt((p(0) - x) * (p(0) - x) + (p(1) - yu) * (p(1) - yu));
        }
    }

    //    double computeDistanceSegmentArc(const Point2d& p1, const Point2d& p2, double radius, double angleBeg, double angleEnd) {
    //        double rho1, alpha1, rho2, alpha2;
    //
    //        rho1 = p1.norm();
    //        alpha1 = normalizeAngle2Pi(atan2(p1.y(), p1.x()));
    //        rho2 = p1.norm();
    //        alpha2 = normalizeAngle2Pi(atan2(p1.y(), p1.x()));
    //
    //    }

    bool intersectCircleHalfPlane(double theta, double rho, double radius, double& angleBeg, double& angleEnd, Point2d& nearestLine, Point2d& nearestCircle) {
        double thetaNorm, rhoNorm, delta;

        // Check if radius is positive
        if (radius <= 0.0) {
            std::cerr << __FILE__ << "," << __LINE__ << ": negative radius " << radius << std::endl;
            return false;
        }

        // Two sets of parameters to represent the same line: 
        // (theta, rho) and (theta + M_PI, -rho).
        // We use the two sets to distinguish among the half-planes divided 
        // by the line. 
        // The normalized parameters are s.t. rho > 0.0. 
        if (rho < 0.0) {
            thetaNorm = normalizeAngle2Pi(theta + M_PI);
            rhoNorm = -rho;
        } else {
            thetaNorm = theta;
            rhoNorm = rho;
        }
        //GLORES_PRINT_MSG(" hplane original: x * " << cos(theta) << " + y * " << sin(theta) << " >= " << rho);
        //GLORES_PRINT_MSG(" hplane normalized: x * " << cos(thetaNorm) << " + y * " << sin(thetaNorm) << " - (" << rhoNorm << ") >= 0.0");
        // If the rho is greater than radius, no intersection occurs
        //        GLORES_PRINT_VARIABLE(radius);
        //        GLORES_PRINT_VARIABLE(rho);
        //        GLORES_PRINT_VARIABLE(rhoNorm);
        //        GLORES_PRINT_VARIABLE(180.0 / M_PI * theta);
        //        GLORES_PRINT_VARIABLE(180.0 / M_PI * thetaNorm);
        if (rhoNorm > radius) {
            if (rho > 0.0) {
                angleBeg = thetaNorm;
                angleEnd = thetaNorm;
                nearestLine(0) = rho * cos(theta);
                nearestLine(1) = rho * sin(theta);
                nearestCircle(0) = radius * cos(theta);
                nearestCircle(1) = radius * sin(theta);
                return false;
            } else {
                angleBeg = 0.0;
                angleEnd = 2.0 * M_PI - 1e-6;
                return true;
            }
        }

        // Case with intersection: computed the interval angle
        delta = acos(rho / radius);
        angleBeg = theta - delta;
        angleEnd = theta + delta;

        //        GLORES_PRINT_VARIABLE(180.0 / M_PI * delta);
        //        GLORES_PRINT_VARIABLE(180.0 / M_PI * angleBeg);
        //        GLORES_PRINT_VARIABLE(180.0 / M_PI * angleEnd);
        return true;
    }

    bool intersectCircularArcHalfPlane(double theta, double rho, double radius, double arcBeg, double arcEnd, Point2d& nearestLine, Point2d& nearestArc) {
        std::vector<AngleInterval> arcIntervals;
        return intersectCircularArcHalfPlane(theta, rho, radius, arcBeg, arcEnd, nearestLine, nearestArc, arcIntervals);
    }

    bool intersectCircularArcHalfPlane(double theta, double rho, double radius, double arcBeg, double arcEnd, Point2d& nearestLine, Point2d& nearestArc, std::vector<AngleInterval>& arcIntervals) {
        double intersBeg, intersEnd, distBeg, distEnd;
        double thetaNorm, rhoNorm;
        bool intersectionOn, arcBegIn, arcEndIn;
        Point2d ptest;

        // Finds the intersection of the half-plane and the supporting circle.
        // If there is an intersection in the form of an intersection arc 
        // (between angles [intersBeg, intersEnd]), it finds overlap between 
        // the intersection arc and the given arc. 
        intersectionOn = intersectCircleHalfPlane(theta, rho, radius, intersBeg, intersEnd, nearestLine, nearestArc);
        //arcBegIn = insideAngleInterval(arcBeg, intersBeg, intersEnd);
        //arcEndIn = insideAngleInterval(arcEnd, intersBeg, intersEnd);
        int intersNum = intersectAngleInterval(std::make_pair(intersBeg, intersEnd), std::make_pair(arcBeg, arcEnd), arcIntervals);

        //        GLORES_PRINT_VARIABLE(intersectionOn);
        //        GLORES_PRINT_VARIABLE(intersNum);
        //        GLORES_PRINT_VARIABLE(180.0 / M_PI * intersBeg);
        //        GLORES_PRINT_VARIABLE(180.0 / M_PI * intersEnd);
        //        GLORES_PRINT_VARIABLE(180.0 / M_PI * arcBeg);
        //        GLORES_PRINT_VARIABLE(180.0 / M_PI * arcEnd);

        if (intersectionOn && intersNum > 0) {
            return true;
        } else {
            // Two sets of parameters to represent the same line: 
            // (theta, rho) and (theta + M_PI, -rho).
            // We use the two sets to distinguish among the half-planes divided 
            // by the line. 
            // The normalized parameters are s.t. rho > 0.0. 
            if (rho < 0.0) {
                thetaNorm = normalizeAngle2Pi(theta + M_PI);
                rhoNorm = -rho;
            } else {
                thetaNorm = theta;
                rhoNorm = rho;
            }
            // With no intersection, it finds the point on the line which is nearest to the circular arc. 
            // If the normal point, arc angular span cover []
            //            GLORES_PRINT_VARIABLE(180.0 / M_PI * arcBeg);
            //            GLORES_PRINT_VARIABLE(180.0 / M_PI * arcEnd);
            //            GLORES_PRINT_VARIABLE(180.0 / M_PI * theta);
            //            GLORES_PRINT_VARIABLE(insideAngleInterval(theta, arcBeg, arcEnd));
            if (insideAngleInterval(theta, arcBeg, arcEnd)) {
                nearestLine(0) = rho * cos(theta);
                nearestLine(1) = rho * sin(theta);
                nearestArc(0) = radius * cos(theta);
                nearestArc(1) = radius * sin(theta);
            } else {
                distBeg = fabs(radius * cos(arcBeg - theta) - rho);
                distEnd = fabs(radius * cos(arcEnd - theta) - rho);
                if (distBeg < distEnd) {
                    nearestArc(0) = radius * cos(arcBeg);
                    nearestArc(1) = radius * sin(arcBeg);
                    projectPointOnLine(nearestArc, theta, rho, nearestLine);
                } else {
                    nearestArc(0) = radius * cos(arcEnd);
                    nearestArc(1) = radius * sin(arcEnd);
                    projectPointOnLine(nearestArc, theta, rho, nearestLine);
                }
            }
        }
        return false;
    }

    bool intersectCircularArcHalfPlaneSegment(const Point2d& seg1, const Point2d& seg2, double radius, double arcBeg, double arcEnd, Point2d& nearestLine, Point2d& nearestArc) {
        std::vector<AngleInterval> arcIntervals;
        std::vector<DistancePointPair> candidates;
        return intersectCircularArcHalfPlaneSegment(seg1, seg2, radius, arcBeg, arcEnd, nearestLine, nearestArc, arcIntervals, candidates);
    }

    bool intersectCircularArcHalfPlaneSegment(const Point2d& seg1, const Point2d& seg2, double radius, double arcBeg, double arcEnd, Point2d& nearestLine, Point2d& nearestArc, std::vector<AngleInterval>& arcIntervals, std::vector<DistancePointPair>& candidates) {
        double intersBeg, intersEnd;
        double theta, rho, thetaNorm, rhoNorm;
        bool intersectionOn;
        double dist, segLen;
        Point2d pArc, pSeg, dir, t;

        computeLine(seg1, seg2, theta, rho, dir);
        //GLORES_PRINT_VARIABLE(180.0 / M_PI * theta);

        //        GLORES_PRINT_VARIABLE(180.0 / M_PI * arcBeg);
        //        GLORES_PRINT_VARIABLE(180.0 / M_PI * arcEnd);
        //        GLORES_PRINT_VARIABLE(180.0 / M_PI * theta);
        //        GLORES_PRINT_VARIABLE(rho);

        // Finds the intersection of the half-plane and the supporting circle.
        // If there is an intersection in the form of an intersection arc 
        // (between angles [intersBeg, intersEnd]), it finds overlap between 
        // the intersection arc and the given arc. 
        intersectionOn = intersectCircleHalfPlane(theta, rho, radius, intersBeg, intersEnd, nearestLine, nearestArc);
        int intersNum = intersectAngleInterval(std::make_pair(intersBeg, intersEnd), std::make_pair(arcBeg, arcEnd), arcIntervals);

        if (intersectionOn && intersNum > 0) {
            return true;
        }

        // Two sets of parameters to represent the same line: 
        // (theta, rho) and (theta + M_PI, -rho).
        // We use the two sets to distinguish among the half-planes divided 
        // by the line. 
        // The normalized parameters are s.t. rho > 0.0. 
        if (rho < 0.0) {
            thetaNorm = normalizeAngle2Pi(theta + M_PI);
            rhoNorm = -rho;
        } else {
            thetaNorm = theta;
            rhoNorm = rho;
        }
        // With no intersection, it finds the point on the line which is nearest to the circular arc. 
        // If the normal point, arc angular span cover []      
        //        GLORES_PRINT_VARIABLE(180.0 / M_PI * theta);
        //        GLORES_PRINT_VARIABLE(180.0 / M_PI * thetaNorm);
        //        GLORES_PRINT_VARIABLE(180.0 / M_PI * normalizeAngle2Pi(arcBeg));
        //        GLORES_PRINT_VARIABLE(180.0 / M_PI * normalizeAngle2Pi(arcEnd));
        //        GLORES_PRINT_VARIABLE(insideAngleInterval(theta, arcBeg, arcEnd));
        //        GLORES_PRINT_VARIABLE(insideAngleInterval(thetaNorm, arcBeg, arcEnd));
        if (insideAngleInterval(theta, arcBeg, arcEnd)) {
            pArc(0) = radius * cos(theta);
            pArc(1) = radius * sin(theta);
            dist = computeDistancePointSegment(pArc, seg1, seg2, pSeg);
            //            dist = fabs(rhoNorm - radius);
            //            pSeg(0) = rho * cos(theta);
            //            pSeg(1) = rho * sin(theta);            
            candidates.push_back(std::make_tuple(dist, pSeg, pArc));
            //            GLORES_PRINT_MSG("internal: dist " << dist << " pSeg [" << pSeg.transpose() << "], pArc [" << pArc.transpose() << "]");
        }

        // Distance between arc endpoint arcBeg and segment (seg1, seg2)
        pArc(0) = radius * cos(arcBeg);
        pArc(1) = radius * sin(arcBeg);
        dist = computeDistancePointSegment(pArc, seg1, seg2, pSeg);
        candidates.push_back(std::make_tuple(dist, pSeg, pArc));
        //        GLORES_PRINT_MSG("arcBeg-segment: dist " << dist << " pSeg [" << pSeg.transpose() << "], pArc [" << pArc.transpose() << "]");

        // Distance between arc endpoint arcEnd and segment (seg1, seg2)
        pArc(0) = radius * cos(arcEnd);
        pArc(1) = radius * sin(arcEnd);
        dist = computeDistancePointSegment(pArc, seg1, seg2, pSeg);
        candidates.push_back(std::make_tuple(dist, pSeg, pArc));
        //        GLORES_PRINT_MSG("arcEnd-segment: dist " << dist << " pSeg [" << pSeg.transpose() << "], pArc [" << pArc.transpose() << "]");

        // Distance segment endpoint seg1 and arc
        dist = computeDistancePointArc(seg1, radius, arcBeg, arcEnd, pArc);
        candidates.push_back(std::make_tuple(dist, seg1, pArc));
        //        GLORES_PRINT_MSG("seg1-arc: dist " << dist << " pSeg [" << seg1.transpose() << "], pArc [" << pArc.transpose() << "]");

        // Distance segment endpoint seg2 and arc
        dist = computeDistancePointArc(seg2, radius, arcBeg, arcEnd, pArc);
        candidates.push_back(std::make_tuple(dist, seg2, pArc));
        //        GLORES_PRINT_MSG("seg2-arc: dist " << dist << " pSeg [" << seg2.transpose() << "], pArc [" << pArc.transpose() << "]");

        dist = 1e+6;
        for (auto& cand : candidates) {
            if (std::get<0>(cand) < dist) {
                dist = std::get<0>(cand);
                nearestLine = std::get<1>(cand);
                nearestArc = std::get<2>(cand);
            }
        }
        //        GLORES_PRINT_VARIABLE(dist);
        //        GLORES_PRINT_VARIABLE(nearestLine.transpose());
        //        GLORES_PRINT_VARIABLE(nearestArc.transpose());
        return false;
    }

    //    bool intersectSegmentArc(const Point2d& seg1, const Point2d& seg2, double radius, double arcBeg, double arcEnd) {
    //        double dist1, dist2, t, alpha;
    //        Point2d pInters;
    //        dist1 = seg1.norm() - radius;
    //        dist2 = seg2.norm() - radius;
    //        if (dist1 * dist2 < 0.0) {
    ////            angSeg1 = atan2(seg1(1), seg1(0));
    ////            angSeg2 = atan2(seg2(1), seg2(0));
    //            dist1 = fabs(dist1);
    //            dist2 = fabs(dist2);
    //            t = dist1 / (dist1 - dist2);
    //            pInters = (1.0 - t) * seg1 + t * seg2;
    //            alpha = atan2(pInters(1), pInters(0));
    //            return insideAngleInterval(alpha, arcBeg, arcEnd);
    //        }
    //        return false;
    //    }

    bool intersectSegmentArcHor(double xl, double xu, double y, double radius, const Point2d pBeg, const Point2d pEnd) {
        double xInt1, xInt2;
        if (y * y > radius * radius) {
            return false;
        } else {
            xInt1 = sqrt(radius * radius - y * y);
            xInt2 = -sqrt(radius * radius - y * y);
            if (xl <= xInt1 && xInt1 <= xu && insideAngleInterval(xInt1, y, pBeg, pEnd)) {
                return true;
            } else if (xl <= xInt2 && xInt2 <= xu && insideAngleInterval(xInt2, y, pBeg, pEnd)) {
                return true;
            } else {
                return false;
            }
        }
    }

    bool intersectSegmentArcVer(double x, double yl, double yu, double radius, Point2d pBeg, Point2d pEnd) {
        double yInt1, yInt2, alpha1, alpha2;
        if (x * x > radius * radius) {
            return false;
        } else {
            yInt1 = sqrt(radius * radius - x * x);
            yInt2 = -sqrt(radius * radius - x * x);
            if (yl <= yInt1 && yInt1 <= yu && insideAngleInterval(x, yInt1, pBeg, pEnd)) {
                return true;
            } else if (yl <= yInt2 && yInt2 <= yu && insideAngleInterval(x, yInt2, pBeg, pEnd)) {
                return true;
            } else {
                return false;
            }
        }
    }

    bool intersectSegmentArcHor(double xl, double xu, double y, double radius, double arcBeg, double arcEnd) {
        double xInt1, xInt2, alpha1, alpha2;
        if (y * y > radius * radius) {
            return false;
        } else {
            xInt1 = sqrt(radius * radius - y * y);
            xInt2 = -sqrt(radius * radius - y * y);
            alpha1 = atan2(y, xInt1);
            alpha2 = atan2(y, xInt2);
            if (xl <= xInt1 && xInt1 <= xu && insideAngleInterval(alpha1, arcBeg, arcEnd)) {
                return true;
            } else if (xl <= xInt2 && xInt2 <= xu && insideAngleInterval(alpha2, arcBeg, arcEnd)) {
                return true;
            } else {
                return false;
            }
        }
    }

    bool intersectSegmentArcVer(double x, double yl, double yu, double radius, double arcBeg, double arcEnd) {
        double yInt1, yInt2, alpha1, alpha2;
        if (x * x > radius * radius) {
            return false;
        } else {
            yInt1 = sqrt(radius * radius - x * x);
            yInt2 = -sqrt(radius * radius - x * x);
            alpha1 = atan2(yInt1, x);
            alpha2 = atan2(yInt2, x);
            if (yl <= yInt1 && yInt1 <= yu && insideAngleInterval(alpha1, arcBeg, arcEnd)) {
                return true;
            } else if (yl <= yInt2 && yInt2 <= yu && insideAngleInterval(alpha2, arcBeg, arcEnd)) {
                return true;
            } else {
                return false;
            }
        }
    }

    void computeVoronoiDiagram(const VectorPoint2d& points, Voronoi2d& vd) {
        // Empty vector of segments
        //std::vector<Segment2d> segments;
        //boost::polygon::construct_voronoi(points.begin(), points.end(), segments.begin(), segments.end(), &vd);
        boost::polygon::construct_voronoi(points.begin(), points.end(), &vd);

        // See example https://www.boost.org/doc/libs/1_54_0/libs/polygon/example/voronoi_basic_tutorial.cpp

        //        // Checks whether all the voronoi cell contain a point and if the 
        //        // index of the point is the same in original input set
        //        for (auto& cell : vd.cells()) {
        //            if (cell.contains_point()) {
        //                std::size_t index = cell.source_index();
        //            }
        //        }
    }

    // --------------------------------------------------------------------
    // REGISTRATION
    // --------------------------------------------------------------------

    void associateNN(const VectorPoint2d& pointsSrc, const VectorPoint2d& pointsDst, const Transformation2d& dstTsrc, double distThres, std::vector<std::pair<int, int> >& associations) {
        double dist, distMin;
        int jmin;

        associations.clear();
        for (int i = 0; i < pointsSrc.size(); ++i) {
            distMin = 1e+6;
            jmin = -1;
            for (int j = 0; j < pointsDst.size(); ++j) {
                dist = (pointsDst[j] - dstTsrc * pointsSrc[i]).norm();
                if (jmin < 0 || dist < distMin) {
                    jmin = j;
                    distMin = dist;
                }
            }

            if (distMin < distThres && 0 <= jmin && jmin < pointsDst.size()) {
                associations.push_back(std::make_pair(i, jmin));
            }
        }
    }

    int computeTransform(const VectorPoint2d& pointsSrc, const VectorPoint2d& pointsDst, const std::vector<std::pair<int, int> >& associations, Transformation2d& dstTsrc) {
        Eigen::Vector2d tSrc = Eigen::Vector2d::Zero();
        Eigen::Vector2d tDst = Eigen::Vector2d::Zero();
        Eigen::Matrix2d S = Eigen::Matrix2d::Zero();
        int n = 0;

        // Computes mean values on each point set for valid associations
        for (auto& a : associations) {
            if (0 <= a.first && a.first < pointsSrc.size() && 0 <= a.second && a.second < pointsDst.size()) {
                tSrc += pointsSrc[a.first];
                tDst += pointsDst[a.second];
                n++;
            }
        }

        // Checks minimum number of valid points
        if (n < 2) {
            std::cerr << __FILE__ << "," << __LINE__ << ": not enough associations between the two point sets! only " << n << " < 2" << std::endl;
            return n;
        }

        // Computes transformation
        tSrc = (1.0 / n) * tSrc;
        tDst = (1.0 / n) * tDst;
        for (auto& a : associations) {
            if (0 <= a.first && a.first < pointsSrc.size() && 0 <= a.second && a.second < pointsDst.size()) {
                S += (pointsSrc[a.first] - tSrc) * (pointsDst[a.second] - tDst).transpose();
            }
        }
        double theta = atan2(S(0, 1) - S(1, 0), S(0, 0) + S(1, 1));
        Eigen::Rotation2Dd rot(theta);
        Eigen::Vector2d transl = tDst - (rot * tSrc);
        dstTsrc = Transformation2d::Identity();
        dstTsrc.prerotate(rot);
        dstTsrc.pretranslate(transl);

        return n;
    }

} // end of namespace

