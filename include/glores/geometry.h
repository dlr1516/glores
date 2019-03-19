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
#pragma once

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Geometry>

#include <vector>
#include <deque>

// In order to use glores custom points with the libraries boost::geometry and 
// boost::polygon
#include <boost/polygon/polygon.hpp>
#include <boost/polygon/voronoi.hpp>

#include <boost/geometry.hpp>
#include <boost/geometry/core/cs.hpp>
#include <boost/geometry/geometries/segment.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/register/point.hpp>
#include <boost/geometry/geometries/register/segment.hpp>



namespace glores {

    /**
     * Types for a 2D point in the form of Euclidean vector. 
     */
    typedef Eigen::Vector2d Point2d;
    typedef Eigen::Vector2f Point2f;
    typedef Eigen::Vector2i Point2i;

    /**
     * Types for a 2D matrix in the form of Euclidean vector. 
     */
    typedef Eigen::Matrix2d Matrix2d;
    typedef Eigen::Matrix2f Matrix2f;
    typedef Eigen::Matrix2i Matrix2i;

    /**
     * Types for STL vectors of different point types. 
     */
    typedef std::vector<Point2d, Eigen::aligned_allocator<Point2d> > VectorPoint2d;
    typedef std::vector<Point2f, Eigen::aligned_allocator<Point2f> > VectorPoint2f;
    typedef std::vector<Point2i, Eigen::aligned_allocator<Point2i> > VectorPoint2i;
    typedef std::deque<Point2d, Eigen::aligned_allocator<Point2d> > DequePoint2d;
    typedef std::deque<Point2f, Eigen::aligned_allocator<Point2f> > DequePoint2f;
    typedef std::deque<Point2i, Eigen::aligned_allocator<Point2i> > DequePoint2i;

    /**
     *  Type for reference frames
     */
    typedef Eigen::Isometry2d Transformation2d;
    typedef Eigen::Isometry2f Transformation2f;


    typedef std::pair<double, double> AngleInterval; 
    typedef std::tuple<double, Point2d, Point2d> DistancePointPair;

} // end of namespace glores

// ------------------------------------------------------------------
// BOOST::GEOMETRY LIBARY TRAITS AND DEFINITIONS
// ------------------------------------------------------------------

BOOST_GEOMETRY_REGISTER_POINT_2D(glores::Point2d, double, cs::cartesian, glores::Point2d::x(), glores::Point2d::y());
BOOST_GEOMETRY_REGISTER_POINT_2D(glores::Point2f, float, cs::cartesian, glores::Point2f::x(), glores::Point2f::y());
BOOST_GEOMETRY_REGISTER_POINT_2D(glores::Point2i, int, cs::cartesian, glores::Point2i::x(), glores::Point2i::y());

namespace glores {

    /**
     * Types for segments using boost::geometry library. 
     */
    typedef boost::geometry::model::segment<Point2d> Segment2d;
    typedef boost::geometry::model::segment<Point2f> Segment2f;
    typedef boost::geometry::model::segment<Point2i> Segment2i;

    /**
     * Types for polygon using boost::geometry library. 
     */
    typedef boost::geometry::model::polygon<Point2d> Polygon2d;
    typedef boost::geometry::model::polygon<Point2f> Polygon2f;
    typedef boost::geometry::model::polygon<Point2i> Polygon2i;
}

// ------------------------------------------------------------------
// BOOST::POLYGON LIBARY TRAITS AND DEFINITIONS
// ------------------------------------------------------------------

namespace boost {
    namespace polygon {

        // Traits for Point2d

        template <>
        struct geometry_concept<glores::Point2d> {
            typedef point_concept type;
        };

        template <>
        struct point_traits<glores::Point2d> {
            typedef double coordinate_type;

            static inline coordinate_type get(const glores::Point2d& point,
                    orientation_2d orient) {
                if (orient == HORIZONTAL)
                    return point(0);
                return point(1);
            }
        };

        // Traits for Segment2d

        template <>
        struct point_mutable_traits<glores::Point2d> {
            typedef int coordinate_type;

            static inline void set(glores::Point2d& point, orientation_2d orient, double value) {
                if (orient == HORIZONTAL)
                    point(0) = value;
                else
                    point(1) = value;
            }

            static inline glores::Point2d construct(double x_value, double y_value) {
                glores::Point2d retval;
                retval(0) = x_value;
                retval(1) = y_value;
                return retval;
            }
        };

        template <>
        struct geometry_concept<glores::Segment2d> {
            typedef segment_concept type;
        };

        template <>
        struct segment_traits<glores::Segment2d> {
            typedef double coordinate_type;
            typedef glores::Point2d point_type;

            static inline point_type get(const glores::Segment2d& segment, direction_1d dir) {
                return dir.to_int() ? segment.second : segment.first;
            }
        };
    }
}

namespace glores {

    typedef boost::polygon::voronoi_diagram<double> Voronoi2d;

}

namespace glores {

    // --------------------------------------------------------------
    // TYPE CONVERSION
    // --------------------------------------------------------------
    
    void poseToIsometry2d(double x, double y, double theta, glores::Transformation2d& transf);
    
    void poseToIsometry2f(float x, float y, float theta, glores::Transformation2f& transf);

    void poseToIsometry2d(const Eigen::Vector3d& pose, glores::Transformation2d& transf);

    void poseToIsometry2f(const Eigen::Vector3f& pose, glores::Transformation2f& transf);

    // --------------------------------------------------------------
    // FUNCTIONS
    // --------------------------------------------------------------

    /**
     * Returns the angle in the interval [0.0, 2*M_PI[
     * @param angle input angle in radians
     */
    double normalizeAngle2Pi(double angle);

    /**
     * Says if a given angle falls inside an angle interval with limits
     * [intervBeg, intervEnd] (where intervEnd is reached in CCW turn from 
     * intervEnd).
     * @param angle
     * @param intervBeg
     * @param intervEnd
     * @return 
     */
    bool insideAngleInterval(double angle, double intervBeg, double intervEnd);

    /**
     * Says if a given angle falls inside an angle interval interv. 
     * @param angle
     * @param interv
     * @return 
     */
    bool insideAngleInterval(double angle, const AngleInterval& interv);

    /**
     * Computes the intersection of two angular intervals in1 and in2, 
     * consisting of a vector with 0, 1 or 2 intersection intervals.
     * @param in1 interval to be intersected
     * @param in2 interval to be intersected
     * @param out vector of intersection intervals
     * @return the number of intervals in the intersection
     */
    int intersectAngleInterval(const AngleInterval& in1, const AngleInterval& in2, std::vector<AngleInterval>& out);

    void projectPointOnLine(const Point2d& p, double theta, double rho, Point2d& q);

    void computeLine(const Point2d& p1, const Point2d& p2, double& theta, double& rho);

    void computeLine(const Point2d& p1, const Point2d& p2, double& theta, double& rho, Point2d& d12);

    bool projectArcBound(double arcBeg, double arcEnd, VectorPoint2d& polygonBound, int vnum = 2);
    
    /**
     * Computes the projection of a given 2D point p into the Euclidean plane 
     * when the rigid transformation is givel in the interval [poseLow, poseUpp], 
     * i.e. translation x s.t. poseLow(0) <= x <= poseUpp(0),
     *      translation y s.t. poseLow(1) <= y <= poseUpp(1),
     *      rotation theta y s.t. poseLow(2) <= theta <= poseUpp(2). 
     * It returns a convex bounding set containing the projection. 
     * 
     * @param p the point to be projected
     * @param poseLow the lower bounds of [x, y, theta] interval
     * @param poseUpp the upper bounds of [x, y, theta] interval
     * @param convexBound the convex polygon that contains all the points
     * @param vnum the number of edges approximating the arc segment of projection
     * @return true if the input interval is correct
     */
    bool projectSE2Interval(const Point2d& p, const Eigen::Vector3d& poseLow, const Eigen::Vector3d& poseUpp, Polygon2d& convexBound, int vnum = 10);

    /**
     * Computes the min/max distance between points pSrc and pDst given the bounds 
     * on the transformation T=(R, t) where T lies in interval I = [poseLow, poseUpp]
     * (R is the rotation matrix, t the translation vector). 
     * More formally:
     * 
     *    dmin = min_{T \in I} \| R * pSrc + t - pDst  \|
     *    dmax = max_{T \in I} \| R * pSrc + t - pDst  \|
     * 
     * The distances are estimated by observing that:
     * 1) R * pSrc is an arc depending on the rotation angle interval [poseLow(2), poseUpp(2)];
     * 2) (pDst - t) is a rectangle depending on x and y interval of translation t.
     * The min/max distance is the distance between the angle and the rectangle
     * 
     * @param pSrc source point to be transformed
     * @param pDst destination point
     * @param poseLow the lower bounds of [x, y, theta] interval
     * @param poseUpp the upper bounds of [x, y, theta] interval
     * @param dmin minimum distance between pSrc and transformed pDst
     * @param dmax maximum distance between pSrc and transformed pDst
     * @return true if the input interval of transformation is correct
     */
    bool computeDistanceMinMax(const Point2d& pSrc, const Point2d& pDst, const Eigen::Vector3d& poseLow, const Eigen::Vector3d& poseUpp, double& dmin, double& dmax);

    /**
     * Computes the min/max distance between points pSrc and pDst given the bounds 
     * on the transformation T=(R, t) where T lies in interval I = [poseLow, poseUpp]
     * (R is the rotation matrix, t the translation vector). 
     * More formally:
     * 
     *    dmin = min_{T \in I} \| R * pSrc + t - pDst  \|
     *    dmax = max_{T \in I} \| R * pSrc + t - pDst  \|
     * 
     * The distances are estimated by observing that:
     * 1) R * pSrc is an arc depending on the rotation angle interval [poseLow(2), poseUpp(2)];
     * 2) (pDst - t) is a rectangle depending on x and y interval of translation t.
     * The min/max distance is the distance between the angle and the rectangle
     * 
     * @param pSrc source point to be transformed
     * @param pDst destination point
     * @param poseLow the lower bounds of [x, y, theta] interval
     * @param poseUpp the upper bounds of [x, y, theta] interval
     * @param dmin minimum distance between pSrc and transformed pDst
     * @param dmax maximum distance between pSrc and transformed pDst
     * @param pMin point of rectangle closest to the arc
     * @param pMax point of rectangle farthest to the arc
     * @return true if the input interval of transformation is correct
     */
    //bool computeDistanceMinMax(const Point2d& pSrc, const Point2d& pDst, const Eigen::Vector3d& poseLow, const Eigen::Vector3d& poseUpp, double& dmin, double& dmax, Point2d& pMin, Point2d& pMax);

    bool computeDistanceMinMax(const Point2d& pSrc, const Point2d& pDst, const Eigen::Vector3d& poseLow, const Eigen::Vector3d& poseUpp, double& dmin, double& dmax, Point2d& pMin, Point2d& pMax, VectorPoint2d& keypoints);

    bool computeDistanceMinMax2(const Point2d& pSrc, const Point2d& pDst, const Eigen::Vector3d& poseLow, const Eigen::Vector3d& poseUpp, double& dmin, double& dmax);

    bool computeDistanceMinMax2(const Point2d& pSrc, const Point2d& pDst, const Eigen::Vector3d& poseLow, const Eigen::Vector3d& poseUpp, double& dmin, double& dmax, Point2d& pMin, Point2d& pMax, VectorPoint2d& keypoints);
 
    
    bool computeDistanceMinMax3(const Point2d& pBeg, const Point2d& pEnd, const double beta, const double radius, const Point2d& pDst, const Eigen::Vector3d& poseLow, const Eigen::Vector3d& poseUpp, double& dmin, double& dmax);
    
    bool computeDistanceMinMax3(const Point2d& pSrc, const Point2d& pDst, const Eigen::Vector3d& poseLow, const Eigen::Vector3d& poseUpp, double& dmin, double& dmax);

    /**
     * Computes the distance of a point p from a circular arc defined by 
     * the radius and the angle interval [angleBeg, angleEnd] (angleEnd is reached
     * from angleBeg in counter-clockwise order). 
     * @param p the point 
     * @param radius the radius
     * @param angleBeg the begin angle of arc
     * @param angleEnd the end angle of arc
     * @return the distance 
     */
    double computeDistancePointArc(const Point2d& p, double radius, double angleBeg, double angleEnd);

    /**
     * Computes the distance of a point p from a circular arc defined by 
     * the radius and the angle interval [angleBeg, angleEnd] (angleEnd is reached
     * from angleBeg in counter-clockwise order). 
     * @param p the point 
     * @param radius the radius
     * @param angleBeg the begin angle of arc
     * @param angleEnd the end angle of arc
     * @param nearest the point of circular arc nearest to point p
     * @return the distance 
     */
    double computeDistancePointArc(const Point2d& p, double radius, double angleBeg, double angleEnd, Point2d& nearest);
    
    double computeDistancePointArc(const Point2d& p, double radius, const Point2d pBeg, const Point2d pEnd);
    
    /**
     * Computes the distance of a point p from a circular arc defined by
     * the radius and the angle interval [angleBeg, angleEnd] (angleEnd is reached
     * from angleBeg in counter-clockwise order).
     * @param p the point
     * @param radius the radius
     * @param angleBeg the begin angle of arc
     * @param angleEnd the end angle of arc
     * @param nearest the point of circular arc nearest to point p
     * @return the distance
     */
    double computeDistancePointArc(const Point2d& p, double radius, const Point2d pBeg, const Point2d pEnd, Point2d& nearest);
    
    
    
    double computeMaxDistancePointArc(const Point2d& p, double radius, double angleBeg, double angleEnd);
    
    double computeMaxDistancePointArc(const Point2d& p, double radius, double angleBeg, double angleEnd, Point2d& farthest);

    
    double computeMaxDistancePointArc(const Point2d& p, double radius, const Point2d pBeg, const Point2d pEnd);
    
    double computeMaxDistancePointArc(const Point2d& p, double radius, const Point2d pBeg, const Point2d pEnd, Point2d& farthest);
    
    
    
    
    double computeDistancePointSegment(const Point2d& p, const Point2d& seg1, const Point2d& seg2, Point2d& nearest);
    
    double computeDistancePointSegmentHor(const Point2d& p, double xl, double xu, double y);
    
    double computeDistancePointSegmentVer(const Point2d& p, double x, double yl, double yu);

    /**
     * Intersects the half-plane x*cos(theta) + y*sin(theta) - rho >= 0 and 
     * the circle with given radius and centred in the origin. 
     * @param theta the angle parameters of the line
     * @param rho the range parameters of the line
     * @param radius the radius of the circle
     * @param angleBeg the initial angle
     * @param angleEnd the final angle
     * @param nearest the point of line nearest to the circle, when no intersection occur 
     * @return true if the circle intersects the half-plane
     */
    bool intersectCircleHalfPlane(double theta, double rho, double radius, double& angleBeg, double& angleEnd, Point2d& nearestLine, Point2d& nearestCircle);

    bool intersectCircularArcHalfPlane(double theta, double rho, double radius, double arcBeg, double arcEnd, Point2d& nearestLine, Point2d& nearestArc);

    bool intersectCircularArcHalfPlane(double theta, double rho, double radius, double arcBeg, double arcEnd, Point2d& nearestLine, Point2d& nearestArc, std::vector<AngleInterval>& arcIntervals);

    bool intersectCircularArcHalfPlaneSegment(const Point2d& seg1, const Point2d& seg2, double radius, double arcBeg, double arcEnd, Point2d& nearestLine, Point2d& nearestArc);

    bool intersectCircularArcHalfPlaneSegment(const Point2d& seg1, const Point2d& seg2, double radius, double arcBeg, double arcEnd, Point2d& nearestLine, Point2d& nearestArc, std::vector<AngleInterval>& arcIntervals, std::vector<DistancePointPair>& candidates);
    
    bool intersectSegmentArc(const Point2d& seg1, const Point2d& seg2, double radius, double arcBeg, double arcEnd);
    
    bool intersectSegmentArcHor(double xl, double xu, double y, double radius, double arcBeg, double arcEnd);
    
    bool intersectSegmentArcVer(double x, double yl, double yu, double radius, double arcBeg, double arcEnd);

    
    bool intersectSegmentArcHor(double xl, double xu, double y, double radius, const Point2d pBeg, const Point2d pEnd);
    
    bool intersectSegmentArcVer(double x, double yl, double yu, double radius, const Point2d pBeg, const Point2d pEnd);

    
    
    /**
     * Computes Voronoi diagram over a set of points. It is just a wrapper of 
     * for the corresponding function in boost::polygon.
     * @param points the input point set
     * @param vd the Voronoi disgram
     */
    void computeVoronoiDiagram(const VectorPoint2d& points, Voronoi2d& vd);

    // --------------------------------------------------------------------
    // REGISTRATION
    // --------------------------------------------------------------------

    /**
     * 
     * @param features
     * @param landmarks
     * @param lTw
     * @param distThres
     * @param associations
     */
    void associateNN(const VectorPoint2d& pointsSrc, const VectorPoint2d& pointsDst, const Transformation2d& dstTsrc, double distThres, std::vector<std::pair<int, int> >& associations);

    /**
     * Computes the transformation from source point sets to destination point set. 
     * @param pointsSrc
     * @param pointsDst
     * @param associations
     * @param dstTsrc
     * @return 
     */
    int computeTransform(const VectorPoint2d& pointsSrc, const VectorPoint2d& pointsDst, const std::vector<std::pair<int, int> >& associations, Transformation2d& dstTsrc);

} // end of namespace glores

