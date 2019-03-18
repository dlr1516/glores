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
#pragma once

#include <iostream>
#include <vector>
#include <queue>
#include <cmath>
#include <Eigen/Dense>

namespace glores {

    const double PNEBI_ARG_MAX = 600.0;
    const double BIG_NUM = 1.0e+10;
    const double SMALL_NUM = 1.0e-10;

    // --------------------------------------------------------
    // COS-SIN FAST EVALUATION
    // --------------------------------------------------------

    /** Computes sine and cosine values using a parabolic approximation of sine.
     *    y = 16.0 * xn * (abs(xn) - 0.5)
     * where xn = x/(2*M_PI) is the normalized value of x over the period. 
     *
     * Code suggested here:
     *  http://forum.devmaster.net/t/fast-and-accurate-sine-cosine/9648/6
     *  http://stackoverflow.com/questions/18662261/fastest-implementation-of-sine-cosine-and-square-root-in-c-doesnt-need-to-b
     */
    void fastCosSin(double x, double& c, double& s);

    /**
     * Computes atan() using a polynomial approximation on interval [-1,1]. See:
     * 
     *  Abramowitz, Stegun, "Handbook of Mathematical Functions", 1965
     * 
     * @param x the argument that must be in interval [-1.0, 1.0]
     * @return the value of atan
     */
    double fastAtan(double x);

    /**
     * Computes atan2() using fastAtan(). It uses clever interval as suggested by:
     *    https://www.dsprelated.com/showarticle/1052.php
     * (one comment in particular improved the discussed code). 
     * @param x
     * @param y
     * @return the approximate atan2. 
     */
    double fastAtan2(double x, double y);

    /**
     * Computes the value of the given Fourier series at the given point theta. 
     * The user must provide the vector of serie coefficients:
     *   S(x) = \sum_{i=0}^{n} ( coeffs[2*i] * cos(i*x) + coeffs[2*i+1] * sin(i*x) )
     * @param coeffs the vector of coefficiens (vector size must be an even number!)
     * @param theta the point where to compute the Fourier serie value
     * @return the value of the function
     */
    double evaluateFourier(const std::vector<double>& coeffs, double theta);

    /**
     * Conputes the value of Rastrigin function. 
     * The dimension is given by the dimension of the vector argument.
     */
    double rastrigin(double a, const Eigen::VectorXd& x);

    // --------------------------------------------------------
    // INTERVAL FUNCTIONS
    // --------------------------------------------------------

    /**
     * Computes the lower and upper bounds, respectively ymin and ymax, of 
     * function y=x^2 on interval [xmin, xmax]
     * @param xmin the domain min
     * @param xmax the domain max
     * @param ymin the function min
     * @param ymax the function max
     */
    void findLUSquare(double xmin, double xmax, double& ymin,double& ymax);

    /** Computes lower and upper bounds of cosine function on a given interval.
     */
    void findLUCos(double a, double b, double& cmin, double& cmax);

    /** Computes lower and upper bounds of Fourier Series (represented by its coefficients)
     * on a given interval.
     * The vector of coefficients coeffs[i] are used in Fourier series:
     *   S(x) = \sum_{i=0}^{n} ( coeffs[2*i] * cos(2*i*x) + coeffs[2*i+1] * sin(2*i*x) )
     */
    void findLUFourier(const std::vector<double>& coeffs, double theta0, double theta1, double& fourierMin, double& fourierfMax);

    void findLURastrigin(double a, const Eigen::VectorXd& xmin, const Eigen::VectorXd& xmax, double& ymin, double& ymax);


} // end of namespace

