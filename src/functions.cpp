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
#include <glores/functions.h>
#include <cassert> 
#include <algorithm>


namespace glores {

    // --------------------------------------------------------
    // COS-SIN FAST
    // --------------------------------------------------------

    void fastCosSin(double x, double& c, double& s) {
        constexpr double factor = 1.0 / (2.0 * M_PI);
        x *= factor;

        c = x - (0.25 + floor(x + 0.25));
        c *= 16.0 * (std::abs(c) - 0.5);
        c += 0.225 * c * (std::abs(c) - 1.0);

        s = x - floor(x + 0.5);
        s *= 16.0 * (0.5 - std::abs(s));
        s += 0.225 * s * (std::abs(s) - 1.0);
    }

    double fastAtan(double x) {
        const double a1 = 0.9998660;
        const double a3 = -0.3302995;
        const double a5 = 0.1801410;
        const double a7 = -0.0851330;
        const double a9 = 0.0208351;
        double x2 = x * x;
        double x4 = x2 * x2;
        return x * (a1 + x2 * (a3 + a7 * x4) + x4 * (a5 + a9 * x4));
    }

    double fastAtan2(double x, double y) {
        double ay = fabs(y);
        double ax = fabs(x);
        bool invert = ay > ax;
        double z = invert ? ax / ay : ay / ax; // z in range [0,1]
        double th = fastAtan(z); // th in range [0,M_PI/4] 
        if (invert) th = M_PI_2 - th; // th in range [0,M_PI/2]
        if (x < 0) th = M_PI - th; // th in range [0,M_PI]
        th = copysign(th, y); // th in range [-M_PI,M_PI]
    }

    double evaluateFourier(const std::vector<double>& coeffs, double theta) {
        double val, cth2, sth2, cth, sth, ctmp, stmp;
        int n;

        if (coeffs.size() % 2 != 0) {
            std::cerr << __FILE__ << "," << __LINE__ << ": the number of coefficients must be even: found " << coeffs.size() << std::endl;
        }
        n = (coeffs.size() / 2) - 1;

        cth2 = cos(theta);
        sth2 = sin(theta);
        cth = 1.0;
        sth = 0.0;
        val = 0.0;

        for (int k = 0; k <= n; ++k) {
            val += coeffs[2 * k] * cth + coeffs[2 * k + 1] * sth;
            ctmp = cth2 * cth - sth2 * sth;
            stmp = sth2 * cth + cth2 * sth;
            cth = ctmp;
            sth = stmp;
        }
        return val;
    }

    double rastrigin(double a, const Eigen::VectorXd& x) {
        int n = x.size();
        double val = a * n;

        for (int i = 0; i < n; ++i) {
            val += x(i) * x(i) - a * cos(2 * M_PI * x(i));
        }

        return val;
    }

    // --------------------------------------------------------
    // INTERVAL FUNCTIONS
    // --------------------------------------------------------

    void findLUSquare(double xmin, double xmax, double& ymin, double& ymax) {
        double xmin2, xmax2;

        if (xmin > xmax) std::swap(xmin, xmax);

        // Finds the minimum and maximum (currently disregarding if 0.0 lies in [xmin,xmax])
        xmin2 = xmin * xmin;
        xmax2 = xmax * xmax;
        if (xmin2 < xmax2) {
            ymin = xmin2;
            ymax = xmax2;
        } else {
            ymin = xmax2;
            ymax = xmin2;
        }

        // Adjust the minimum value when  if 0.0 lies in [xmin,xmax]
        if (xmin < 0.0 && 0.0 < xmax) {
            ymin = 0.0;
        }
    }

    void findLUCos(double a, double b, double& cmin, double& cmax) {
        double amod, bmod;

        if (a > b) std::swap(a, b);

        if (b - a >= 2.0 * M_PI) {
            cmin = -1.0;
            cmax = +1.0;
        } else {
            // Normalizes to circular interval [0, 2*M_PI[
            amod = fmod(a, 2.0 * M_PI);
            if (amod < 0.0) amod += 2.0 * M_PI;
            bmod = fmod(b, 2.0 * M_PI);
            if (bmod < 0.0) bmod += 2.0 * M_PI;
            // Case bmod < amod: for example [300,30[ deg: angle 0 is included.
            if (bmod < amod) {
                cmax = +1.0;
                if (bmod < M_PI && M_PI < amod) {
                    cmin = std::min(cos(amod), cos(bmod));
                } else {
                    cmin = -1.0;
                }
                //      if (M_PI < bmod || amod < M_PI) {
                //        cmin = -1.0;
                //      }
                //      else {
                //        cmin = std::min(cos(amod),cos(bmod));
                //      }
            } else {
                cmax = std::max(cos(amod), cos(bmod));
                if (amod < M_PI && M_PI < bmod) {
                    cmin = -1.0;
                } else {
                    cmin = std::min(cos(amod), cos(bmod));
                }
            }
        }
    }

    void findLUFourier(const std::vector<double>& coeffs, double theta0, double theta1, double& fourierMin, double& fourierMax) {
        double amplitude, phase, sinusoidMin, sinusoidMax;
        int n, i0, i1;

        if (coeffs.size() % 2 != 0) {
            std::cerr << __FILE__ << "," << __LINE__ << ": the number of coefficients must be even: found " << coeffs.size() << std::endl;
        }
        n = (coeffs.size() / 2) - 1;

        if (theta1 < theta0) {
            std::cerr << __FILE__ << "," << __LINE__ << ": invalid interval [" << theta0 << "," << theta1 << "]: swapping endpoints to continue" << std::endl;
            std::swap(theta0, theta1);
        }

        // fourierMin and fourierMax initialized with constant component
        fourierMin = coeffs[0];
        fourierMax = coeffs[0];
        for (int k = 1; k <= n; ++k) {
            // t_k = a_k * cos(2*k*theta) + b_k * sin(2*k*theta) = amplitude * cos(2*k*theta - phase)
            // Period of k-th terms is M_PI / k. 
            amplitude = sqrt(coeffs[2 * k] * coeffs[2 * k] + coeffs[2 * k + 1] * coeffs[2 * k + 1]);
            phase = atan2(coeffs[2 * k + 1], coeffs[2 * k]);
            //std::cout << "k " << k << ", amplitude " << amplitude << ", phase[deg] " << (180.0/M_PI*phase) << std::endl;
            // If the [theta0,theta1] interval is larger than period, then the whole sinusoid amplitude is considered.
            // Otherwise, a more refined evaluation is performed.
            findLUCos(2 * k * theta0 - phase, 2 * k * theta1 - phase, sinusoidMin, sinusoidMax);
            fourierMin += amplitude * sinusoidMin;
            fourierMax += amplitude * sinusoidMax;
        }
    }

    void findLURastrigin(double a, const Eigen::VectorXd& xmin, const Eigen::VectorXd& xmax, double& ymin, double& ymax) {
        double squareMin, squareMax;
        double cosMin, cosMax;
        int n = xmin.size();
        if (xmax.size() != n) {
            std::cerr << __FILE__ << "," << __LINE__ << ": invalid dimensions: xmin.size() " << xmin.size() << " != xmax.size() " << xmax.size() << std::endl;
            return;
        }

        ymin = a * n;
        ymax = a * n;
        for (int i = 0; i < n; ++i) {
            // Finds the interval of x(i)^2 function
            findLUSquare(xmin(i), xmax(i), squareMin, squareMax);
            // Finds the interval of -a * cos(2 * M_PI * x(i))
            findLUCos(2.0 * M_PI * xmin(i), 2.0 * M_PI * xmax(i), cosMin, cosMax);
            cosMin = a * cosMin;
            cosMax = a * cosMax;
            // Updates ymin and ymax
            ymin += squareMin - cosMax;
            ymax += squareMax - cosMin;
        }
    }

} // end of namespace

