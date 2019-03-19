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
#include <Eigen/Dense>
//#include <random>
#include <glores/PointSetGlobalRegistration.h>
#include <glores/geometry.h>
//#include <glores/ParamMap.h>
//#include <glores/fileutils.h>
//#include <glores/thirdparty/gnuplot-iostream.h>
//#include <glores/Profiler.h>

#include "mex.h"

void convertMxArrayToVectorPoint2(const mxArray *Parray, glores::VectorPoint2d& points) {
    int nrows, ncols;
    double *P;
    glores::Point2d p;

    if (mxGetNumberOfDimensions(Parray) != 2) {
        mexErrMsgIdAndTxt("MATLAB:mxgetnzmax:inputNotDouble",
                "Input argument must be a double matrix.");
    }
    nrows = mxGetM(Parray);
    ncols = mxGetN(Parray);
    if (nrows != 2) {
        mexErrMsgIdAndTxt("MATLAB:mxgetnzmax:inputNotDouble",
                "Input argument must have 2 rows.");
    }

    P = mxGetPr(Parray);
    points.reserve(ncols);
    for (int c = 0; c < ncols; ++c) {
        p << P[c * nrows], P[c * nrows + 1];
        points.push_back(p);
    }
}

void convertMxArrayToVector3(const mxArray *Varray, Eigen::VectorXd& vec) {
    int nrows, ncols;
    double *V;
    glores::Point2d p;

    //  if (mxGetNumberOfDimensions(Varray) != 1) {
    //      mexErrMsgIdAndTxt("MATLAB:mxgetnzmax:inputNotDouble",
    //              "Input argument must be a double vector.");
    //  }
    nrows = mxGetM(Varray);
    ncols = mxGetN(Varray);
    if (nrows != 1 && ncols != 1) { //not sure if it is a column or row vector
        mexErrMsgIdAndTxt("MATLAB:mxgetnzmax:inputNotDouble",
                "Input argument must have 2 rows.");
    }
    if (nrows * ncols != 3) {
        mexErrMsgIdAndTxt("MATLAB:mxgetnzmax:inputNotDouble",
                "Input argument must have 3 elements.");
    }

    V = mxGetPr(Varray);
    vec.resize(3);
    vec << V[0], V[1], V[2];
}

void convertVectorToMxArray(Eigen::VectorXd& vec, mxArray *&result) {
    result = mxCreateDoubleMatrix((mwSize) vec.rows(), (mwSize) 1, mxREAL);
    mxDouble *rd = mxGetDoubles(result);
    for (int i = 0; i < vec.rows(); ++i) {
        rd[i] = vec(i);
    }
}

void convertDoubleToMxArray(double v, mxArray *&result) {
    result = mxCreateDoubleMatrix((mwSize) 1, (mwSize) 1, mxREAL);
    mxDouble *rd = mxGetDoubles(result);
    rd[0] = v;
}

void convertStatToMxArray(const std::vector<glores::PointSetGlobalRegistration::BoundStats>& stats, mxArray *&result, int fieldId) {
    result = mxCreateDoubleMatrix((mwSize) stats.size(), (mwSize) 1, mxREAL);
    mxDouble *rd = mxGetDoubles(result);
    for (int i = 0; i < stats.size(); ++i) {
        if (fieldId == 0) {
            rd[i] = stats[i].lowerBound;
        } else if (fieldId == 1) {
            rd[i] = stats[i].upperBound;
        } else if (fieldId == 2) {
            rd[i] = (stats[i].relaxationBoundOn ? 1.0 : 0.0);
        }
    }
}

void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[]) {
    /* If the mex file was built using interleaved complex flag, display
     * an error and exit.
     */

    /* Declare variable */

    const mxArray *P, *Q, *betamin, *betamax, *opt, *relaxationThresMat;
    double tolerance, inlierRatio;
    glores::PointSetGlobalRegistration optimizer;
    glores::VectorPoint2d pointsSrc;
    glores::VectorPoint2d pointsDst;
    Eigen::VectorXd poseLow(3), poseUpp(3), poseOpt(3), relaxationThres;
    double funcLow, funcUpp;
    bool enableCheapBoundQueue;

    //    mwSize m,n;
    //    mwSize nzmax;
    //    mwIndex *irs,*jcs,j,k;
    //    int cmplx;
    //    double *pr,*pi,*si,*sr;
    //    double percent_sparse;

    /* Check for proper number of input and output arguments */
    if (nrhs != 8) {
        mexErrMsgIdAndTxt("MATLAB:fulltosparse:invalidNumInputs",
                "Seven input arguments required.");
    }
    if (nlhs > 6) {
        mexErrMsgIdAndTxt("MATLAB:fulltosparse:maxlhs",
                "Too many output arguments.");
    }

    /* Check data type of input argument2  */

    if (!mxIsDouble(prhs[0])) {
        mexErrMsgIdAndTxt("MATLAB:mxgetnzmax:inputNotDouble",
                "Input argument 1 must be a double matrix.");
    }
    P = prhs[0];
    convertMxArrayToVectorPoint2(P, pointsSrc);

    if (!mxIsDouble(prhs[1])) {
        mexErrMsgIdAndTxt("MATLAB:mxgetnzmax:inputNotDouble",
                "Input argument 2 must be a 2 x N double matrix.");
    }
    Q = prhs[1];
    convertMxArrayToVectorPoint2(Q, pointsDst);

    if (!mxIsDouble(prhs[2])) {
        mexErrMsgIdAndTxt("MATLAB:mxgetnzmax:inputNotDouble",
                "Input argument 3 must be a double matrix.");
    }
    betamin = prhs[2];
    convertMxArrayToVector3(betamin, poseLow);

    if (!mxIsDouble(prhs[3])) {
        mexErrMsgIdAndTxt("MATLAB:mxgetnzmax:inputNotDouble",
                "Input argument 4 must be a double matrix.");
    }
    betamax = prhs[3];
    convertMxArrayToVector3(betamax, poseUpp);

    if (!(mxIsDouble(prhs[4]))) {
        mexErrMsgIdAndTxt("MATLAB:fulltosparse:inputNotDouble",
                "Input argument 5 must be of type double.");
    }
    tolerance = *mxGetPr(prhs[4]);

    if (!(mxIsDouble(prhs[5]))) {
        mexErrMsgIdAndTxt("MATLAB:fulltosparse:inputNotDouble",
                "Input argument 6 must be of type double.");
    }
    inlierRatio = *mxGetPr(prhs[5]);

    if (!(mxIsDouble(prhs[6]))) {
        mexErrMsgIdAndTxt("MATLAB:fulltosparse:inputNotDouble",
                "Input argument 7 must be a double matrix.");
    }
    relaxationThresMat = prhs[6];
    convertMxArrayToVector3(relaxationThresMat, relaxationThres);

    if (!(mxIsLogical(prhs[7]))) {
        mexErrMsgIdAndTxt("MATLAB:fulltosparse:inputNotDouble",
                "Input argument 8 must be a boolean flag.");
    }
    enableCheapBoundQueue = mxIsLogicalScalarTrue(prhs[7]);

    // Converts from mxArray to glores::VectorPoint2d

    optimizer.setPointSource(pointsSrc);
    optimizer.setPointDestination(pointsDst);
    optimizer.setYTolerance(tolerance);
    optimizer.enableYTolerance(true);
    optimizer.setInlierRatio(inlierRatio);
    optimizer.setRelaxationThres(relaxationThres);
    optimizer.enableCheapBoundQueueOn(enableCheapBoundQueue);
    optimizer.setVerboseLevel(0);

    optimizer.findGlobalMin(poseLow, poseUpp, poseOpt, funcLow, funcUpp);

    convertVectorToMxArray(poseOpt, plhs[0]);

    if (nlhs > 1) {
        convertDoubleToMxArray(funcLow, plhs[1]);
    }

    if (nlhs > 2) {
        convertDoubleToMxArray(funcUpp, plhs[2]);
    }

    if (nlhs > 3) {
        convertDoubleToMxArray(optimizer.getNodeNum(), plhs[3]);
    }

    if (nlhs > 4) {
        convertStatToMxArray(optimizer.getBoundStats(), plhs[4], 0); // lower bound
    }

    if (nlhs > 5) {
        convertStatToMxArray(optimizer.getBoundStats(), plhs[5], 1); // upper bound
    }

    if (nlhs > 6) {
        convertStatToMxArray(optimizer.getBoundStats(), plhs[6], 2); // relaxationBoundUsed
    }
}
