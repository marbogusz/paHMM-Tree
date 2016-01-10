//==============================================================================
// Pair-HMM phylogenetic tree estimator
// 
// Copyright (c) 2015 Marcin Bogusz.
// 
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
//     This program is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//     GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses>.
//==============================================================================

/*
 * BrentOptimizer.cpp
 *
 *  Created on: Aug 5, 2015
 *      Author: marcin
 *
 *  Partially based on John D. Cook's implementation
 *  Copyright (c) John D. Cook
 *  http://www.codeproject.com/Articles/30201/Optimizing-a-Function-of-One-Variable
 */

#include <cmath>
#include <cfloat>
#include "core/BrentOptimizer.hpp"

namespace EBC {


BrentOptimizer::BrentOptimizer(OptimizedModelParameters* mp,
		IOptimizable* opt, double accuracy) : accuracy(accuracy), omp(mp), target(opt)
{

DEBUG("Brent numerical optimizer with 1" << " parameter created");
//cerr << "DLIB optimizer init with " << paramsCount << " parameters" << endl;
}

double BrentOptimizer::objectiveFunction(double x)
{
	omp->setSingleDivergenceParam(0,x);
	return target->runIteration();
}


double BrentOptimizer::optimize()
{
    double leftEnd = leftBound;     // left bound should be about 0
    double rightEnd = rightBound;
    double epsilon = accuracy;
    double minLoc = omp->getDivergenceTime(0);

    double d, e, m, p, q, r, tol, t2, u, v, w, fu, fv, fw, fx;
    double golden_ratio = 0.5*(3.0 - sqrt(5.0));
    //double ZEPS = sqrt(DBL_EPSILON);
    double ZEPS = numeric_limits<double>::epsilon() * 0.001;

    double& a = leftEnd; double& b = rightEnd; double& x = minLoc;

    //v = w = x = a + c*(b - a);
    v = w = x;
    d = e = 0.0;
    fv = fw = fx = this->objectiveFunction(x);
    int counter = 0;

    for (counter=0; counter < Definitions::BrentMaxIter; counter++){
		m = 0.5*(a + b);
		tol = ZEPS + (fabs(x)*epsilon); t2 = 2.0*tol;
		// Check stopping criteria
		if (fabs(x - m) > t2 - 0.5*(b - a))
		{
			p = q = r = 0.0;
			if (fabs(e) > tol)
			{
				// fit parabola
				r = (x - w)*(fx - fv);
				q = (x - v)*(fx - fw);
				p = (x - v)*q - (x - w)*r;
				q = 2.0*(q - r);
				(q > 0.0) ? p = -p : q = -q;
				r = e; e = d;
			}
			if (fabs(p) < fabs(0.5*q*r) && p < q*(a - x) && p < q*(b - x))
			{
				// A parabolic interpolation step
				d = p/q;
				u = x + d;
				// f must not be evaluated too close to a or b
				if (u - a < t2 || b - u < t2)
					d = (x < m) ? tol : -tol;
			}
			else
			{
				// A golden section step
				e = (x < m) ? b : a;
				e -= x;
				d = golden_ratio*e;
			}
			// f must not be evaluated too close to x
			if (fabs(d) >= tol)
				u = x + d;
			else if (d > 0.0)
				u = x + tol;
			else
				u = x - tol;
			fu = this->objectiveFunction(u);
			// Update a, b, v, w, and x
			if (fu <= fx)
			{
				(u < x) ? b = x : a = x;
				v = w; fv = fw;
				w = x; fw = fx;
				x = u; fx = fu;
			}
			else
			{
				(u < x) ? a = u : b = u;
				if (fu <= fw || w == x)
				{
					v = w; fv = fw;
					w = u; fw = fu;
				}
				else if (fu <= fv || v == x || v == w)
				{
					v = u; fv = fu;
				}
			}
		}
		else
			break;
    }
    omp->setSingleDivergenceParam(0,minLoc);
    return  fx;
}


void BrentOptimizer::setTarget(IOptimizable* opt) {
	target = opt;
}

} /* namespace EBC */
