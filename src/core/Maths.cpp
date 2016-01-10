//==============================================================================
// Pair-HMM phylogenetic tree estimator
// 
// Copyright (c) 2015 Marcin Bogusz.
// 
// Most of the methods of discerte gamma and eigen decomposition taken from
// Ziheng Yang's PAML package, Copyright (c) Ziheng Yang
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

#include "core/Maths.hpp"
#include <iostream>

namespace EBC
{

Maths::Maths()
{
	srand((unsigned int)time(0));
	z_rndu = rand();

}

double Maths::logSum(double a, double b, double c)
{
	double partialSum = logSum(a,b);
	return logSum(partialSum,c);
}

double Maths::logSum(double a, double b)
{
	if (a > b)
	{
		return a + log(1+exp(b-a));
	}
	else
	{
		return b + log(1+exp(a-b));
	}
}


double* Maths::matrixMultiply(double *matA, double *matB, int size)
{
	double* matResult = new double[size*size];

	for (int i=0; i<size; i++)
	{
		for (int j=0; j<size; j++)
		{
			matResult[i*size + j]=0;
			for (int k=0; k<size; k++)
			{
				matResult[i*size + j] += matA[i*size+ k] * matB[k *size + j];
			}
		}
	}
	return matResult;
}

double* Maths::matrixByDiagonalMultiply(double *matA, double *matDiag, int size)
{
	double* res = new double[size*size];
	for (int i=0; i<size; i++)
		for (int j=0; j<size; j++)
			res[i*size+j]  = matA[i*size+j] *  matDiag[j];
	return res;
}

void Maths::matrixByDiagonalMultiplyMutable(double *matA, double *matDiag, int size)
{
	for (int i=0; i<size; i++)
		for (int j=0; j<size; j++)
			matA[i*size+j] *=  matDiag[j];
}

void Maths::vectorMultiply(double* vector, double factor, int size)
{
	for(int i=0; i<size; i++)
		vector[i] *= factor;
}

void Maths::matrixAppend(double* matA, double* matB, int size)
{
	for (int i = 0; i < size*size; i++)
	{
		matA[i] += matB[i];
	}
}

void Maths::matrixScale(double* matrix, double factor, int size)
{
	this->vectorMultiply(matrix, factor, size*size);
}

double* Maths::expLambdaT(double* lambda, double t, int size)
{
	double *res = new double[size];
	for(int i=0; i<size; i++)
		res[i] = exp(lambda[i]*t);

	return res;
}

double Maths::getRandom(double lo=0.0, double hi=1.0)
{
	double divider = RAND_MAX/hi;
	return lo +  (rand()/divider);
}

void Maths::revMatrixExponentiation(double* Q, double* P)
{

}

int Maths::eigenRealSym(double A[], int n, double Root[], double work[])
{
/* This finds the eigen solution of a real symmetrical matrix A[n*n].  In return,
   A has the right vectors and Root has the eigenvalues.
   work[n] is the working space.
   The matrix is first reduced to a tridiagonal matrix using HouseholderRealSym(),
   and then using the QL algorithm with implicit shifts.

   Adapted from routine tqli in Numerical Recipes in C, with reference to LAPACK
   Ziheng Yang, 23 May 2001
*/
   int status=0;
   HouseholderRealSym(A, n, Root, work);
   status = EigenTridagQLImplicit(Root, work, n, A);
   EigenSort(Root, A, n);

   return(status);
}

int Maths::eigenQREV (double Q[], double pi[], int n,
               double Root[], double U[], double V[], double spacesqrtpi[])
{
/*
   This finds the eigen solution of the rate matrix Q for a time-reversible
   Markov process, using the algorithm for a real symmetric matrix.
   Rate matrix Q = S * diag{pi} = U * diag{Root} * V,
   where S is symmetrical, all elements of pi are positive, and U*V = I.
   space[n] is for storing sqrt(pi).

   [U 0] [Q_0 0] [U^-1 0]    [Root  0]
   [0 I] [0   0] [0    I]  = [0     0]

   Ziheng Yang, 25 December 2001 (ref is CME/eigenQ.pdf)
*/
   int i,j, inew, jnew, nnew, status;
   double *pi_sqrt=spacesqrtpi, small=1e-100;

   for(j=0,nnew=0; j<n; j++)
      if(pi[j]>small)
         pi_sqrt[nnew++] = sqrt(pi[j]);

   /* store in U the symmetrical matrix S = sqrt(D) * Q * sqrt(-D) */

   if(nnew==n) {
      for(i=0; i<n; i++)
         for(j=0,U[i*n+i] = Q[i*n+i]; j<i; j++)
            U[i*n+j] = U[j*n+i] = (Q[i*n+j] * pi_sqrt[i]/pi_sqrt[j]);

      status=eigenRealSym(U, n, Root, V);
      for(i=0; i<n; i++) for(j=0; j<n; j++)  V[i*n+j] = U[j*n+i] * pi_sqrt[j];
      for(i=0; i<n; i++) for(j=0; j<n; j++)  U[i*n+j] /= pi_sqrt[i];
   }
   else {
      for(i=0,inew=0; i<n; i++) {
         if(pi[i]>small) {
            for(j=0,jnew=0; j<i; j++)
               if(pi[j]>small) {
                  U[inew*nnew+jnew] = U[jnew*nnew+inew]
                                    = Q[i*n+j] * pi_sqrt[inew]/pi_sqrt[jnew];
                  jnew++;
               }
            U[inew*nnew+inew] = Q[i*n+i];
            inew++;
         }
      }

      status=eigenRealSym(U, nnew, Root, V);

      for(i=n-1,inew=nnew-1; i>=0; i--)   /* construct Root */
         Root[i] = (pi[i]>small ? Root[inew--] : 0);
      for(i=n-1,inew=nnew-1; i>=0; i--) {  /* construct V */
         if(pi[i]>small) {
            for(j=n-1,jnew=nnew-1; j>=0; j--)
               if(pi[j]>small) {
                  V[i*n+j] = U[jnew*nnew+inew]*pi_sqrt[jnew];
                  jnew--;
               }
               else
                  V[i*n+j] = (i==j);
            inew--;
         }
         else
            for(j=0; j<n; j++)  V[i*n+j] = (i==j);
      }
      for(i=n-1,inew=nnew-1; i>=0; i--) {  /* construct U */
         if(pi[i]>small) {
            for(j=n-1,jnew=nnew-1;j>=0;j--)
               if(pi[j]>small) {
                  U[i*n+j] = U[inew*nnew+jnew]/pi_sqrt[inew];
                  jnew--;
               }
               else
                  U[i*n+j] = (i==j);
            inew--;
         }
         else
            for(j=0;j<n;j++)
               U[i*n+j] = (i==j);
      }
   }
   return 0;
}

void Maths::HouseholderRealSym(double a[], int n, double d[], double e[])
{
/* This uses HouseholderRealSym transformation to reduce a real symmetrical matrix
   a[n*n] into a tridiagonal matrix represented by d and e.
   d[] is the diagonal (eigends), and e[] the off-diagonal.
*/
   int m,k,j,i;
   double scale,hh,h,g,f;

   for (i=n-1;i>=1;i--) {
      m=i-1;
      h=scale=0;
      if (m > 0) {
         for (k=0;k<=m;k++)
            scale += fabs(a[i*n+k]);
         if (scale == 0)
            e[i]=a[i*n+m];
         else {
            for (k=0;k<=m;k++) {
               a[i*n+k] /= scale;
               h += a[i*n+k]*a[i*n+k];
            }
            f=a[i*n+m];
            g=(f >= 0 ? -sqrt(h) : sqrt(h));
            e[i]=scale*g;
            h -= f*g;
            a[i*n+m]=f-g;
            f=0;
            for (j=0;j<=m;j++) {
               a[j*n+i]=a[i*n+j]/h;
               g=0;
               for (k=0;k<=j;k++)
                  g += a[j*n+k]*a[i*n+k];
               for (k=j+1;k<=m;k++)
                  g += a[k*n+j]*a[i*n+k];
               e[j]=g/h;
               f += e[j]*a[i*n+j];
            }
            hh=f/(h*2);
            for (j=0;j<=m;j++) {
               f=a[i*n+j];
               e[j]=g=e[j]-hh*f;
               for (k=0;k<=j;k++)
                  a[j*n+k] -= (f*e[k]+g*a[i*n+k]);
            }
         }
      }
      else
         e[i]=a[i*n+m];
      d[i]=h;
   }
   d[0]=e[0]=0;

   /* Get eigenvectors */
   for (i=0;i<n;i++) {
      m=i-1;
      if (d[i]) {
         for (j=0;j<=m;j++) {
            g=0;
            for (k=0;k<=m;k++)
               g += a[i*n+k]*a[k*n+j];
            for (k=0;k<=m;k++)
               a[k*n+j] -= g*a[k*n+i];
         }
      }
      d[i]=a[i*n+i];
      a[i*n+i]=1;
      for (j=0;j<=m;j++) a[j*n+i]=a[i*n+j]=0;
   }
}

void Maths::EigenSort(double d[], double U[], int n)
{
/* this sorts the eigen values d[] and rearrange the (right) eigen vectors U[]
*/
   int k,j,i;
   double p;

   for (i=0; i<n-1; i++) {
      p = d[k=i];
      for (j=i+1; j<n; j++)
         if (d[j] >= p) p = d[k=j];
      if (k != i) {
         d[k] = d[i];
         d[i] = p;
         for (j=0;j<n;j++) {
            p = U[j*n+i];
            U[j*n+i] = U[j*n+k];
            U[j*n+k] = p;
         }
      }
   }
}

int Maths::EigenTridagQLImplicit(double d[], double e[], int n, double z[])
{
/* This finds the eigen solution of a tridiagonal matrix represented by d and e.
   d[] is the diagonal (eigenvalues), e[] is the off-diagonal
   z[n*n]: as input should have the identity matrix to get the eigen solution of the
   tridiagonal matrix, or the output from HouseholderRealSym() to get the
   eigen solution to the original real symmetric matrix.
   z[n*n]: has the orthogonal matrix as output

   Adapted from routine tqli in Numerical Recipes in C, with reference to
   LAPACK fortran code.
   Ziheng Yang, May 2001
*/
   int m,j,iter,niter=30, status=0, i,k;
   double s,r,p,g,f,dd,c,b, aa,bb;

   for (i=1;i<n;i++) e[i-1]=e[i];  e[n-1]=0;
   for (j=0;j<n;j++) {
      iter=0;
      do {
         for (m=j;m<n-1;m++) {
            dd=fabs(d[m])+fabs(d[m+1]);
            if (fabs(e[m])+dd == dd) break;  /* ??? */
         }
         if (m != j) {
            if (iter++ == niter) {
               status=-1;
               break;
            }
            g=(d[j+1]-d[j])/(2*e[j]);

            /* r=pythag(g,1); */

            if((aa=fabs(g))>1)  r=aa*sqrt(1+1/(g*g));
            else                r=sqrt(1+g*g);

            g=d[m]-d[j]+e[j]/(g+((g) >= 0.0 ? fabs(r) : -fabs(r)));
            s=c=1;
            p=0;
            for (i=m-1;i>=j;i--) {
               f=s*e[i];
               b=c*e[i];

               /*  r=pythag(f,g);  */
               aa=fabs(f); bb=fabs(g);
               if(aa>bb)       { bb/=aa;  r = aa*sqrt(1+bb*bb); }
               else if(bb==0)             r = 0;
               else            { aa/=bb;  r = bb*sqrt(1+aa*aa); }

               e[i+1]=r;
               if (r == 0) {
                  d[i+1] -= p;
                  e[m]=0;
                  break;
               }
               s=f/r;
               c=g/r;
               g=d[i+1]-p;
               r=(d[i]-g)*s+2*c*b;
               d[i+1]=g+(p=s*r);
               g=c*r-b;
               for (k=0;k<n;k++) {
                  f=z[k*n+i+1];
                  z[k*n+i+1]=s*z[k*n+i]+c*f;
                  z[k*n+i]=c*z[k*n+i]-s*f;
               }
            }
            if (r == 0 && i >= j) continue;
            d[j]-=p; e[j]=g; e[m]=0;
         }
      } while (m != j);
   }
   return(status);
}

/*From PAML by Ziheng Yang*/
double Maths::rndu (void)
{
/* 32-bit integer assumed.
   From Ripley (1987) p. 46 or table 2.4 line 2.
   This may return 0 or 1, which can be a problem.
*/
   z_rndu = z_rndu*69069 + 1;
   if(z_rndu==0 || z_rndu==4294967295)  z_rndu = 13;
   return z_rndu/4294967295.0;
}

int  Maths::DiscreteGamma (double freqK[], double rK[], double alpha, double beta, int K, int UseMedian)
{
/* discretization of G(alpha, beta) with equal proportions in each category.
*/
   int i;
   double t, mean=alpha/beta, lnga1;

   if(UseMedian) {   /* median */
      for(i=0; i<K; i++) rK[i] = QuantileChi2((i*2.+1)/(2.*K), (2.0*alpha))/ (2.0*beta);
      for(i=0,t=0; i<K; i++) t += rK[i];
      for(i=0; i<K; i++) rK[i] *= mean*K/t;   /* rescale so that the mean is alpha/beta. */
   }
   else {            /* mean */
      lnga1 = LnGamma(alpha+1);
      for (i=0; i<K-1; i++) /* cutting points, Eq. 9 */
         freqK[i] = QuantileChi2((i+1.0)/K, 2.0 * alpha)/ (2.0*beta);
      for (i=0; i<K-1; i++) /* Eq. 10 */
         freqK[i] = IncompleteGamma(freqK[i]*beta, alpha+1, lnga1);
      rK[0] = freqK[0]*mean*K;
      for (i=1; i<K-1; i++)  rK[i] = (freqK[i] - freqK[i-1])*mean*K;
      rK[K-1] = (1-freqK[K-2])*mean*K;
   }

   for (i=0; i<K; i++) freqK[i] = 1.0/K;

   return (0);
}

double Maths::LnGamma (double x)
{
/* returns ln(gamma(x)) for x>0, accurate to 10 decimal places.
   Stirling's formula is used for the central polynomial part of the procedure.

   Pike MC & Hill ID (1966) Algorithm 291: Logarithm of the gamma function.
   Communications of the Association for Computing Machinery, 9:684
*/
   double f=0, fneg=0, z, lng;
   int nx=(int)x;

   if((double)nx==x && nx>=0 && nx<=11)
      lng = log((double)factorial(nx-1));
   else {
      if(x<=0) {
         if((int)x-x==0) {return(-1); }
         for (fneg=1; x<0; x++) fneg /= x;
         if(fneg<0)
            std::cerr << "strange!! check lngamma" << std::endl;
         fneg = log(fneg);
      }
      if (x<7) {
         f = 1;
         z = x-1;
         while (++z<7)
            f *= z;
         x = z;
         f = -log(f);
      }
      z = 1/(x*x);
      lng = fneg+ f + (x-0.5)*log(x) - x + .918938533204673
             + (((-.000595238095238*z+.000793650793651)*z-.002777777777778)*z
                  +.083333333333333)/x;
   }
   return  lng;
}

long Maths::factorial (int n)
{
   long f=1, i;
   for (i=2; i<=(long)n; i++) f *= i;
   return (f);
}

double Maths::QuantileChi2 (double prob, double v)
{
/* returns z so that Prob{x<z}=prob where x is Chi2 distributed with df=v
   returns -1 if in error.   0.000002<prob<0.999998
   RATNEST FORTRAN by
       Best DJ & Roberts DE (1975) The percentage points of the
       Chi2 distribution.  Applied Statistics 24: 385-388.  (AS91)
   Converted into C by Ziheng Yang, Oct. 1993.
*/
   double e=.5e-6, aa=.6931471805, p=prob, g, small=1e-6;
   double xx, c, ch, a=0,q=0,p1=0,p2=0,t=0,x=0,b=0,s1,s2,s3,s4,s5,s6;

   if (p<small)   return(0);
   if (p>1-small) return(9999);
   if (v<=0)      return (-1);

   g = LnGamma (v/2);
   xx=v/2;   c=xx-1;
   if (v >= -1.24*log(p)) goto l1;

   ch=pow((p*xx*exp(g+xx*aa)), 1/xx);
   if (ch-e<0) return (ch);
   goto l4;
l1:
   if (v>.32) goto l3;
   ch=0.4;   a=log(1-p);
l2:
   q=ch;  p1=1+ch*(4.67+ch);  p2=ch*(6.73+ch*(6.66+ch));
   t=-0.5+(4.67+2*ch)/p1 - (6.73+ch*(13.32+3*ch))/p2;
   ch-=(1-exp(a+g+.5*ch+c*aa)*p2/p1)/t;
   if (fabs(q/ch-1)-.01 <= 0) goto l4;
   else                       goto l2;

l3:
   x = QuantileNormal(p);
   p1 = 0.222222/v;
   ch = v*pow((x*sqrt(p1)+1-p1), 3.0);
   if (ch>2.2*v+6)
      ch = -2*(log(1-p)-c*log(.5*ch)+g);
l4:
   q=ch;   p1=.5*ch;
   if ((t=IncompleteGamma (p1, xx, g))<0)
      std::cerr << "\nIncompleteGamma";
   p2=p-t;
   t=p2*exp(xx*aa+g+p1-c*log(ch));
   b=t/ch;  a=0.5*t-b*c;

   s1=(210+a*(140+a*(105+a*(84+a*(70+60*a))))) / 420;
   s2=(420+a*(735+a*(966+a*(1141+1278*a))))/2520;
   s3=(210+a*(462+a*(707+932*a)))/2520;
   s4=(252+a*(672+1182*a)+c*(294+a*(889+1740*a)))/5040;
   s5=(84+264*a+c*(175+606*a))/2520;
   s6=(120+c*(346+127*c))/5040;
   ch+=t*(1+0.5*t*s1-b*c*(s1-b*(s2-b*(s3-b*(s4-b*(s5-b*s6))))));
   if (fabs(q/ch-1) > e) goto l4;

   return (ch);
}
double Maths::QuantileNormal (double prob)
{
/* returns z so that Prob{x<z}=prob where x ~ N(0,1) and (1e-12)<prob<1-(1e-12)
   returns (-9999) if in error
   Odeh RE & Evans JO (1974) The percentage points of the normal distribution.
   Applied Statistics 22: 96-97 (AS70)

   Newer methods:
     Wichura MJ (1988) Algorithm AS 241: the percentage points of the
       normal distribution.  37: 477-484.
     Beasley JD & Springer SG  (1977).  Algorithm AS 111: the percentage
       points of the normal distribution.  26: 118-121.
*/
   double a0=-.322232431088, a1=-1, a2=-.342242088547, a3=-.0204231210245;
   double a4=-.453642210148e-4, b0=.0993484626060, b1=.588581570495;
   double b2=.531103462366, b3=.103537752850, b4=.0038560700634;
   double y, z=0, p=prob, p1;

   p1 = (p<0.5 ? p : 1-p);
   if (p1<1e-20) z=999;
   else {
      y = sqrt (log(1/(p1*p1)));
      z = y + ((((y*a4+a3)*y+a2)*y+a1)*y+a0) / ((((y*b4+b3)*y+b2)*y+b1)*y+b0);
   }
   return (p<0.5 ? -z : z);
}

double Maths::IncompleteGamma (double x, double alpha, double ln_gamma_alpha)
{
/* returns the incomplete gamma ratio I(x,alpha) where x is the upper
           limit of the integration and alpha is the shape parameter.
   returns (-1) if in error
   ln_gamma_alpha = ln(Gamma(alpha)), is almost redundant.
   (1) series expansion,     if (alpha>x || x<=1)
   (2) continued fraction,   otherwise
   RATNEST FORTRAN by
   Bhattacharjee GP (1970) The incomplete gamma integral.  Applied Statistics,
   19: 285-287 (AS32)
*/
   int i;
   double p=alpha, g=ln_gamma_alpha;
   double accurate=1e-10, overflow=1e60;
   double factor, gin=0, rn=0, a=0,b=0,an=0,dif=0, term=0, pn[6];

   if (x==0) return (0);
   if (x<0 || p<=0) return (-1);

   factor=exp(p*log(x)-x-g);
   if (x>1 && x>=p) goto l30;
   /* (1) series expansion */
   gin=1;  term=1;  rn=p;
 l20:
   rn++;
   term *= x/rn;   gin += term;
   if (term > accurate) goto l20;
   gin *= factor/p;
   goto l50;
 l30:
   /* (2) continued fraction */
   a = 1-p;   b = a+x+1;  term = 0;
   pn[0] = 1;  pn[1] = x;  pn[2] = x+1;  pn[3] = x*b;
   gin = pn[2]/pn[3];
 l32:
   a++;
   b += 2;
   term++;
   an = a*term;
   for (i=0; i<2; i++)
      pn[i+4] = b*pn[i+2] - an*pn[i];
   if (pn[5] == 0) goto l35;
   rn = pn[4]/pn[5];
   dif = fabs(gin-rn);
   if (dif > accurate) goto l34;
   if (dif <= accurate*rn) goto l42;
 l34:
   gin = rn;
 l35:
   for (i=0; i<4; i++) pn[i] = pn[i+2];
   if (fabs(pn[4]) < overflow) goto l32;
   for (i=0; i<4; i++) pn[i] /= overflow;
   goto l32;
 l42:
   gin = 1-factor*gin;

 l50:
   return (gin);
}

} /* namespace EBC */
