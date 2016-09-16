#ifndef BFGS_H_INCLUDED
#define BFGS_H_INCLUDED

#include "MeritFunction.h"
#include <iomanip>

class BFGS_Optimization{

    public:

        BFGS_Optimization(double tolerance,double maxStepSize,double optEps,std::vector<int>& outBasis);
        double minimize();
        double bestResult;
        MeritFunction meritFunction;

    private:

        double tol,eps;

        double stepMonitor,rho,alphaMax;

        double alpha0,alpha1,alpha2;
        double phi0,phi1,phi2;
        double phiPrime0,phiPrime1,phiPrime2;

        Eigen::VectorXd position,alphaPosition;
        Eigen::VectorXd gradient;
        Eigen::VectorXd p,s,y;

        Eigen::MatrixXd H,I;
        void setInitialH();
        void setInverseHessian();
        void setGradient();
        double alpha();
        void printResultReport();
        void printStepMonitor();
        double phi(double& a);
        int counter;
        double phiPrime(double& a);
        double zoom(double alphaLow,double alphaHigh,double phiLow,double phiHigh,double phiLowPrime);
        void setAlphaJ(double& alphaj,double& alphaLow,double& alphaHigh,double& phiLow,double& phiHigh,double& phiLowPrime);
        bool zoomGuard,quadInterpolationFailure,wolfeConditionFailure,maxStepSize,maxIterationGuard;
        double secondDerivativeTest;

};


#endif