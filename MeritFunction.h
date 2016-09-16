#ifndef MERITFUNCTION_H_INCLUDED
#define MERITFUNCTION_H_INCLUDED


#include "AncillaAugment.h"
#include "PUA.h"
#include "LOTransform.h"
#include <fstream>
#include <iomanip>

class MeritFunction{

    public:

        MeritFunction();
        void setMeritFunction(double optEps,std::vector<int>& outBasis);
        double f(Eigen::VectorXd& position);
        int funcDimension;
        void printReport(Eigen::VectorXd& position);
        Eigen::VectorXd setInitialPosition();
        bool validBasis;

    private:

        int photons,modes,ancillaPhotons,ancillaModes,measModes,measOutcome,compSubspaceDim,numbMeasOutcomes;
        Eigen::MatrixXi computationalBasisIn,computationalBasisOut;
        PUA MeasAssistOp1;
        Eigen::MatrixXcd IdealOp,IdealOpConverted,totalOpConverted;
        double IdealFidelityNorm,totalOpFidelityNorm,fidelity,successProbability,optimizationEpsilon;
        int nonZeroCoord1,nonZeroCoord2;
        int g(int n,int m);
        double doublefactorial(int x);
        void setIdealOpConverted();
        Eigen::MatrixXi generateSubBasisVector(int subPhotons, int subModes);
        Eigen::MatrixXi generateBasisVector(int subPhotons,int subModes, int subMeasureModes);
        bool isEqual(Eigen::ArrayXi a1,Eigen::ArrayXi a2);
        Eigen::ArrayXi compBasisAddressIn,compBasisAddressOut;
        void setCompBasisAddress(Eigen::MatrixXi& basisVector);
        void orderIdealOpConvertedColumns();
        void setTotalOpConverted();
        void setNumbMeasOutcomes();
        void setIdealFidelityNorm();
        void setTotalOpFidelityNorm();
        void setFidelity();
        void setSuccessProbability();
        void setNonZeroIdealCoordinates();
        Eigen::ArrayXXcd genUnitary(Eigen::ArrayXd& a);
        Eigen::MatrixXcd genHermitian(Eigen::ArrayXd& a);
        Eigen::ArrayXXcd matrixExp(Eigen::MatrixXcd X);
        Eigen::ArrayXXcd matrixLog(Eigen::MatrixXcd X);
        Eigen::ArrayXd convertHermittoA(Eigen::ArrayXXcd& H);
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> saes;
        Eigen::MatrixXi genBasisOut(std::vector<int>& outBasis);


};

#endif // MERITFUNCTION_H_INCLUDED
