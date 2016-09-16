#include "MeritFunction.h"

#define INITIAL_CONDITION_RANDOM_DEGREE 200
#define PI 3.14159265359


/** ====== Restrict Unitarity of Linear Operators (also change in PUA.cpp and LOTransform.cpp  ======================= */

#define RESTRICT_TO_UNITARY

//#define ALLOW_ARBITRARY_VACUUM_MODES

/** ================================================================================================================== */


void MeritFunction::setMeritFunction(double optEps,std::vector<int>& outBasis){

    photons = 2;
    modes = 4;

    ancillaPhotons = 0;
    ancillaModes = 0;

    measModes = 0;
    measOutcome = 0;

    compSubspaceDim = 4;

    computationalBasisIn.resize(compSubspaceDim,modes);
    computationalBasisOut.resize(compSubspaceDim,modes);

    std::ifstream infile("InBasis.dat");

    for(int i=0;i<compSubspaceDim;i++){

        for(int j=0;j<modes;j++){

            infile >> computationalBasisIn(i,j);

        }

    }

    infile.close();

    computationalBasisOut = genBasisOut(outBasis);


    std::ofstream outfile("BasisCheck.dat");
    outfile << computationalBasisIn << std::endl;
    outfile << "TO" << std::endl;
    outfile << computationalBasisOut << std::endl << std::endl;
    outfile.close();

    IdealOp.resize(compSubspaceDim,compSubspaceDim);

    std::complex<double> I(0.0,1.0);

    IdealOp <<  0.520441 + 0.801184 * I ,0.0,-0.260429 + 0.139366*I,0.0,
                    0.0,            0.520441 + 0.801184 * I, 0.0, -0.260429 + 0.139366*I,
                0.123947 + 0.268111*I,      0.0,        0.898012 - 0.326081*I, 0.0,
                0.0,0.123947 + 0.268111*I,      0.0,        0.898012 - 0.326081*I;


    MeasAssistOp1.setPUA(photons,modes,ancillaPhotons,ancillaModes,measModes,measOutcome);

    setIdealOpConverted();

    setNonZeroIdealCoordinates();

    setIdealFidelityNorm();

    #ifdef RESTRICT_TO_UNITARY

        funcDimension = (modes + ancillaModes) * (modes + ancillaModes) + 2*g(ancillaPhotons,ancillaModes);

    #endif // RESTRICT_TO_UNITARY

    #ifdef ALLOW_ARBITRARY_VACUUM_MODES

        funcDimension = 2 * (modes + ancillaModes) * (modes + ancillaModes) + 2*g(ancillaPhotons,ancillaModes);

    #endif // ALLOW_ARBITRARY_VACUUM_MODES

    setNumbMeasOutcomes();

    optimizationEpsilon = optEps;

    validBasis = false;

    return;

}


void MeritFunction::printReport(Eigen::VectorXd& position){

    MeasAssistOp1.setQuantumOperator(position);

    std::cout << std::endl << std::endl;
    std::cout << "RESULT FIDELITY: \n" << fidelity << std::endl << std::endl;


    if(fidelity > 1.0 - 1.0e-6){

        std::ofstream outfile("Successful Basis Change.dat",std::ofstream::app);
        outfile << computationalBasisIn << std::endl << std::endl;
        outfile << "TO" << std::endl << std::endl;
        outfile << computationalBasisOut << std::endl << std::endl;
        outfile.close();

        std::ofstream outfile2("Successful Outbasis List.dat",std::ofstream::app);
        outfile2 << computationalBasisOut << std::endl << std::endl;
        outfile2.close();

        validBasis = true;

    }

    return;

}


double MeritFunction::f(Eigen::VectorXd& position){

    MeasAssistOp1.setQuantumOperator(position);

    setTotalOpConverted();

    setTotalOpFidelityNorm();

    setFidelity();

    setSuccessProbability();

/** Checking trace-preservingness of the Kraus Operators =========================================== */

//    Eigen::MatrixXcd TraceTest = Eigen::MatrixXcd::Zero(compSubspaceDim,compSubspaceDim);
//
//    for(int i=0;i<numbMeasOutcomes;i++){
//        std::cout << "i: " << i << std::endl;
//        MeasAssistOp1.setPUA(photons,modes,ancillaPhotons,ancillaModes,measModes,i);
//        MeasAssistOp1.setQuantumOperator(position);
//        totalOpConverted.resize(MeasAssistOp1.TotalOp.rows(),compSubspaceDim);
//        setTotalOpConverted();
//        TraceTest += totalOpConverted.conjugate().transpose() * totalOpConverted;
//        std::cout << "TraceTest: \n" <<  TraceTest << std::endl << std::endl;
//
//    }

/** ================================================================================================ */

    return - (fidelity + optimizationEpsilon * successProbability);

}


Eigen::MatrixXi MeritFunction::genBasisOut(std::vector<int>& outBasis){

    Eigen::MatrixXi output(compSubspaceDim,modes);

    Eigen::MatrixXi basisVector = generateBasisVector(photons,modes,1);

    for(int i=0;i<compSubspaceDim;i++){

        output.row(i) = basisVector.row(outBasis.at(i));

    }

    return output;

}

void MeritFunction::setTotalOpConverted(){

    for(int i=0;i<compSubspaceDim;i++){

        totalOpConverted.col(i) = MeasAssistOp1.TotalOp.col(compBasisAddressIn(i));

    }

    return;

}

void MeritFunction::setSuccessProbability(){

    successProbability = norm( totalOpConverted(nonZeroCoord1,nonZeroCoord2)/IdealOpConverted(nonZeroCoord1,nonZeroCoord2) );

    return;

}

void MeritFunction::setFidelity(){

    fidelity = norm((1.0/compSubspaceDim) * ((IdealOpConverted.conjugate().transpose() * totalOpConverted).trace()));

    fidelity /= totalOpFidelityNorm;

    fidelity /= IdealFidelityNorm;

    return;

}

void MeritFunction::setTotalOpFidelityNorm(){

    totalOpFidelityNorm = sqrt(norm((1.0/compSubspaceDim) * (totalOpConverted.conjugate().transpose() *totalOpConverted).trace()));

    return;

}

void MeritFunction::setNonZeroIdealCoordinates(){

    for(int i=0;i<IdealOpConverted.rows();i++){

        for(int j=0;j<IdealOpConverted.cols();j++){

            if(sqrt(norm(IdealOpConverted(i,j))) > 1e-10){

                nonZeroCoord1 = i;
                nonZeroCoord2 = j;

                return;

            }

        }

    }

}

Eigen::VectorXd MeritFunction::setInitialPosition(){

    std::complex<double> I(0.0,1.0);

    #ifdef RESTRICT_TO_UNITARY

        Eigen::VectorXd output = Eigen::VectorXd::Random(funcDimension);
        output *= PI;
        Eigen::MatrixXcd U = Eigen::MatrixXcd::Identity(modes + ancillaModes,modes + ancillaModes);

        for(int i=0;i<INITIAL_CONDITION_RANDOM_DEGREE;i++){
            Eigen::ArrayXd a = Eigen::ArrayXd::Random((modes + ancillaModes)*(modes + ancillaModes));
            a *= 2000*PI;
            Eigen::MatrixXcd Utemp((modes + ancillaModes),(modes + ancillaModes));
            Utemp = genUnitary(a);
            U *= Utemp;
        }

        Eigen::ArrayXXcd H(modes + ancillaModes,modes + ancillaModes);
        H = matrixLog(U)/I;

        Eigen::ArrayXd a = convertHermittoA(H);

        for(int i=0;i<(modes + ancillaModes)*(modes + ancillaModes);i++){
            output(i) = a(i);
        }

    #endif // RESTRICT_TO_UNITARY

    #ifdef ALLOW_ARBITRARY_VACUUM_MODES

        Eigen::VectorXd output = Eigen::VectorXd::Random(funcDimension);
        output *= PI;
        Eigen::MatrixXcd U = Eigen::MatrixXcd::Random(modes + ancillaModes,modes + ancillaModes);

        saes.compute(U.conjugate().transpose() * U);

        U *= 1.0/sqrt(saes.eigenvalues().maxCoeff());

        int k=0;

        for(int i=0;i<modes+ancillaModes;i++){

            for(int j=0;j<modes+ancillaModes;j++){

                output(k) = real(U(i,j));
                k++;

                output(k) = imag(U(i,j));
                k++;

            }

        }

    #endif // ALLOW_ARBITRARY_VACUUM_MODES

    return output;

}

Eigen::ArrayXd MeritFunction::convertHermittoA(Eigen::ArrayXXcd& H){
    int Hsize = H.rows();
    Eigen::ArrayXd output(Hsize*Hsize);
    int rowIndex = 0;
    int outputIndex =0;
    for(int i=0;i<Hsize;i++){
        output(outputIndex) = real(H(rowIndex,rowIndex));
        outputIndex++;
        for(int j=rowIndex+1;j<Hsize;j++){
            output(outputIndex) = sqrt(norm(H(rowIndex,j)));
            outputIndex++;
            output(outputIndex) = arg(H(rowIndex,j));
            outputIndex++;
        }
        rowIndex++;
    }
    return output;
}


Eigen::ArrayXXcd MeritFunction::genUnitary(Eigen::ArrayXd& a){
    return matrixExp(genHermitian(a));
}

Eigen::MatrixXcd MeritFunction::genHermitian(Eigen::ArrayXd& a){

    std::complex<double> I(0.0,1.0);
    int Hsize = sqrt(a.size());
    Eigen::MatrixXcd m(Hsize,Hsize);
    int extractIndex=0;                                     //REWRITE THIS FUNCTION IT NEEDS TO BE EFFICIENT- EIGEN SHOULD HAVE A STANDARD ONE
    for(int i=0;i<Hsize;i++){
        m(i,i)=a(extractIndex);
        extractIndex++;
        for(int j=i;j<Hsize;j++){
            if(i!=j){
                //m(i,j)=a(extractIndex)+I*a(extractIndex+1);
                //m(j,i)=a(extractIndex)-I*a(extractIndex+1);
                m(i,j) = a(extractIndex) * exp(I*a(extractIndex+1));
                m(j,i) = a(extractIndex) * exp(-I*a(extractIndex+1));
                extractIndex++;
                extractIndex++;
            }
        }
    }

    return m;
}

Eigen::ArrayXXcd MeritFunction::matrixExp(Eigen::MatrixXcd X){

    int matrixSize;
    matrixSize = X.rows();
    std::complex<double> I(0.0,1.0);

                                                                //THIS NEEDS TO BE AS EFFICIENT AS POSSIBLE
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> ces;
    ces.compute(X);
    Eigen::VectorXd evalues=ces.eigenvalues();
    Eigen::MatrixXcd evectors=(ces.eigenvectors());
    Eigen::MatrixXcd cevectors=evectors.conjugate();
    Eigen::ArrayXXcd sylvester[matrixSize];

    for(int i=0;i < matrixSize;i++){
        sylvester[i].resize(matrixSize,matrixSize);
        for(int m=0; m<matrixSize;m++){
            for(int n=0;n<matrixSize;n++){
                sylvester[i](n,m)=evectors(n,i)*cevectors(m,i);
            }
        }
    }

    Eigen::ArrayXXcd result(matrixSize,matrixSize);
    result = exp(I*evalues(0))*sylvester[0];
    for(int j=1;j<matrixSize;j++){
        result=result+exp(I*evalues(j))*sylvester[j];
    }

    return result;
}

Eigen::ArrayXXcd MeritFunction::matrixLog(Eigen::MatrixXcd X){

    int matrixSize;
    matrixSize = X.rows();

    Eigen::ComplexEigenSolver<Eigen::MatrixXcd> ces;
    ces.compute(X);
    Eigen::VectorXcd evalues=ces.eigenvalues();
    Eigen::MatrixXcd evectors=(ces.eigenvectors());
    Eigen::MatrixXcd cevectors=evectors.conjugate();
    Eigen::ArrayXXcd sylvester[matrixSize];

    for(int i=0;i < matrixSize;i++){
        sylvester[i].resize(matrixSize,matrixSize);
        for(int m=0; m<matrixSize;m++){
            for(int n=0;n<matrixSize;n++){
                sylvester[i](n,m)=evectors(n,i)*cevectors(m,i);
            }
        }
    }

    Eigen::ArrayXXcd result(matrixSize,matrixSize);
    result = log(evalues(0))*sylvester[0];
    for(int j=1;j<matrixSize;j++){
        result = result + log(evalues(j))*sylvester[j];
    }

    return result;
}

void MeritFunction::setIdealFidelityNorm(){

    IdealFidelityNorm = sqrt(norm((1.0/compSubspaceDim) * (IdealOpConverted.conjugate().transpose() * IdealOpConverted).trace()));

    return;

}


void MeritFunction::setNumbMeasOutcomes(){

    numbMeasOutcomes = 0;

    if(measModes == 0){

        numbMeasOutcomes = 1;

        return;

    }

    for(int i=0;i<=photons+ancillaPhotons;i++){

        numbMeasOutcomes += g(i,measModes);

    }

    return;

}


void MeritFunction::setIdealOpConverted(){

    Eigen::MatrixXi basisVector = generateBasisVector(photons,modes,1);


    setCompBasisAddress(basisVector);


    IdealOpConverted = Eigen::MatrixXcd::Zero(g(photons,modes),compSubspaceDim);
    totalOpConverted = Eigen::MatrixXcd::Zero(g(photons,modes),compSubspaceDim);


    for(int i=0;i<compSubspaceDim;i++){

        IdealOpConverted.row(compBasisAddressOut(i)) = IdealOp.row(i);

    }

    orderIdealOpConvertedColumns();


    return;

}

void MeritFunction::orderIdealOpConvertedColumns(){

    while(true){

        for(int i=0;i<compSubspaceDim-1;i++){

            if(compBasisAddressIn(i) > compBasisAddressIn(i+1)){

                int tempVal = compBasisAddressIn(i);
                compBasisAddressIn(i) = compBasisAddressIn(i+1);
                compBasisAddressIn(i+1) = tempVal;

                IdealOpConverted.col(i).swap(IdealOpConverted.col(i+1));

                break;

            }

            if(i==compSubspaceDim-2) return;

        }

    }

}

void MeritFunction::setCompBasisAddress(Eigen::MatrixXi& basisVector){

    compBasisAddressIn.resize(compSubspaceDim);

    for(int i=0;i<compSubspaceDim;i++){

        for(int j=0;j<g(photons,modes);j++){

            if(isEqual(computationalBasisIn.row(i),basisVector.row(j))){

                compBasisAddressIn(i) = j;

            }

        }

    }

    compBasisAddressOut.resize(compSubspaceDim);

    for(int i=0;i<compSubspaceDim;i++){

        for(int j=0;j<g(photons,modes);j++){

            if(isEqual(computationalBasisOut.row(i),basisVector.row(j))){

                compBasisAddressOut(i) = j;

            }

        }

    }

    return;

}



MeritFunction::MeritFunction(){

}


int MeritFunction::g(int n,int m){
    if(n==0 && m==0){
        return 0;
    }
    else if(n==0 && m>0){
        return 1;
    }

    else{
        return (int)(doublefactorial(n+m-1)/(doublefactorial(n)*doublefactorial(m-1))+0.5);
    }
}


double MeritFunction::doublefactorial(int x){
    double total=1.0;
    if (x>=0){
        for(int i=x;i>0;i--){
            total=i*total;
        }
    }
    else{
        std::cout << "invalid factorial" << std::endl;
        total=-1;
    }
    return total;
}


Eigen::MatrixXi MeritFunction::generateSubBasisVector(int subPhotons, int subModes){

    int markers = subPhotons + subModes - 1;
    int myints[markers];
    int i = 0;
    while(i<subPhotons){
        myints[i]=1;
        i++;
    }
    while(i<markers){
        myints[i]=0;
        i++;
    }
    Eigen::MatrixXi nv = Eigen::MatrixXi::Zero(g(subPhotons,subModes),subModes);
    i = 0;
    int j,k = 0;
    do {
        j = 0;
        k = 0;
        while(k<markers){
        if(myints[k]==1){
            nv(i,j)=nv(i,j)+1;
        }
        else if(myints[k]==0){
            j++;
        }

        k++;
        }
        i++;
    } while ( std::prev_permutation(myints,myints+markers) );

    return nv;

}



Eigen::MatrixXi MeritFunction::generateBasisVector(int subPhotons,int subModes, int subMeasureModes){
    Eigen::MatrixXi output(0,subModes);
    for(int i=subPhotons;i>=0;i--){

        Eigen::MatrixXi measModes = generateSubBasisVector(i,subMeasureModes);
        if(subMeasureModes==subModes){
            return measModes;
        }

        Eigen::MatrixXi othaModes = generateSubBasisVector(subPhotons-i,subModes-subMeasureModes);

        int numbRows = measModes.rows()*othaModes.rows();
        int numbMeasRows = measModes.rows();
        int numbOthaRows = othaModes.rows();
        int outputrows = output.rows();
        output.conservativeResize(outputrows+numbRows,subModes);

        for(int k=0;k<numbMeasRows;k++){
            for(int j=0;j<numbOthaRows;j++){
                for(int l=0;l<subMeasureModes;l++){
                    output(outputrows+j+numbOthaRows*k,l) = measModes(k,l);
                }
                for(int l=subMeasureModes;l<subModes;l++){
                    output(outputrows+j+numbOthaRows*k,l) = othaModes(j,l-subMeasureModes);
                }
            }
        }
    }

    return output;

}

bool MeritFunction::isEqual(Eigen::ArrayXi a1,Eigen::ArrayXi a2){

    assert(a1.size() == a2.size() && "Error: arrays of different sizes are being compared.");

    for(int i=0;i<a1.size();i++){

        if(a1(i) != a2(i)) return false;

    }

    return true;

}
