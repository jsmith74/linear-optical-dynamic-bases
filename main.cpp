#include "BFGS_Optimization.h"
#include "ChaseAlgorithm.h"
#include <fstream>
#include <stdlib.h>
#include <string>
#include <sstream>

std::string command;

void setInbasisCommand(int& lines){

    std::stringstream ss;
    ss << lines + 2;
    ss >> command;

    command = "tail InBasis.dat -n +" + command;

    command = command + " > InBasisNew.dat && mv InBasisNew.dat InBasis.dat";

    return;

}

int main(){

    remove("Successful Basis Change.dat");

    remove("Successful Outbasis List.dat");

    remove("BasisCheck.dat");

    int compBasisDim = 4;

    setInbasisCommand(compBasisDim);

    while(true){

        ChaseAlgorithm NChooseKStates(10,4);

        for(int i=0;i<5040;i++){

            BFGS_Optimization optimizer(0.00001,2.0,0.0,NChooseKStates.subset);

            for(int j=0;j<40;j++){

                if(optimizer.meritFunction.validBasis == true) continue;

                optimizer.minimize();


            }

            NChooseKStates.iterate();

        }

        compBasisDim = system(command.c_str());

        std::ifstream eofTest("InBasis.dat");
        eofTest >> compBasisDim;
        if(eofTest.eof() == true){

            eofTest.close();
            return 0;

        }

        eofTest.close();

    }

    return 1;

}
