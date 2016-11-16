#include "BFGS_Optimization.h"
#include "ChaseAlgorithm.h"
#include <fstream>
#include <stdlib.h>
#include <string>
#include <sstream>


int main(){

    std::cout << "RAND_MAX check: "  << RAND_MAX << std::endl << std::endl;

    srand(611*time(NULL));

    while(true){

        int randomInBasisChoice = rand() % 32760 + 1;

        int randomOutBasisChoice = rand() % 32760 + 1;

        ChaseAlgorithm NChooseKStatesIn(15,4);

        ChaseAlgorithm NChooseKStatesOut(15,4);

        for(int r=0;r<randomInBasisChoice;r++){

            NChooseKStatesIn.iterate();

        }

        for(int r=0;r<randomOutBasisChoice;r++){

            NChooseKStatesOut.iterate();

        }


        BFGS_Optimization optimizer(0.00001,2.0,0.0,NChooseKStatesIn.subset,NChooseKStatesOut.subset);

        for(int j=0;j<40;j++){

            if(optimizer.meritFunction.validBasis == true) continue;

            optimizer.minimize();

        }

    }

    return 1;

}
