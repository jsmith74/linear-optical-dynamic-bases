#include "BFGS_Optimization.h"
#include "ChaseAlgorithm.h"
#include <fstream>

int main(){

    remove("Successful Basis Change.dat");

    remove("outbasisCheck.dat");

    ChaseAlgorithm NChooseKStates(10,4);

    for(int i=0;i<5040;i++){

        BFGS_Optimization optimizer(0.00001,2.0,0.0,NChooseKStates.subset);

        for(int j=0;j<40;j++){

            if(optimizer.meritFunction.validBasis == true) continue;

            optimizer.minimize();


        }

        std::ofstream outfile("outbasisCheck.dat",std::ofstream::app);
        for(int j=0;j<4;j++) outfile << NChooseKStates.subset.at(j) << " ";
        outfile << std::endl;
        outfile.close();

        NChooseKStates.iterate();

    }


    return 0;

}
