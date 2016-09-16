#ifndef CHASEALGORITHM_H_INCLUDED
#define CHASEALGORITHM_H_INCLUDED

#include <vector>
#include <iostream>
#include <algorithm>
#include <assert.h>

class ChaseAlgorithm{

    public:

        ChaseAlgorithm(int n,int k);
        std::vector<int> totalSet,subset;
        void iterate();

    private:

        int subsetDone,subsetSize,totalSize;
        int factorial(int x);
        int subsetPermutations;
        void iterateSubset();
        void iterateSubsetElements();

};


#endif // CHASEALGORITHM_H_INCLUDED
