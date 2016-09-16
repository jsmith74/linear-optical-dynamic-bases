#include "ChaseAlgorithm.h"


ChaseAlgorithm::ChaseAlgorithm(int n,int k){

    totalSet.resize(n);

    subset.resize(k);

    for(int i=0;i<n;i++) totalSet[i] = i;

    for(int i=0;i<k;i++) subset[i] = i;

    subsetDone = 0;

    subsetPermutations = factorial(k);

    subsetSize = k;

    totalSize = n;

}

void ChaseAlgorithm::iterate(){

    if(subsetDone < subsetPermutations - 1){

        iterateSubset();

    }

    else{

        iterateSubset();

        iterateSubsetElements();

       subsetDone = 0;

    }

    return;

}

void ChaseAlgorithm::iterateSubsetElements(){

    subset[subsetSize-1]++;

    int minAdjustedIndex;

    bool adjusted = false;

    int j=0;
    for(int i=subsetSize-1;i>0;i--){

        if(subset.at(i) >= totalSize - j){

            subset.at(i-1)++;
            minAdjustedIndex = i;
            adjusted = true;

        }

        j++;

    }

    if(adjusted){

        for(int i=minAdjustedIndex;i<subsetSize;i++){

            subset.at(i) = subset.at(i-1) + 1;

        }

    }

    return;

}

void ChaseAlgorithm::iterateSubset(){

        const int asize = subset.size();
        int a[asize];
        for(int i=0;i<asize;i++) a[i] = subset[i];
        std::next_permutation(a,a+asize);
        for(int i=0;i<asize;i++) subset[i] = a[i];
        subsetDone++;
        return;

}


int ChaseAlgorithm::factorial(int x){

    double total=1.0;

    assert(x<13 && "ERROR: INTEGER OVERFLOW IN FACTORIAL FUNCTION");

    if (x>=0){
        for(int i=x;i>0;i--){
            total=i*total;
        }
    }
    else{
        std::cout << "invalid factorial" << std::endl;
        total=-1;
    }
    return (int) (total+0.5);
}
