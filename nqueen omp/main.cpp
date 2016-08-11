#include <iostream>
#include "Genetic.h"
#include <time.h>
using namespace std;

static int SIZE;
int main()
{
    int size, thread_num, gra;
    cout<<"Enter size of the queen problem: ";
    cin>>size;
    cout<<"\nEnter thread num: ";
    cin>>thread_num;
    cout<<"\nEnter  granularity: ";
    cin>>gra;
    Genetic gen(size, thread_num, gra);
    gen.initializePopulation();


    return 0;
}

