#include <iostream>
#include "Geneticmpi.h"
#include <time.h>
using namespace std;

static int SIZE;
int main()
{

    int size, thread_num;
    cout<<"Enter size of the queen problem: ";
    cin>>size;
    cout<<"\nthread size: ";
    cin>>thread_num;
    Genetic gen(size, thread_num);
    gen.initializePopulation();


    return 0;
}
