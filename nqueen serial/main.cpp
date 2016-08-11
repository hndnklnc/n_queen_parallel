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
    Genetic gen(size);
    gen.initializePopulation();


    return 0;
}
