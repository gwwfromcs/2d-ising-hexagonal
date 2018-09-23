#include "ran1.h"
#include <iostream>

int main()
{
    long int seed = 1345559;
    for(int i = 1; i<100; i++)
    {
        std::cout << ran1(&seed) << std::endl;
    }
    return 0;
}
