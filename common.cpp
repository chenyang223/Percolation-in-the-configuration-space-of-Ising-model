#include "head.h"

double RandDouble(int i)
{
    return rand() / (double)(RAND_MAX)*i;
}

int RandInt(int i)
{
    return (rand() % i);
}
