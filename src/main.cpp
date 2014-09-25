#include <iostream>
#include <assert.h>
#include "pll/pll.h"

int main (int argc, char** argv)
{
    std::cout << "Hello World!\n";

    pllInstance * inst;
    pllInstanceAttr attr;
    pllAlignmentData * alignmentData;

    /* set PLL instance attributes */
    attr.rateHetModel = PLL_GAMMA;
    attr.fastScaling  = PLL_FALSE;
    attr.saveMemory   = PLL_FALSE;
    attr.useRecom     = PLL_FALSE;
    attr.randomNumberSeed = 0x12345;

    /* Create the PLL instance */
    inst = pllCreateInstance (&attr);

    /* Destroy the PLL instance */
    //pllDestroyInstance (inst);

    std::cout << "Goodbye.\n";
    return 0;
}
