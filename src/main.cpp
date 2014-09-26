#include <string>
#include <iostream>
#include <fstream>
#include <assert.h>

#include "pll/pll.h"

using namespace std;

int main (int argc, char* argv[])
{
    // check if PHYLIP alignment file is given
    if (argc == 1) {
        cout << "No PHYLIP input file given, aborting." << endl;
        return 0;
    }

    // parse the alignment file
    cout << "Using PHYLIP input file '" << argv[1] << "'." << endl;
    string alignmentFilename = argv[1];
    pllAlignmentData* alignmentData = pllParseAlignmentFile (PLL_FORMAT_PHYLIP, alignmentFilename.c_str());

    // create a partition file containing only one partition
    // (PLL cannot take the partition via string parameter)
    string partitionFilename = alignmentFilename + ".partition";
    ofstream partitionFile (partitionFilename.c_str());
    if (!partitionFile.is_open()) {
        cout << "Unable to write partition file '" << partitionFilename << "', aborting." << endl;
        return 0;
    }
    cout << "Writing partition file '" << partitionFilename << "'." << endl;
    partitionFile << "DNA, Partition = 1-" << alignmentData->sequenceLength << "\n";
    partitionFile.close();

    // parse the partiotion file and connect it to the alignment
    pllQueue* parts = pllPartitionParse (partitionFilename.c_str());
    if (pllPartitionsValidate (parts, alignmentData) == PLL_FALSE) {
        cout << "Something is wrong with the partitions, aborting." << endl;
        return 0;
    }
    partitionList* partitions = pllPartitionsCommit (parts, alignmentData);
    pllAlignmentRemoveDups (alignmentData, partitions);



    // set PLL instance attributes
    pllInstanceAttr attr;
    attr.rateHetModel     = PLL_GAMMA;
    attr.fastScaling      = PLL_FALSE;
    attr.saveMemory       = PLL_FALSE;
    attr.useRecom         = PLL_FALSE;
    attr.randomNumberSeed = 0x12345;

    // create the PLL instance
    pllInstance* inst;
    inst = pllCreateInstance (&attr);

    // load the alignment
    pllTreeInitTopologyForAlignment (inst, alignmentData);
    pllLoadAlignment (inst, alignmentData, partitions, PLL_DEEP_COPY);

    cout << partitions << endl;
    cout << partitions->numberOfPartitions << endl;
    cout << partitions->partitionData[0]->lower << endl;
    cout << partitions->partitionData[0]->upper << endl;

    // create a tree
    pllComputeRandomizedStepwiseAdditionParsimonyTree (inst, partitions);

    partitions->partitionData[0]->optimizeBaseFrequencies   = true;
    partitions->partitionData[0]->optimizeAlphaParameter    = true;
    partitions->partitionData[0]->optimizeSubstitutionRates = true;

    pllOptimizeModelParameters (inst, partitions, 0.1);

    //pllSetSubstitutionRateMatrixSymmetries (char *string, partitionList * pr, int model);

    // Destroy the PLL instance
    pllAlignmentDataDestroy (alignmentData);
    pllQueuePartitionsDestroy (&parts); // for some reason, PLL expects **pQueue here...
    //pllDestroyInstance (inst);
    return 0;
}
