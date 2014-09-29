#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <assert.h>
#include "pll/pll.h"
#include "main.h"

using namespace std;

void analyzeModel(string alignmentFilename, string partitionFilename, string model)
{
    //string model = "1,1,1,1,1,2";
    cout << endl << "=============================" << endl;
    cout << " Model " << model << endl;

    // parse the alignment file
    pllAlignmentData* alignmentData = pllParseAlignmentFile (PLL_FORMAT_PHYLIP, alignmentFilename.c_str());

    // parse the partiotion file and connect it to the alignment
    pllQueue* parts = pllPartitionParse (partitionFilename.c_str());
    if (pllPartitionsValidate (parts, alignmentData) == PLL_FALSE) {
        cout << "Something is wrong with the partitions, aborting." << endl;
        exit(0);
    }

    // set PLL instance attributes
    pllInstanceAttr attr;
    attr.rateHetModel = PLL_GAMMA;
    attr.fastScaling = PLL_FALSE;
    attr.saveMemory = PLL_FALSE;
    attr.useRecom = PLL_FALSE;
    attr.randomNumberSeed = 0x12345;
    pllInstance* inst;
    partitionList* partitions;
    partitions = pllPartitionsCommit (parts, alignmentData);
    pllAlignmentRemoveDups (alignmentData, partitions);

    // optimize everything
    partitions->partitionData[0]->optimizeBaseFrequencies = true;
    partitions->partitionData[0]->optimizeAlphaParameter = true;
    partitions->partitionData[0]->optimizeSubstitutionRates = true;

    // load the alignment
    inst = pllCreateInstance (&attr);
    pllTreeInitTopologyForAlignment (inst, alignmentData);
    pllLoadAlignment (inst, alignmentData, partitions, PLL_DEEP_COPY);

    // create a tree
    pllComputeRandomizedStepwiseAdditionParsimonyTree (inst, partitions);
    pllInitModel (inst, partitions, alignmentData);
    pllSetSubstitutionRateMatrixSymmetries (const_cast<char*>(model.c_str()), partitions, 0);
    pllOptimizeModelParameters (inst, partitions, 0.1);

    // Destroy the PLL instance
    pllAlignmentDataDestroy (alignmentData);
    pllQueuePartitionsDestroy (&parts); // for some reason, PLL expects **pQueue here...
    pllPartitionsDestroy (inst, &partitions); // again, expects a **partitionList here
    pllDestroyInstance (inst);
}

int main (int argc, char* argv[])
{
    // check if PHYLIP alignment file is given
    if (argc == 1) {
        cout << "No PHYLIP input file given, aborting." << endl;
        exit(0);
    }
    string alignmentFilename = argv[1];
    cout << "Using PHYLIP input file '" << alignmentFilename << "'." << endl;

    // create a partition file containing only one partition
    // (PLL does not provide a default "use all" partition and
    // cannot take the partition via string parameter)
    pllAlignmentData* alignmentData = pllParseAlignmentFile (PLL_FORMAT_PHYLIP, alignmentFilename.c_str());
    string partitionFilename = alignmentFilename + ".partition";
    ofstream partitionFile (partitionFilename.c_str());
    if (!partitionFile.is_open()) {
        cout << "Unable to write partition file '" << partitionFilename << "', aborting." << endl;
        exit(0);
    }
    cout << "Writing partition file '" << partitionFilename << "'." << endl;
    partitionFile << "DNA, Partition = 1-" << alignmentData->sequenceLength << "\n";
    partitionFile.close();
    pllAlignmentDataDestroy (alignmentData);

    // try out every substitution model to find the best one
    // (PLL cannot reuse the data, so we have to do all the input file reading
    // anew for every substitution model)
    for (string model : MODEL_STRINGS) {
        analyzeModel(alignmentFilename, partitionFilename, model);
    }
    return 0;
}
