#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <assert.h>
#include <unistd.h>

#include "pll/pll.h"
#include "main.h"

using namespace std;
bool verbose = false; // for user output

double analyzeModel(string alignmentFilename, string partitionFilename, string model, BaseFrequencyType bf)
{
    //model = "0,0,1,2,0,0";

    // user output
    if (verbose) {
        cout << endl << "=============================" << endl;
        cout << "     Model " << model << endl << endl;
    }

    // parse the input files
    pllAlignmentData* alignmentData = pllParseAlignmentFile (PLL_FORMAT_PHYLIP, alignmentFilename.c_str());
    pllQueue* parts = pllPartitionParse (partitionFilename.c_str());
    if (pllPartitionsValidate (parts, alignmentData) == PLL_FALSE) {
        cout << "Something is wrong with the partitions, aborting." << endl;
        exit (0);
    }

    // create partitions object, optimizing everything
    partitionList* partitions = pllPartitionsCommit (parts, alignmentData);
    pllAlignmentRemoveDups (alignmentData, partitions);
    partitions->partitionData[0]->optimizeBaseFrequencies   = true;
    partitions->partitionData[0]->optimizeAlphaParameter    = true;
    partitions->partitionData[0]->optimizeSubstitutionRates = true;

    // set PLL instance attributes
    pllInstanceAttr attr;
    attr.rateHetModel       = PLL_GAMMA;
    attr.fastScaling        = PLL_FALSE;
    attr.saveMemory         = PLL_FALSE;
    attr.useRecom           = PLL_FALSE;
    attr.randomNumberSeed   = 0x12345;

    // create the actual PLL instance using the alignment
    pllInstance* inst = pllCreateInstance (&attr);
    pllTreeInitTopologyForAlignment (inst, alignmentData);
    pllLoadAlignment (inst, alignmentData, partitions, PLL_DEEP_COPY);

    // create a tree
    pllComputeRandomizedStepwiseAdditionParsimonyTree (inst, partitions);
    pllInitModel (inst, partitions, alignmentData);

    // if user wishes, we take empirical base freqs instead of ML estimates
    if (bf == EMPIRICAL) {
        pllSetFixedBaseFrequencies(partitions->partitionData[0]->empiricalFrequencies, 4, 0, partitions, inst);
        // (this also automatically resets optimizeBaseFrequencies)
    }

    // use the substutition model to optimize the model params
    if (!pllSetSubstitutionRateMatrixSymmetries (const_cast<char*>(model.c_str()), partitions, 0)) {
        cout << "Cannot set substitution rate matrix to " << model << ", aborting.";
        exit (0);
    }
    initReversibleGTR (inst, partitions, 0);
    pllOptimizeModelParameters (inst, partitions, 0.1);
    double likelihood = inst->likelihood;

    if (verbose) {
        cout << "Substitution rates: " << endl;
        for (int i = 0; i < 6; i++) {
            cout << "Rate " << i << ": " << partitions->partitionData[0]->substRates[i] << endl;
        }
        cout << endl;
    }

    // destroy the PLL instance
    pllAlignmentDataDestroy     (alignmentData);
    pllQueuePartitionsDestroy   (&parts); // for some reason, PLL expects **pQueue here...
    pllPartitionsDestroy        (inst, &partitions); // again, expects a **partitionList here
    pllDestroyInstance          (inst);

    return likelihood;
}

int main (int argc, char* argv[])
{
    // option for using empirical vs ML base freqs
    BaseFrequencyType bf = ML;

    // parse command line arguments
    int c;
    while ((c = getopt (argc, argv, "e::m::v::")) != -1) {
        switch (c) {
            case 'e':
                bf = EMPIRICAL;
                break;

            case 'm':
                bf = ML;
                break;

            case 'v':
                verbose = true;
                break;

            case '?':
                if (isprint (optopt))
                    fprintf (stderr, "Unknown option `-%c'.\n", optopt);
                else
                    fprintf (stderr, "Unknown option character `\\x%x'.\n", optopt);
            default:
                exit (0);
        }
    }

    // check if PHYLIP alignment file is given in arguments
    if (optind >= argc) {
        cout << "No PHYLIP input file given, aborting." << endl;
        exit (0);
    }

    // some user output
    string alignmentFilename = argv[optind];
    cout << "Using PHYLIP input file '" << alignmentFilename;
    switch (bf) {
        case EMPIRICAL:
            cout << " with empirical base frequencies." << endl;
            break;

        case ML:
            cout << " with ML estimated base frequencies." << endl;
            break;
    }

    // create a partition file containing only one partition
    // (PLL does not provide a default "use all" partition and
    // cannot take the partition via string parameter)
    pllAlignmentData* alignmentData = pllParseAlignmentFile (PLL_FORMAT_PHYLIP, alignmentFilename.c_str());
    string partitionFilename = alignmentFilename + ".partition";
    ofstream partitionFile (partitionFilename.c_str());
    if (!partitionFile.is_open()) {
        cout << "Unable to write partition file '" << partitionFilename << "', aborting." << endl;
        exit (0);
    }
    cout << "Writing partition file '" << partitionFilename << "'." << endl;
    partitionFile << "DNA, Partition = 1-" << alignmentData->sequenceLength << "\n";
    int sequenceLength = alignmentData->sequenceLength;
    partitionFile.close();
    pllAlignmentDataDestroy (alignmentData);

    // prepare values for main loop
    double likelihood;

    // try out every substitution model to find the best one
    // (PLL cannot reuse the data, so we have to do all the input file reading
    // anew for every substitution model)
    for (string model : MODEL_STRINGS) {
        // first, I broke down the process into smaller steps, but the objects
        // are so closely interlinked and the init procedures appear all over
        // the code, so that is was even a greater mess then this long func
        likelihood = analyzeModel(alignmentFilename, partitionFilename, model, bf);

        if (verbose) {
            cout << "Likelihood " << likelihood << endl;
            cout << "AIC        " << AIC(likelihood, 6) << endl;
            cout << "BIC        " << BIC(likelihood, 6, sequenceLength) << endl;
        }
    }

    return 0;
}
