#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <assert.h>
#include <unistd.h>

#include "pll/pll.h"
#include "main.h"

using namespace std;

/**
 * Calculate the likelihood given a particular substitution model.
 *
 * @param   alignmentFilename
 * @param   parts
 * @param   model
 * @param   bf
 * @param   likelihood
 * @param   substRates
 */
void analyzeModel (
    string alignmentFilename, pllQueue* parts, string model, BaseFrequencyType bf,
    double* likelihood, double** substRates
    )
{
    //model = "0,0,1,2,0,0";

    // parse the alignment file (PLL cannot reuse the alignment data...)
    pllAlignmentData* alignmentData = pllParseAlignmentFile (PLL_FORMAT_PHYLIP, alignmentFilename.c_str());

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
    *likelihood = inst->likelihood;
    *substRates = partitions->partitionData[0]->substRates;

    // destroy the PLL instance
    pllAlignmentDataDestroy (alignmentData);
    pllPartitionsDestroy    (inst, &partitions); // again, expects a **partitionList here
    pllDestroyInstance      (inst);
}

/**
 * Prints a 4x4 matrix of NT substitution rates
 *
 * @param sr    Subst rate matrix of type double[6]
 */
void printSubstMatrix (double* sr)
{
    cout << fixed << setprecision(3) << setw(12) << -(sr[0]+sr[1]+sr[2]);
    cout << fixed << setprecision(3) << setw(12) << sr[0];
    cout << fixed << setprecision(3) << setw(12) << sr[1];
    cout << fixed << setprecision(3) << setw(12) << sr[2];
    cout << endl;

    cout << fixed << setprecision(3) << setw(12) << sr[0];
    cout << fixed << setprecision(3) << setw(12) << -(sr[0]+sr[3]+sr[4]);
    cout << fixed << setprecision(3) << setw(12) << sr[3];
    cout << fixed << setprecision(3) << setw(12) << sr[4];
    cout << endl;

    cout << fixed << setprecision(3) << setw(12) << sr[1];
    cout << fixed << setprecision(3) << setw(12) << sr[3];
    cout << fixed << setprecision(3) << setw(12) << -(sr[1]+sr[3]+sr[5]);
    cout << fixed << setprecision(3) << setw(12) << sr[5];
    cout << endl;

    cout << fixed << setprecision(3) << setw(12) << sr[2];
    cout << fixed << setprecision(3) << setw(12) << sr[4];
    cout << fixed << setprecision(3) << setw(12) << sr[5];
    cout << fixed << setprecision(3) << setw(12) << -(sr[2]+sr[4]+sr[5]);
    cout << endl;
}

void evaluateTree()
{
    //pllTreeEvaluate ( pllInstance *tr, partitionList *pr, int maxSmoothIterations );
}

int main (int argc, char* argv[])
{
    BaseFrequencyType   bf      = ML;       // option for using empirical vs ML base freqs
    bool                verbose = false;    // for user output

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
    partitionFile.close();

    // create a partition queue from the alignment
    int sequenceLength = alignmentData->sequenceLength;
    pllQueue* parts = pllPartitionParse (partitionFilename.c_str());
    if (pllPartitionsValidate (parts, alignmentData) == PLL_FALSE) {
        cout << "Something is wrong with the partitions, aborting." << endl;
        exit (0);
    }
    pllAlignmentDataDestroy (alignmentData);

    // prepare values for main loop
    double  likelihood;
    double* substRates = NULL;
    double  bestL = 1.0;    // best likelihood
    string  bestM;          // best model
    double  bestS[6];          // best substitution matrix

    // user output
    cout << "Analysing all " << MODEL_STRINGS.size() << " substitution models..." << endl;
    if (verbose) {
        cout << endl << "Model         Likelihood          AIC          BIC" << endl;
        cout         << "----------- ------------ ------------ ------------" << endl;
    }

    // try out every substitution model to find the best one
    // (PLL cannot reuse the data, so we have to do all the input file reading
    // anew for every substitution model)
    for (string model : MODEL_STRINGS) {
        // first, I broke down the process into smaller steps, but the objects
        // are so closely interlinked and the init procedures appear all over
        // the code, so that is was even a greater mess then this long func
        analyzeModel(alignmentFilename, parts, model, bf, &likelihood, &substRates);

        // user output
        if (verbose) {
            cout << model << " ";
            cout << fixed << setprecision(3) << setw(12) << likelihood << " ";
            cout << fixed << setprecision(3) << setw(12) << AIC(likelihood, 6) << " ";
            cout << fixed << setprecision(3) << setw(12) << BIC(likelihood, 6, sequenceLength) << endl;
        }

        // store best model (AIC and BIC are monotonous over the likelihood
        // for fixed input and parameter count, so we do not need them here)
        if (bestL > 0.0  ||  likelihood > bestL) {
            bestL = likelihood;
            bestM = model;
            for (int i = 0; i < 6; i++) {
                bestS[i] = substRates[i];
            }
        }
    }

    // user output
    cout << "...done." << endl << endl;
    cout << "Best model: " << bestM << endl;
    cout << "Likelihood: " << bestL << endl;
    cout << "AIC:        " << AIC(bestL, 6) << endl;
    cout << "BIC:        " << BIC(bestL, 6, sequenceLength) << endl;
    cout << endl << "Substitution matrix:" << endl;
    printSubstMatrix(bestS);

    pllQueuePartitionsDestroy (&parts); // for some reason, PLL expects **pQueue here...
    return 0;
}
