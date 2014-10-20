// necessary includes
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cstring>
#include <vector>
#include <algorithm>
#include <assert.h>
#include <unistd.h>

// optional include defined by makefile
#ifdef USE_MPI
    #include "mpi.h"
#endif

// our own headers
#include "pll/pll.h"
#include "main.h"

// define macros to make code more readable
#ifdef USE_MPI
    #define IF_MAIN_THREAD if (MPI::COMM_WORLD.Get_rank() == 0)
#else
    #define IF_MAIN_THREAD if (true)
#endif

using namespace std;

#ifdef USE_MPI
MPI_Datatype ModelScoreMPI;

void initMpiTypes ()
{
    int blocks[2]         = {9, 12};
    MPI_Datatype types[2] = {MPI_DOUBLE, MPI_CHAR};

    MPI_Aint displacements[2];
    MPI_Aint dblex;

    MPI_Type_extent(MPI_DOUBLE, &dblex);

    displacements[0] = static_cast<MPI_Aint>(0);
    displacements[1] = dblex*9;

    MPI_Type_struct(2, blocks, displacements, types, &ModelScoreMPI);
    MPI_Type_commit (&ModelScoreMPI);
}
#endif

/*
 * Returns the number of free parameters for a particular model.
 */
int countFreeParameters (string model, int taxaCount)
{
    // parameters:
    // frequencies -> 3
    // alpha param -> 1
    // num of taxa -> 2 * taxaCount - 3
    // subst model -> in interval [0-5] (added at the end of this function)
    int r = 3 + 1 + 2 * taxaCount - 3;

    // sanity check of the subst model
    if (model[0] != '0' || model.length() != 11) {
        cout << "Wrong  model string '" << model << "', ignoring those parameters." << endl;
        return r;
    }

    // count how many parameters the subst model has.
    // (assuming a model in the form of '0,1,2,3,2,0',
    // this is just the highest number appearing in the string)
    int h = 0;
    char* p = strtok(const_cast<char*>(model.c_str()), ",");
    while (p != NULL) {
        h = max(h, atoi(p));
        p = strtok(NULL, ",");
    }

    return r+h;
}

/*
 * Prints a 4x4 matrix of NT substitution rates.
 *
 * @param sr    Subst rate matrix of type double[6]
 */
void printSubstMatrix (double* sr, ostream &os)
{
    auto pf = [](double v, ostream &os) {
        os << fixed << setprecision(3) << setw(12) << v;
    };

    pf (-(sr[0]+sr[1]+sr[2]), os);
    pf (sr[0], os);
    pf (sr[1], os);
    pf (sr[2], os);
    os << endl;

    pf (sr[0], os);
    pf (-(sr[0]+sr[3]+sr[4]), os);
    pf (sr[3], os);
    pf (sr[4], os);
    os << endl;

    pf (sr[1], os);
    pf (sr[3], os);
    pf (-(sr[1]+sr[3]+sr[5]), os);
    pf (sr[5], os);
    os << endl;

    pf (sr[2], os);
    pf (sr[4], os);
    pf (sr[5], os);
    pf (-(sr[2]+sr[4]+sr[5]), os);
    os << endl;
}

/*
 * Creates a PLL instance given an alignment and a partition file.
 */
void createPLL (
    const string alignmentFilename, pllQueue* parts,
    pllAlignmentData* &alignmentData, partitionList* &partitions, pllInstance* &inst
    )
{
    // parse the alignment file (PLL cannot reuse the alignment data...)
    alignmentData = pllParseAlignmentFile (PLL_FORMAT_PHYLIP, alignmentFilename.c_str());

    // create partitions object, optimizing everything
    partitions = pllPartitionsCommit (parts, alignmentData);
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
    inst = pllCreateInstance (&attr);
    pllTreeInitTopologyForAlignment (inst, alignmentData);
    pllLoadAlignment (inst, alignmentData, partitions);
}

/*
 * Destroys the instance and other objects that createPLL yielded.
 */
void destroyPLL (pllAlignmentData* alignmentData, partitionList* partitions, pllInstance* inst)
{
    // destroy everything created by createPLL
    pllAlignmentDataDestroy (alignmentData);
    pllPartitionsDestroy    (inst, &partitions);
    pllDestroyInstance      (inst);
}

/*
 * Calculate the likelihood given a particular substitution model.
 *
 * @param   alignmentFilename
 * @param   parts
 * @param   model
 * @param   bf
 * @param   doSearch            If set, runs a ML evaluation to find the best tree
 * @param   likelihood
 * @param   substRates
 * @param   treeString
 */
void analyzeModel (
    const string alignmentFilename, pllQueue* parts, const string model, const BaseFrequencyType bf, const bool doSearch,
    double &likelihood, double* &substRates, string &treeString
    )
{
    // create a PLL instance
    pllAlignmentData*   alignmentData;
    partitionList*      partitions;
    pllInstance*        inst;
    createPLL (alignmentFilename, parts, alignmentData, partitions, inst);

    // create a tree
    //pllTreeInitTopologyRandom (inst, 10, alignmentData->sequenceLabels);
    pllComputeRandomizedStepwiseAdditionParsimonyTree (inst, partitions);
    pllInitModel (inst, partitions);

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
    pllOptimizeModelParameters (inst, partitions, 0.1);

    if (doSearch) {
        // find the best tree
        pllRaxmlSearchAlgorithm (inst, partitions, true);
        pllTreeToNewick (inst->tree_string, inst, partitions, inst->start->back, true, true, false, false, false, PLL_SUMMARIZE_LH, false, false);
        treeString = inst->tree_string;
    }

    // hand over results to referenced vars
    likelihood = inst->likelihood;
    for (int i = 0; i < 6; i++) {
        substRates[i] = partitions->partitionData[0]->substRates[i];
    }

    // destroy the PLL instance
    destroyPLL (alignmentData, partitions, inst);
}

/*
 *
 */
int main (int argc, char* argv[])
{
    // command line dependend options
    BaseFrequencyType   bf      = ML;       // option for using empirical vs ML base freqs
    bool                verbose = false;    // for more user output

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

    // init mpi
#ifdef USE_MPI
    MPI::Init (argc, argv);
    initMpiTypes();
#endif

    // check if PHYLIP alignment file is given in arguments
    if (optind >= argc) {
        cout << "No PHYLIP input file given, aborting." << endl;
        exit (0);
    }
    string alignmentFilename = argv[optind];

    IF_MAIN_THREAD {
        // user output
        cout << "Using PHYLIP input file '" << alignmentFilename;
        switch (bf) {
            case EMPIRICAL:
                cout << " with empirical base frequencies." << endl;
                break;

            case ML:
                cout << " with ML estimated base frequencies." << endl;
                break;
        }
    }

    // create a partition file containing only one partition
    // (PLL does not provide a default "use all" partition and
    // cannot take the partition via string parameter)
    pllAlignmentData* alignmentData = pllParseAlignmentFile (PLL_FORMAT_PHYLIP, alignmentFilename.c_str());
    string partitionFilename = alignmentFilename + ".partition";
    IF_MAIN_THREAD {
        ofstream partitionFile (partitionFilename.c_str());
        if (!partitionFile.is_open()) {
            cout << "Unable to write partition file '" << partitionFilename << "', aborting." << endl;
            exit (0);
        }
        cout << "Writing partition file '" << partitionFilename << "'." << endl;
        partitionFile << "DNA, Partition = 1-" << alignmentData->sequenceLength << "\n";
        partitionFile.close();
    }

    // make sure the partition file is written before other threads continue
#ifdef USE_MPI
    MPI_Barrier (MPI_COMM_WORLD);
#endif

    // create a partition queue from the alignment
    int sequenceLength = alignmentData->sequenceLength;
    int sequenceCount  = alignmentData->sequenceCount;
    pllQueue* parts = pllPartitionParse (partitionFilename.c_str());
    if (pllPartitionsValidate (parts, alignmentData) == PLL_FALSE) {
        cout << "Something is wrong with the partitions, aborting." << endl;
        exit (0);
    }
    pllAlignmentDataDestroy (alignmentData);

    // prepare values for main loop
    int     paramCount;
    double  likelihood, aic, bic;
    double* substRates;
    string  treeString;
    vector<std::string> model_strings;

    ModelScore best;
    best.likelihood =  1.0;
    best.AIC        = -1.0;
    best.BIC        = -1.0;

#ifdef USE_MPI
    // PLL does some output when reading the sequence above, and we don't
    // want that ouput to appear in between the output of the main loop,
    // so we wait until all threads are done with this
    MPI_Barrier (MPI_COMM_WORLD);

    // split up the models for mpi threads
    for (unsigned int i = 0; i < MODELS.size(); i++) {
        if ((i - MPI::COMM_WORLD.Get_rank()) % MPI::COMM_WORLD.Get_size() == 0) {
            model_strings.push_back (MODELS[i]);
        }
    }

    // user output
    if (MPI::COMM_WORLD.Get_rank() == 0) {
        cout << endl << "Analysing all " << MODELS.size() << " substitution models..." << endl;
        if (verbose) {
            cout << endl;
            cout << "MPI Model         Likelihood          AIC          BIC" << endl;
            cout << "--- ----------- ------------ ------------ ------------" << endl;
        }
    }

    // the threads shall not print output before the table header is printed,
    // so we wait here for the main thread
    MPI_Barrier (MPI_COMM_WORLD);
#else
    // do all models
    model_strings = MODELS;

    // user output
    cout << endl << "Analysing all " << MODELS.size() << " substitution models..." << endl;
    if (verbose) {
        cout << endl;
        cout << "Model         Likelihood          AIC          BIC" << endl;
        cout << "----------- ------------ ------------ ------------" << endl;
    }
#endif

    // try out every substitution model to find the best one
    // (PLL cannot reuse the data, so we have to do the input file reading
    // anew for every substitution model)
    for (string model : model_strings) {
        // get the likelihood for the current model
        substRates = new double[6];
        analyzeModel(alignmentFilename, parts, model, bf, false, likelihood, substRates, treeString);

        // calculate our criteria
        paramCount = countFreeParameters(model, sequenceCount);
        aic = AIC(likelihood, paramCount);
        bic = BIC(likelihood, paramCount, sequenceLength);

        // user output
        if (verbose) {
#ifdef USE_MPI
            cout << fixed << setw(3) << MPI::COMM_WORLD.Get_rank() << " ";
#endif
            cout << model << " ";
            cout << fixed << setprecision(3) << setw(12) << likelihood << " ";
            cout << fixed << setprecision(3) << setw(12) << aic << " ";
            cout << fixed << setprecision(3) << setw(12) << bic << endl;
        }

        // if AIC and BIC are better for the current model than for a previous one,
        // store this model
        if (best.AIC < 0.0 || best.BIC < 0.0 || ( aic < best.AIC && bic < best.BIC) ) {
            best.likelihood = likelihood;
            best.AIC        = aic;
            best.BIC        = bic;
            strncpy(best.model, model.c_str(), 12);
            for (int i = 0; i < 6; i++) {
                best.rates[i] = substRates[i];
            }
        }

        // clean up
        delete [] substRates;
    }

#ifdef USE_MPI
    // main thread reduces results to actual best
    if (MPI::COMM_WORLD.Get_rank() == 0) {
        ModelScore nbest;
        MPI_Status status;

        for (int i = 1; i < MPI::COMM_WORLD.Get_size(); i++) {
            MPI_Recv(&nbest, 1, ModelScoreMPI, i, 0, MPI_COMM_WORLD, &status);

            if (nbest.AIC < best.AIC && nbest.BIC < best.BIC) {
                best.likelihood = nbest.likelihood;
                best.AIC        = nbest.AIC;
                best.BIC        = nbest.BIC;
                strncpy(best.model, nbest.model, 12);
                for (int i = 0; i < 6; i++) {
                    best.rates[i] = nbest.rates[i];
                }
            }
        }
    } else {
        // send results to main thread
        MPI_Send(&best, 1, ModelScoreMPI, 0, 0, MPI_COMM_WORLD);
    }
#endif

    IF_MAIN_THREAD {
        // user output
        cout << "...done." << endl << endl;
        cout << "Best model: " << best.model << endl;
        cout << "Likelihood: " << best.likelihood << endl;
        cout << "AIC:        " << best.AIC << endl;
        cout << "BIC:        " << best.BIC << endl;
        cout << endl << "Substitution matrix:" << endl;
        printSubstMatrix(best.rates, cout);

        // write subst matrix to file
        string matrixFilename = alignmentFilename + ".matrix";
        ofstream matrixFile (matrixFilename.c_str());
        if (!matrixFile.is_open()) {
            cout << endl << "Unable to write matrix file '" << matrixFilename << "'." << endl;
        } else {
            cout << endl << "Writing matrix file '" << matrixFilename << "'." << endl;
            printSubstMatrix(best.rates, matrixFile);
            matrixFile.close();
        }

        // run the ML tree search
        cout << endl << "Running ML tree search with best model..." << endl;
        analyzeModel(alignmentFilename, parts, best.model, bf, true, likelihood, substRates, treeString);
        cout << "...done." << endl;

        // write tree to file
        string treeFilename = alignmentFilename + ".tree";
        ofstream treeFile (treeFilename.c_str());
        if (!treeFile.is_open()) {
            cout << "Unable to write tree file '" << treeFilename << "'." << endl;
        } else {
            cout << endl << "Writing tree file '" << treeFilename << "'." << endl;
            treeFile << treeString;
            treeFile.close();
        }
    }

    // clean up
    pllQueuePartitionsDestroy (&parts);
#ifdef USE_MPI
    MPI_Type_free (&ModelScoreMPI);
    MPI_Finalize();
#endif
    return 0;
}
