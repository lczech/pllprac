pllprac
=======

**NAME**  
    pllprac - Programming Practical of the Exelicis Lab
    using the Phylogenetic Likelihood Library (PLL).

**SYNOPSIS**  
    `pllprac [-e|m] [-v] phylip_file`  

**DESCRIPTION**  
    This program takes a phylip multiple sequence alignment file
    as input and evaluates all possible 203 time-reversible
    substitution matrices on this data.
    It then takes the best matrix as the substitution model
    to run a full ML tree search for the taxa in the file.
    The resulting tree is outputted to a file.

**OPTIONS**  
    `-e|m`  
        Use (E)mpirical base frequencies or (M)L estimates.  

    `-v`  
        Verbose user output.  

**ENVIRONMENT**  
    * PLL  
        In order to run this program, you need an installation
        of the Phylogenetic Likelihood Library (PLL), which can be
        obtained from `http://www.libpll.org/`  
    * MPI  
        The program can be compiled with MPI support via `make mpi`.  

**AUTHOR**  
    Lucas Czech <lucas.czech@h-its.org>
