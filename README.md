## Introduction
A Haskell implementation of the problem described
in Aruoba-Villaverde's paper "A Comparison of Programming
Languages in Economics" 
http://economics.sas.upenn.edu/~jesusfv/comparison_languages.pdf

The original code in various programming languages is at
https://github.com/jesusfv/Comparison-Programming-Languages-Economics

We use the Haskell library REPA for our arrays
https://hackage.haskell.org/package/repa

## Speed

Right now own my maching the Haskell code compiled with GHC's LLVM 
backend runs about three times slower than Auroba-Villaverde's C++ code
compiled with GCC.

## Building

The easiest way to build is using Cabal

    cabal install --only-dependencies
    cabal build

will build the package. Our default configuration uses GHC's LLVM
backend. If you do not have LLVM installed edit the '.cabal' file 
to remove the option '-fllvm'.

To run the program

    cabal run

## Contact

Author: Jyotirmoy Bhattacharya, jyotirmoy@jyotirmoy.net

