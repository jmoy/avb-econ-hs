## Introduction
A Haskell implementation of the problem described
in Aruoba and Fernández-Villaverde's paper "A Comparison of Programming
Languages in Economics" 
http://economics.sas.upenn.edu/~jesusfv/comparison_languages.pdf

The original code in various programming languages is at
https://github.com/jesusfv/Comparison-Programming-Languages-Economics

We use the Haskell library REPA and Data.Vector.Unboxed for our arrays
https://hackage.haskell.org/package/repa
The ST monad allows us to use mutable vectors in performance-critical
parts while keeping the rest of the program mutation free.

## Speed

Right now on my machine the Haskell code compiled with GHC's LLVM 
backend takes about 80% more time than Auroba-Villaverde's C++ code
compiled with GCC.

## Building

The easiest way to build is using Cabal

    cabal install --only-dependencies
    cabal build

will build the package. 

Our default configuration uses GHC's LLVM
backend. If you do not have LLVM installed run

    cabal configure -f-llvm

before building.

To run the program

    cabal run

## Contact

Author: Jyotirmoy Bhattacharya, jyotirmoy@jyotirmoy.net

