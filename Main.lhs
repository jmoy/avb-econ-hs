\documentclass{article}
%include polycode.fmt
\begin{document}
\section{Imports}

\begin{code}
{-# LANGUAGE BangPatterns, ExplicitForAll #-}
{- (c) Jyotirmoy Bhattacharya, 2014, jyotirmoy@@jyotirmoy.net -}
{- Licensed under GPL v3 -}

module Main (main) where

import Control.Monad
import Control.Monad.ST
import Data.Array.Repa as R
import Data.Array.Repa.Algorithms.Matrix
import Prelude as P
import Text.Printf
import qualified Data.Vector.Unboxed as V
import qualified Data.Vector.Unboxed.Mutable as M

\end{code}

\section{Type Definitions}

\begin{code}
type Vector = V.Vector Double
type Matrix = Array U DIM2 Double
\end{code}

\section{Parameters}

\begin{code}
ninfnty::Double
ninfnty=read "-Infinity"

aalpha::Double
aalpha = (1.0/3.0)     -- Elasticity of output w.r.t. capital

bbeta::Double
bbeta  = 0.95          -- Discount factor

vProductivity::Vector
vProductivity = V.fromList l
  where
    l = [0.9792, 0.9896, 1.0000, 1.0106, 1.0212]

mTransition::Matrix
mTransition   = fromListUnboxed (ix2 5 5) 
  [0.9727, 0.0273, 0.0000, 0.0000, 0.0000,
   0.0041, 0.9806, 0.0153, 0.0000, 0.0000,
   0.0000, 0.0082, 0.9837, 0.0082, 0.0000,
   0.0000, 0.0000, 0.0153, 0.9806, 0.0041,
   0.0000, 0.0000, 0.0000, 0.0273, 0.9727]


capitalSteadyState::Double
capitalSteadyState = (aalpha*bbeta)**(1/(1-aalpha))

outputSteadyState::Double
outputSteadyState = capitalSteadyState**aalpha

consumptionSteadyState::Double
consumptionSteadyState = outputSteadyState-capitalSteadyState               

nGridCapital::Int
nGridCapital = 17800

vGridCapital::Vector
vGridCapital = vec
  where
    start = 0.5*capitalSteadyState
    step = capitalSteadyState/(fromIntegral nGridCapital)
    vec = V.enumFromStepN start step nGridCapital

nGridProductivity::Int
nGridProductivity = V.length vProductivity
\end{code}

All the two-dimensional matrices we build have \texttt{nGridCapital}
rows and \texttt{nGridProductivity} columns.

\begin{code}
matShape::DIM2
matShape = ix2 nGridCapital nGridProductivity
\end{code}

We precompute output for all capital and productivity levels.

\begin{code}
mOutput::Matrix
mOutput = computeS $ fromFunction matShape fn
  where
    fn (Z:.i:.j) = (k**aalpha)*p
      where
        k = vGridCapital `V.unsafeIndex` i  
        p = vProductivity `V.unsafeIndex` j
\end{code}

\section{Optimization}
\subsection{The value function}

\begin{code}
{-# INLINE computeVf #-}
computeVf::Matrix       -- The expected value function
                        --  indexed by capital stock and
                        --  and previous period's productivity        
           ->Int        -- Index of this period's capital
           ->Int        -- Index of this period's productivity
           ->Int        -- Index of next period's capital
           ->Double     -- Value
computeVf evf cap prod nxt 
  = (1-bbeta)*(log c)+bbeta*ev
  where
    y = mOutput `unsafeIndex` (ix2 cap prod)
    k' = vGridCapital `V.unsafeIndex` nxt
    c = y - k'
    ev = evf `unsafeIndex` (ix2 nxt prod)
\end{code}

\subsection{Maximization}

Given a function over integers and an integer range,
find the maximum of the function assuming that it is
single-peaked.

\begin{code}
{-# INLINE findPeak #-}
findPeak::(Int->Double)         -- function by which indices are ranked
          ->Int                 -- starting index for search
          ->Int                 -- 1+the last index to be searched
          ->(Int,Double)        -- the index at which the function peaks
                                --   and the value fo the function at the
                                --   peak
findPeak keyfn start end = go (keyfn start) start
  where
    go !oldv !olds  =  
      if olds==end-1 then
        (olds,oldv)
      else 
        let 
          news = olds+1
          newv = keyfn news in 
        if newv<=oldv then
          (olds,oldv)
        else
          go newv news
\end{code}

\subsection{Policies}

Find the best policies for each level of capital for a 
given level of productivity. 

We make use of the monotonicity of the policy function to
begin our search for each level of capital at the point
where the search for the previous level of capital
succeeded.

The results are stored in a mutable vector. The use of a
mutable vector here is essential for good performance since
the answers are generated for one value of productivity at a
time whereas the final policy and value function are indexed
by the level of capital first. So using immutable vectors
not only incurs the cost of memory allocation and copying
but also the cost of doing something equivalent to a
transpose.

\begin{code}
writePolicy::forall s. Matrix                   -- Expected value function
             -> M.MVector s (Double,Double)     -- Vector to store optimal
                                                --  values and policies, it is
                                                --  actuall a 
                                                --  (nGridCapital*nGridProductivity)
                                                --  matrix in row-major order
             -> Int                             -- Index of productivity.
             -> ST s ()
writePolicy evf mv prod = loop 0 0
  where
    ix i = i*nGridProductivity+prod

    policy cap start
      = findPeak (computeVf evf cap prod) start nGridCapital

    loop cap _ | cap==nGridCapital = return()
    loop cap start = do
      let (n,v) = policy cap start
      let k = vGridCapital `V.unsafeIndex` n
      M.unsafeWrite mv (ix cap) (k,v)
      loop (cap+1) n 
\end{code}
 
\section{Value function iteration}

\begin{code}
iterDP::Matrix          --the old value function
        ->(Matrix,      --the new value function
           Matrix)      --the new policy function
iterDP !vf = (nvf,npf)
  where
    evf = mmultS vf (transpose2S mTransition)
    bestpv = V.create $ do
      v <- M.new (nGridCapital*nGridProductivity)
      mapM_ (writePolicy evf v) [0..(nGridProductivity-1)]
      return v
    (npf',nvf')= V.unzip bestpv
    npf = fromUnboxed matShape npf'
    nvf = fromUnboxed matShape nvf'
          
supdiff::Matrix->Matrix->Double
supdiff v1 v2 = foldAllS max ninfnty $ R.map abs (v1 -^ v2)
\end{code}

\section{Drivers}

\begin{code}
initstate::Matrix
initstate
  = fromUnboxed matShape zeros
  where
    zeros = V.replicate (nGridCapital*nGridProductivity) 0.0
    
tolerance::Double
tolerance = 1e-7

maxIter::Int
maxIter=1000


main::IO()
main = do
  _ <- printf "Output = %.6g, Capital = %.6g, Consumption = %.6g\n" 
       outputSteadyState capitalSteadyState consumptionSteadyState
  go 1 initstate
  where
    go::Int->Matrix->IO()
    go !count !vf = 
      let (nvf,npf) = iterDP vf
          d = supdiff vf nvf 
          putLog::IO()
          putLog = printf "Iteration = %d, Sup Diff = %.6g\n" 
                   count d 
      in
      if (d <tolerance) || (count>maxIter) then do
        putLog
        printf "My check = %.6g\n" (npf ! ix2 999 2)
      else do
        when (count `mod` 10==0) putLog
        go (count+1) nvf
\end{code}
\end{document}
