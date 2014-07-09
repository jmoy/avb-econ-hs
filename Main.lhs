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

\section{Parameters}
\begin{code}

ninfnty::Double
ninfnty=read "-Infinity"

aalpha::Double
aalpha = (1.0/3.0)     --Elasticity of output w.r.t. capital

bbeta::Double
bbeta  = 0.95    -- Discount factor

vProductivity::V.Vector Double
vProductivity = V.fromList l
  where
    l = [0.9792, 0.9896, 1.0000, 1.0106, 1.0212]

mTransition::Array U DIM2 Double
mTransition   = fromListUnboxed (ix2 5 5) (
  [0.9727, 0.0273, 0.0000, 0.0000, 0.0000,
   0.0041, 0.9806, 0.0153, 0.0000, 0.0000,
   0.0000, 0.0082, 0.9837, 0.0082, 0.0000,
   0.0000, 0.0000, 0.0153, 0.9806, 0.0041,
   0.0000, 0.0000, 0.0000, 0.0273, 0.9727])


capitalSteadyState::Double
capitalSteadyState = (aalpha*bbeta)**(1/(1-aalpha))

outputSteadyState::Double
outputSteadyState = capitalSteadyState**aalpha

consumptionSteadyState::Double
consumptionSteadyState = outputSteadyState-capitalSteadyState               

nGridCapital::Int
nGridCapital = 17800

vGridCapital::V.Vector Double
vGridCapital = vec
               where
                 start = 0.5*capitalSteadyState
                 step = capitalSteadyState/(fromIntegral nGridCapital)
                 vec = V.enumFromStepN start step nGridCapital

nGridProductivity::Int
nGridProductivity = V.length vProductivity

\end{code}

We precompute output for all capital and productivity levels.
\begin{code}

mOutput::Array U DIM2 Double
mOutput = computeS $ fromFunction (ix2 nGridCapital nGridProductivity) fn
  where
    fn (Z:.i:.j) = (k**aalpha)*p
      where
        k = vGridCapital `V.unsafeIndex` i  
        p = vProductivity `V.unsafeIndex` j

\end{code}

\section{Optimization}
\subsection{The value function}

\begin{code}

{-# INLINE compute_vf #-}
compute_vf::Array U DIM2 Double -- The expected value function
                                --  indexed by capital stock and
                                --  and previous period's productivity        
            ->Int               -- Index of this period's capital
            ->Int               -- Index of this period's productivity
            ->Int               -- Index of next period's capital
            ->Double            -- Value
compute_vf evf cap prod nxt 
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
    go !v !s  =  
      if s==end then
        (s-1,v)
      else 
        let ns = s+1 
            ky = keyfn ns in 
        if ky<=v then
          (s,v)
        else
          go ky (s+1)

\end{code}

\subsection{Policies}
Find the optimal policy for a given level of capital and productivity.

\begin{code}

policy::Int                     -- Index of capital stock
        ->Int                   -- Index of productivity
        ->Int                   -- Index of next period's capital from
                                --  where to start search
        ->Array U DIM2 Double   -- Expected value function
        ->(Int,Double)          -- Index of optimal value of next
                                --  period's capital and the corresponding
                                --  value of the value function
policy cap prod start evf = 
  findPeak fn start nGridCapital
  where
    fn nxt = compute_vf evf cap prod nxt 

\end{code}

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

writePolicy::forall s. Array U DIM2 Double      -- Expected value function
             -> M.MVector s (Double,Double)     -- Vector to store optimal
                                                --  values and policies, it is
                                                --  actuall a 
                                                --  (nGridCapital*nGridProductivity)
                                                --  matrix in row-major order
             -> Int                             -- Index of productivity.
             -> ST s ()
writePolicy evf mv prod = update 0 0
  where
    ix i = i*nGridProductivity+prod
    update cap start = do
      let (n,v) = policy cap prod start evf
      let k = vGridCapital `V.unsafeIndex` n
      M.unsafeWrite mv (ix cap) (k,v)
      case cap==(nGridCapital-1) of
        True -> return ()
        False ->update (cap+1) n 
\end{code}
 
\section{Value function iteration}
\begin{code}
data DPState = DPState {vf::Array U DIM2 Double,
                        pf::Array U DIM2 Double}
               
iterDP::DPState->DPState
iterDP s = DPState {vf = nvf,pf =npf}
  where
    evf = mmultS (vf s) (transpose2S mTransition)
    bestpv = V.create $ do
      v <- M.new (nGridCapital*nGridProductivity)
      mapM_ (writePolicy evf v) [0..(nGridProductivity-1)]
      return v
    (npf',nvf')= V.unzip bestpv
    npf = fromUnboxed (Z:.nGridCapital:.nGridProductivity) npf'
    nvf = fromUnboxed (Z:.nGridCapital:.nGridProductivity) nvf'
          
supdiff::Array U DIM2 Double->Array U DIM2 Double->Double
supdiff v1 v2 = foldAllS max ninfnty $ R.map abs (v1 -^ v2)
 
\end{code}

\section{Drivers}
\begin{code}
initstate::DPState
initstate = DPState {vf=z,pf=z}
  where
    z = fromUnboxed (Z:.nGridCapital:.nGridProductivity) v 
    v = V.replicate (nGridCapital*nGridProductivity) 0.0
    
printvf::Array U DIM2 Double->IO()
printvf v = mapM_ go [(i,j)|i<-[0,100..(nGridCapital-1)],
                       j<-[0..(nGridProductivity-1)]]
  where
    go (i,j) = (printf "%g\t%g\t%g\n" 
                (vGridCapital V.! i)
                (vProductivity V.! j)
                (v ! (Z:.i:.j)))
               

tolerance::Double
tolerance = 1e-7

maxIter::Int
maxIter=1000


main::IO()
main = do
  go 1 initstate
  where
    go::Int->DPState->IO()
    go !count !s = 
      let ns = iterDP s
          d = supdiff (vf s) (vf ns) 
          putLog::IO()
          putLog = printf "Iteration = %d, Sup Diff = %.6g\n" count d in
      if (d <tolerance) || (count>maxIter) then do
        putLog
        printf "My check = %.6g\n" (pf ns ! ix2 999 2)
        --printvf (vf s)
      else do
        when (count `mod` 10==0) putLog
        go (count+1) ns
\end{code}
\end{document}
