\documentclass{article}
\begin{document}
\section{Imports}
\begin{code}
{-# LANGUAGE BangPatterns, ExplicitForAll #-}
{-
(c) Jyotirmoy Bhattacharya, 2014
jyotirmoy@jyotirmoy.net

Licensed under GPL v3
-}

module Main where

import Control.Monad
import Control.Monad.ST
import Data.Array.Repa as R
import Data.Array.Repa.Algorithms.Matrix
import Prelude as P
import Text.Printf
import qualified Data.Vector.Unboxed as V
import qualified Data.Vector.Unboxed.Mutable as M

ninfnty::Double
ninfnty=read "-Infinity"
\end{code}
\section{Parameters}

\begin{code}

aalpha::Double
aalpha = (1.0/3.0)     --Elasticity of output w.r.t. capital

bbeta::Double
bbeta  = 0.95    -- Discount factor

-- Productivity values
vProductivity::Array U DIM1 Double
vProductivity = fromListUnboxed (ix1 $ P.length l) l
  where
    l = [0.9792, 0.9896, 1.0000, 1.0106, 1.0212]

-- Transition matrix
mTransition::Array U DIM2 Double
mTransition   = fromListUnboxed (ix2 5 5) (
  [0.9727, 0.0273, 0.0000, 0.0000, 0.0000,
   0.0041, 0.9806, 0.0153, 0.0000, 0.0000,
   0.0000, 0.0082, 0.9837, 0.0082, 0.0000,
   0.0000, 0.0000, 0.0153, 0.9806, 0.0041,
   0.0000, 0.0000, 0.0000, 0.0273, 0.9727])

\end{code}

\section{Steady State}

\begin{code}

capitalSteadyState::Double
capitalSteadyState = (aalpha*bbeta)**(1/(1-aalpha))

outputSteadyState::Double
outputSteadyState = capitalSteadyState**aalpha

consumptionSteadyState::Double
consumptionSteadyState = outputSteadyState-capitalSteadyState               
nGridCapital::Int
nGridCapital = 17800
\end{code}

\section{Working variables}

We generate the grid of capital.

\begin{code}
vGridCapital::Array U DIM1 Double
vGridCapital = fromUnboxed (ix1 nGridCapital) vec
               where
                 start = 0.5*capitalSteadyState
                 step = capitalSteadyState/(fromIntegral nGridCapital)
                 vec = V.enumFromStepN start step nGridCapital

nGridProductivity::Int
(Z:.nGridProductivity) = extent vProductivity
\end{code}
 
We pre-build output for each point in the grid.

\begin{code}
mOutput::Array U DIM2 Double
mOutput = computeS $ (R.zipWith f xK xP)
  where
    f !k !p = (k**aalpha) * p
    xK = extend (Z:.All:.nGridProductivity) vGridCapital
    xP = extend (Z:.nGridCapital:.All) vProductivity
\end{code}

\section{Maximization}

Compute the value function given the level of this period's 
capital, productivity and next period's capital. We make use
of the expected value function \texttt{evf} that is passed to us.

\begin{code}
{-# INLINE compute_vf #-}
compute_vf::Array U DIM2 Double->Int->Int->Int->Double
compute_vf evf cap prod nxt = v
  where
    y = mOutput `unsafeIndex` (ix2 cap prod)
    k' = vGridCapital `unsafeIndex` (ix1 nxt)
    c = y - k'
    ev = evf `unsafeIndex` (ix2 nxt prod)
    v = (1-bbeta)*(log c)+bbeta*ev 
\end{code}
 
A helper function to compute the peak of a single-peaked function.
It is passed the function itself, and the points in its domain to
consider.

\begin{code}
{-# INLINE findPeak #-}
findPeak::(Int->Double)->[Int]->(Int,Double)
findPeak _ [] = error "Empty argument to findPeak"
findPeak keyfn (x:xs) = go (keyfn x) x xs
  where
    go !v !yp l =  
      case l of
        [] -> (yp,v)
        (y:ys) -> 
          let ky = keyfn y in 
          case ky<=v of
            True -> (yp,v) 
            False ->  go ky y ys
\end{code}


Find the best policy for for a cerain stock of capital
and productivity. We are passed a lower bound for
the domain to search in the parameter \texttt{start}.

\begin{code}
policy::Int->Int->Int->Array U DIM2 Double->(Int,Double)
policy cap prod start evf = 
  findPeak fn [start..(nGridCapital-1)]
  where
    fn nxt = compute_vf evf cap prod nxt 
\end{code}

Find the best policies for each level of capital for a 
given level of productivity. We make use of the 
monotonicity of the policy function to begin our search
for each level of capital at the point where the search
for the previous level of capital succeeded.

\begin{code}
writePolicy::forall s. Array U DIM2 Double
             -> M.MVector s (Double,Double)
             -> Int
             -> ST s ()
writePolicy evf mv prod = update 0 0
  where
    ix i = i*nGridProductivity+prod
    update cap start = do
      let (n,v) = policy cap prod start evf
      let k = vGridCapital `unsafeIndex` (ix1 n)
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
                (vGridCapital ! (ix1 i))
                (vProductivity ! (ix1 j))
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
          d = supdiff (vf s) (vf ns) in
      if (d <tolerance) || (count>maxIter) then do
        printf "My check = %.6g\n" (pf ns ! ix2 999 2)
        --printvf (vf s)
      else do
        when (count `mod` 10==0) $ do
          printf "Iteration = %d, Sup Diff = %.6g\n" count d
        go (count+1) ns
\end{code}
\end{document}