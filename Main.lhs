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
import Kernel

import Control.Monad
import Control.Monad.ST
import Data.Array.Repa as R
import Data.Array.Repa.Algorithms.Matrix
import Prelude as P
import Text.Printf
import qualified Data.Vector.Unboxed as V
import qualified Data.Vector.Unboxed.Mutable as M

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
