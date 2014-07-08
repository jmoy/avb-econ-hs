{-# LANGUAGE BangPatterns #-}

module Kernel where


import Data.Array.Repa as R
import qualified Data.Vector.Unboxed as V
import Prelude as P

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

mOutput::Array U DIM2 Double
mOutput = computeS $ fromFunction (ix2 nGridCapital nGridProductivity) fn
  where
    fn (Z:.i:.j) = ((vGridCapital `V.unsafeIndex` i)**aalpha) * (vProductivity `V.unsafeIndex` j)

{-# INLINE compute_vf #-}
compute_vf::Array U DIM2 Double->Int->Int->Int->Double
compute_vf evf cap prod nxt = v
  where
    y = mOutput `unsafeIndex` (ix2 cap prod)
    k' = vGridCapital `V.unsafeIndex` nxt
    c = y - k'
    ev = evf `unsafeIndex` (ix2 nxt prod)
    v = (1-bbeta)*(log c)+bbeta*ev 

{-# INLINE findPeak #-}
findPeak::(Int->Double)         -- function by which indices are ranked
          ->Int                 -- starting index for search
          ->Int                 -- 1+the last index to be searched
          ->(Int,Double)        -- the index at which the function peaks
findPeak keyfn start end = go (keyfn start) start
  where
    go !v !s  =  
      if s==end then
        (s,v)
      else 
        let ns = (s+1)
            ky = keyfn ns in 
        if ky<=v then
          (s,v)
        else
          go ky ns


policy::Int->Int->Int->Array U DIM2 Double->(Int,Double)
policy cap prod start evf = 
  findPeak fn start nGridCapital
  where
    fn nxt = compute_vf evf cap prod nxt 

