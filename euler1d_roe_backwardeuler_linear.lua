local Euler1D = require 'euler1d'
local class = require 'ext.class'
local RoeBackwardEulerLinear = require 'roe_backwardeuler_linear'
local Euler1DRoeBackwardEulerLinear = class(RoeBackwardEulerLinear)
Euler1DRoeBackwardEulerLinear.equation = Euler1D()
return Euler1DRoeBackwardEulerLinear
