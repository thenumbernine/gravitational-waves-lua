local class = require 'ext.class'
local HLL = require 'hll'
local Euler1DHLL = class(HLL)
Euler1DHLL.name = 'Euler 1D HLL'
Euler1DHLL.equation = require 'euler1d'()
return Euler1DHLL

