local class = require 'ext.class'
local Roe = require 'roe'
local Euler1D = require 'euler1d'
local Euler1DRoe = class(Roe)
Euler1DRoe.equation = Euler1D()
return Euler1DRoe
