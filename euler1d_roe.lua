local class = require 'ext.class'
local Roe = require 'roe'
local Euler1DRoe = class(Roe)
Euler1DRoe.equation = require 'euler1d'()
return Euler1DRoe
