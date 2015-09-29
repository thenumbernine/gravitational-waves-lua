local class = require 'ext.class'
local table = require 'ext.table'
local Roe = require 'roe'
local BSSNOK1D = require 'bssnok1d'
local BSSNOK1DRoe = class(Roe)
function BSSNOK1DRoe:init(args)
	args = table(args)
	args.equation = BSSNOK1D(args)
	BSSNOK1DRoe.super.init(self, args)
end
return BSSNOK1DRoe
