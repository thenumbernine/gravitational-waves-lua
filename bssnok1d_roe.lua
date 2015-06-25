local class = require 'ext.class'
local table = require 'ext.table'

local Roe = require 'roe'

local BSSNOK1DRoe = class(Roe)

function BSSNOK1DRoe:init(args)
	args = table(args)
	args.equation = require 'bssnok1d'(args)
	BSSNOK1DRoe.super.init(self, args)
end

function BSSNOK1DRoe:addSourceToDerivCell(dq_dts, i)
	local alpha, phi, A_x, Phi_x, K, ATilde_xx = unpack(self.qs[i])
	local f = self.equation.calc.f(alpha)
	dq_dts[i][1] = dq_dts[i][1] - alpha * alpha * f * K
	dq_dts[i][2] = dq_dts[i][2] - alpha * K / 6
end

return BSSNOK1DRoe

