local class = require 'ext.class'
local table = require 'ext.table'

local Roe = require 'roe'

local ADM1D5VarRoe = class(Roe)

ADM1D5VarRoe.name = 'ADM 1D 5-Var Roe'

function ADM1D5VarRoe:init(args)
	args = table(args)
	args.equation = require 'adm1d5var'(args)
	ADM1D5VarRoe.super.init(self, args)
end

function ADM1D5VarRoe:addSourceToDerivCell(dq_dts, i)
	local alpha, g_xx, A_x, D_xxx, K_xx = unpack(self.qs[i])
	local f = self.equation.calc.f(alpha)
	dq_dts[i][1] = dq_dts[i][1] - alpha * alpha * f * K_xx / g_xx
	dq_dts[i][2] = dq_dts[i][2] - 2 * alpha * K_xx
	dq_dts[i][5] = dq_dts[i][5] + alpha * (A_x * D_xxx - K_xx * K_xx) / g_xx
end

return ADM1D5VarRoe
