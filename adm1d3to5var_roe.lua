local class = require 'ext.class'
local table = require 'ext.table'

local ADM1D3to5Var = require 'adm1d3to5var'

local Roe = require 'roe'

local ADM1D3to5VarRoe = class(Roe)

function ADM1D3to5VarRoe:init(args)
	args = table(args)
	args.equation = ADM1D3to5Var(args)
	ADM1D3to5VarRoe.super.init(self, args)
end

-- TODO move all 'addSourceToDerivCell' to the Equation
function ADM1D3to5VarRoe:addSourceToDerivCell(dq_dts,i)
	local alpha, g_xx, A_x, D_xxx, KTilde_xx = unpack(self.qs[i])
	local f = self.equation.calc.f(alpha)
	local dalpha_f = self.equation.calc.dalpha_f(alpha)
	
	dq_dts[i][1] = dq_dts[i][1] - alpha * alpha * f * KTilde_xx / (g_xx * sqrt(g_xx))
	dq_dts[i][2] = dq_dts[i][2] - 2 * alpha * KTilde_xx / sqrt(g_xx)
end

return ADM1D3to5VarRoe
