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

local function buildField(call)
	return function(self, i, v)
		local v1, v2, v3, v4, v5 = unpack(v)
		
		local avgQ = {}
		for j=1,self.numStates do 
			avgQ[j] = (self.qs[i-1][j] + self.qs[i][j]) / 2
		end
		local alpha, g_xx, A_x, D_xxx, K_xx = unpack(avgQ)
		local f = self.equation.calc.f(alpha)
		return {call(alpha, f, g_xx, A_x, D_xxx, K_xx, v1, v2, v3, v4, v5)}
	end
end

ADM1D5VarRoe.fluxTransform = buildField(function(alpha, f, g_xx, A_x, D_xxx, K_xx, v1, v2, v3, v4, v5)
	return 
		0,
		0,
		v1*f*K_xx/g_xx - v2*alpha*f*K_xx/g_xx^2 + v5*alpha*f/g_xx,
		v1*K_xx + v5*alpha,
		v1*A_x + v3*alpha
end)
--[[fixme
ADM1D5VarRoe.eigenfields = buildField(function(alpha, f, g_xx, A_x, D_xxx, K_xx, v1, v2, v3, v4, v5)
	return 
		v1 * (g_xx * A_x / f - K_xx * sqrt(g_xx / f)) / (2 * alpha) + v3 * g_xx / (2 * f) - v5 * .5 * sqrt(g_xx / f),
		v1 / alpha,
		-v1 * (g_xx * A_x) / (alpha * f) - v3 * g_xx / f + v4,
		v2,
		v1 * (g_xx * A_x / f + K_xx * sqrt(g_xx / f)) / (2 * alpha) + v3 * g_xx / (2 * f) + v5 * .5 * sqrt(g_xx / f) 
end),
--]]

function ADM1D5VarRoe:addSourceToDerivCell(dq_dts, i)
	local alpha, g_xx, A_x, D_xxx, K_xx = unpack(self.qs[i])
	local f = self.equation.calc.f(alpha)
	dq_dts[i][1] = dq_dts[i][1] - alpha * alpha * f * K_xx / g_xx
	dq_dts[i][2] = dq_dts[i][2] - 2 * alpha * K_xx
	dq_dts[i][5] = dq_dts[i][5] + alpha * (A_x * D_xxx - K_xx * K_xx) / g_xx
end

return ADM1D5VarRoe

