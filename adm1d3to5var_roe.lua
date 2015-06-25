local class = require 'ext.class'
local table = require 'ext.table'

local Roe = require 'roe'

local ADM1D3to5VarRoe = class(Roe)

function ADM1D3to5VarRoe:init(args)
	args = table(args)
	args.equation = require 'adm1d3to5var'(args)
	ADM1D3to5VarRoe.super.init(self, args)
end

local function buildField(call)
	return function(self, i, v)
		local avgQ = {}
		for j=1,self.numStates do 
			avgQ[j] = (self.qs[i-1][j] + self.qs[i][j]) / 2
		end
		local alpha, g_xx, A_x, D_xxx, KTilde_xx = unpack(avgQ)		
		local f = self.equation.calc.f(alpha)
		
		return {call(alpha, f, g_xx, A_x, D_xxx, KTilde_xx, unpack(v))}
	end
end

ADM1D3to5VarRoe.fluxTransform = buildField(function(alpha, f, g_xx, A_x, D_xxx, KTilde_xx, ...)
	local v1, v2, v3, v4, v5 = ... 
	return
		0,
		0,
		v5 * alpha * f / sqrt(g_xx),
		v5 * 2 * alpha / sqrt(g_xx),
		v3 * alpha / sqrt(g_xx)
end)
ADM1D3to5VarRoe.eigenfields = buildField(function(alpha, f, g_xx, A_x, D_xxx, KTilde_xx, ...)
	local v1, v2, v3, v4, v5 = ... 
	return
		v3 / (2 * f) - v5 / (2 * sqrt(f)),	-- first column so it lines up with the min eigenvalue
		v1,
		v2,
		-2*v3/f + v4,
		v3 / (2 * f) + v5 / (2 * sqrt(f))
end)
ADM1D3to5VarRoe.eigenfieldsInverse = buildField(function(alpha, f, g_xx, A_x, D_xxx, KTilde_xx, ...)
	local v1, v2, v3, v4, v5 = ...
	return
		v2,
		v3,
		(v1 + v5) * f,
		2 * v1 + v4 + 2 * v5,
		sqrt(f) * (-v1 + v5)
end)

function ADM1D3to5VarRoe:addSourceToDerivCell(dq_dts,i)
	local alpha, g_xx, A_x, D_xxx, KTilde_xx = unpack(self.qs[i])
	local f = self.equation.calc.f(alpha)
	local dalpha_f = self.equation.calc.dalpha_f(alpha)
	
	dq_dts[i][1] = dq_dts[i][1] - alpha * alpha * f * KTilde_xx / (g_xx * sqrt(g_xx))
	dq_dts[i][2] = dq_dts[i][2] - 2 * alpha * KTilde_xx / sqrt(g_xx)
end

return ADM1D3to5VarRoe

