local class = require 'ext.class'
local table = require 'ext.table'

local Roe = require 'roe'
local ADM1D3VarRoe = class(Roe)

local ADM1D3Var = require 'adm1d3var'

function ADM1D3VarRoe:init(args)
	args = table(args)
	args.equation = ADM1D3Var(args)

	ADM1D3VarRoe.super.init(self, args)
end

local State = class(Roe.State)

function State:init(...)
	State.super.init(self, ...)
	for i=1,#self do
		self[i].alpha = 0
		self[i].g_xx = 0
	end
end

function State:clone()
	local dst = State.super.clone(self)
	for i=1,#self do
		dst[i].alpha = self[i].alpha
		dst[i].g_xx = self[i].g_xx
	end
	return dst
end

function State.__add(a,b)
	local c = State(#a, #a[1])
	for i=1,#a do
		for j=1,#a[1] do
			c[i][j] = a[i][j] + b[i][j]
		end
		c[i].alpha = a[i].alpha + b[i].alpha
		c[i].g_xx = a[i].g_xx + b[i].g_xx
	end
	return c
end

function State.__mul(a,b)
	local function is(x) return type(x) == 'table' and x.isa and x:isa(State) end
	local src = is(a) and a or b
	local c = State(#src, #src[1])
	for i=1,#src do
		for j=1,#src[1] do
			local aij = type(a) == 'number' and a or a[i][j]
			local bij = type(b) == 'number' and b or b[i][j]
			c[i][j] = aij * bij
		end
		c[i].alpha = (type(a) == 'number' and a or a[i].alpha) * (type(b) == 'number' and b or b[i].alpha)
		c[i].g_xx = (type(a) == 'number' and a or a[i].g_xx) * (type(b) == 'number' and b or b[i].g_xx)
	end
	return c
end

ADM1D3VarRoe.State = State

function ADM1D3VarRoe:addSourceToDerivCell(dq_dts, i)
	local A_x, D_xxx, KTilde_xx = unpack(self.qs[i])
	local alpha = self.qs[i].alpha
	local g_xx = self.qs[i].g_xx
	local f = self.equation.calc.f(alpha)
	local dalpha_f = self.equation.calc.dalpha_f(alpha)
	
	dq_dts[i].alpha = dq_dts[i].alpha - alpha * alpha * f * KTilde_xx / (g_xx * sqrt(g_xx))
	dq_dts[i].g_xx = dq_dts[i].g_xx - 2 * alpha * KTilde_xx / sqrt(g_xx)
end

return ADM1D3VarRoe
