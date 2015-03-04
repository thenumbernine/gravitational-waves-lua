--[[
3-var system extrapolated to 5-vars so the integration of alpha and g are not separate
--]]

require 'ext'
local Simulation = require 'simulation'

local ADM1D3to5VarSim = class(Simulation)

ADM1D3to5VarSim.numStates = 5

-- initial conditions
function ADM1D3to5VarSim:init(args, ...)
	ADM1D3to5VarSim.super.init(self, args, ...)

	local symmath = require 'symmath'

	local x = assert(args.x)

	local h = symmath.clone(assert(args.h)):simplify()
	self.calc_h = h:compile{x}
	
	local dx_h = h:diff(x):simplify()
	self.calc_dx_h = dx_h:compile{x}
	
	local d2x_h = dx_h:diff(x):simplify()
	self.calc_d2x_h = d2x_h:compile{x}

	local g = symmath.clone(assert(args.g)):simplify()
	self.calc_g = g:compile{x}

	local dx_g = g:diff(x):simplify()
	self.calc_dx_g = dx_g:compile{x}

	local alpha = symmath.clone(assert(args.alpha)):simplify()
	self.calc_alpha = alpha:compile{x}

	local dx_alpha = alpha:diff(x):simplify()
	self.calc_dx_alpha = dx_alpha:compile{x}

	local f_param = assert(args.f_param)

	local f = symmath.clone(assert(args.f)):simplify()
	self.calc_f = f:compile{f_param}

	local dalpha_f = f:diff(f_param):simplify()
	self.calc_dalpha_f = dalpha_f:compile{f_param}

	local get_state = index:bind(self.qs)
	local get_alpha = get_state:index(1)
	local get_g = get_state:index(2)
	local get_A = get_state:index(3)
	local get_D = get_state:index(4)
	local get_KTilde = get_state:index(5)
	local get_K = get_KTilde / sqrt:compose(get_g)
	self.graphInfos = {
		{viewport={0/3, 0/3, 1/3, 1/3}, getter=get_alpha, name='alpha', color={1,0,1}},
		{viewport={0/3, 1/3, 1/3, 1/3}, getter=get_A, name='A', color={0,1,0}},
		{viewport={1/3, 0/3, 1/3, 1/3}, getter=get_g, name='g', color={.5,.5,1}},
		{viewport={1/3, 1/3, 1/3, 1/3}, getter=get_D, name='D', color={1,1,0}},
		{viewport={2/3, 0/3, 1/3, 1/3}, getter=get_K, name='K', color={0,1,1}},
		{viewport={2/3, 1/3, 1/3, 1/3}, getter=get_alpha * sqrt:compose(get_g), name='volume', color={0,1,1}},
		{viewport={0/3, 2/3, 1/3, 1/3}, getter=log:compose(index:bind(self.eigenbasisErrors)), name='log eigenbasis error', color={1,0,0}, range={-30, 30}},
		{viewport={1/3, 2/3, 1/3, 1/3}, getter=log:compose(index:bind(self.fluxMatrixErrors)), name='log reconstuction error', color={1,0,0}, range={-30, 30}},
	}

	local function buildField(call)
		return function(i, v)
			local avgQ = {}
			for j=1,self.numStates do 
				avgQ[j] = (self.qs[i-1][j] + self.qs[i][j]) / 2
			end
			local alpha, g, A, D, KTilde = unpack(avgQ)		
			local f = self.calc_f(alpha)
			
			return {call(alpha, f, g, A, D, KTilde, unpack(v))}
		end
	end

	self.fluxTransform = buildField(function(alpha, f, g, A, D, KTilde, ...)
		local v1, v2, v3, v4, v5 = ... 
		return
			0,
			0,
			v5 * alpha * f / sqrt(g),
			v5 * 2 * alpha / sqrt(g),
			v3 * alpha / sqrt(g)
	end)
	self.eigenfields = buildField(function(alpha, f, g, A, D, KTilde, ...)
		local v1, v2, v3, v4, v5 = ... 
		return
			v3 / (2 * f) - v5 / (2 * sqrt(f)),	-- first column so it lines up with the min eigenvalue
			v1,
			v2,
			-2*v3/f + v4,
			v3 / (2 * f) + v5 / (2 * sqrt(f))
	end)
	self.eigenfieldsInverse = buildField(function(alpha, f, g, A, D, KTilde, ...)
		local v1, v2, v3, v4, v5 = ...
		return
			v2,
			v3,
			(v1 + v5) * f,
			2 * v1 + v4 + 2 * v5,
			sqrt(f) * (-v1 + v5)
	end)
end

function ADM1D3to5VarSim:initCell(i)
	local x = self.xs[i]
	local alpha = self.calc_alpha(x)
	local g = self.calc_g(x)
	local A = self.calc_dx_alpha(x) / self.calc_alpha(x)
	local D = 1/2 * self.calc_dx_g(x)
	local K = -self.calc_d2x_h(x) / sqrt(self.calc_g(x))
	local KTilde = K / sqrt(g)
	return {alpha, g, A, D, KTilde}
end

function ADM1D3to5VarSim:calcInterfaceEigenBasis(i)
	local avgQ = {}
	for j=1,self.numStates do 
		avgQ[j] = (self.qs[i-1][j] + self.qs[i][j]) / 2
	end
	local alpha, g, A, D, KTilde = unpack(avgQ)
	local f = self.calc_f(alpha)
	local lambda = alpha * sqrt(f / g)		
	self.eigenvalues[i] = {-lambda, 0, 0, 0, lambda}
end

function ADM1D3to5VarSim:addSourceToDerivCell(dq_dts, i)
	local alpha, g, A, D, KTilde = unpack(self.qs[i])
	local f = self.calc_f(alpha)
	local dalpha_f = self.calc_dalpha_f(alpha)
	
	dq_dts[i][1] = dq_dts[i][1] - alpha * alpha * f * KTilde / (g * sqrt(g))
	dq_dts[i][2] = dq_dts[i][2] - 2 * alpha * KTilde / sqrt(g)
end

return ADM1D3to5VarSim

