--[[
3-var system extrapolated to 5-vars so the integration of alpha and gamma_xx are not separate
I think this is no longer any different from ADM1D3Var 
--]]

local class = require 'ext.class'
local table = require 'ext.table'
local Equation = require 'equation'

local ADM1D3to5Var = class(Equation)
ADM1D3to5Var.name = 'ADM 1D 3-to-5-Var'

ADM1D3to5Var.numStates = 5

-- initial conditions
function ADM1D3to5Var:init(args, ...)
	
	local symmath = require 'symmath'
	local function makesym(field)
		return symmath.clone(assert(args[field], "expected to find field "..field)):simplify() 
	end
	
	-- parameters that are variables of symbolic functions
	local x = assert(args.x)
	
	-- parameters that are symbolic functions -- based on coordinates 
	local exprs = table{'alpha', 'gamma_xx', 'K_xx'}:map(function(name)
		return makesym(name), name
	end)
	
	-- derived functions
	exprs.dx_alpha = exprs.alpha:diff(x):simplify()
	exprs.dx_gamma_xx = exprs.gamma_xx:diff(x):simplify()
	
	-- convert from symbolic functions to Lua functions
	self.calc = exprs:map(function(expr, name)
		return expr:compile{x}, name
	end)
	
	-- parameters that are symbolic functions -- based on f's alpha
	-- NOTICE: assign these *after* mapping exprs to calc
	local f = makesym'f'
	local f_param = assert(args.f_param)
	self.calc.f = f:compile{f_param}
	
	local dalpha_f = f:diff(f_param):simplify()
	self.calc.dalpha_f = dalpha_f:compile{f_param}
end

do
	local q = function(self,i) return self.qs[i] end
	local alpha = q:_(1)
	local gamma_xx = q:_(2)
	local a_x = q:_(3)
	local D_g = q:_(4)
	local d_xxx = D_g * gamma_xx / 2
	local KTilde = q:_(5)
	local K_xx = KTilde * math.sqrt:o(gamma_xx)
	local K = K_xx / gamma_xx
	local volume = alpha * math.sqrt:o(gamma_xx)
	ADM1D3to5Var:buildGraphInfos{
		{alpha = alpha},
		{a_x = a_x},
		{gamma_xx = gamma_xx},
		{d_xxx = d_xxx},
		{D_g = D_g},
		{K_xx = K_xx},
		{KTilde = KTilde},
		{K = K},
		{volume = volume},
	}
end

function ADM1D3to5Var:initCell(sim,i)
	local x = sim.xs[i]
	local alpha = self.calc.alpha(x)
	local gamma_xx = self.calc.gamma_xx(x)
	local a_x = self.calc.dx_alpha(x) / self.calc.alpha(x)
	local D_g = self.calc.dx_gamma_xx(x) / self.calc.gamma_xx(x)
	local K_xx = self.calc.K_xx(x)
	local KTilde = K_xx / math.sqrt(gamma_xx)
	return {alpha, gamma_xx, a_x, D_g, KTilde}
end

function ADM1D3to5Var:calcRoeValues(qL, qR)
	local alpha = .5 * (qL[1] + qR[1])
	local gamma_xx = .5 * (qL[2] + qR[2])
	local f = self.calc.f(alpha)
	return alpha, gamma_xx, f
end

function ADM1D3to5Var:calcEigenBasis(lambda, evr, evl, dF_dU, alpha, gamma_xx, f)
	fill(lambda, self:calcEigenvalues(alpha, gamma_xx))
	fill(evl, f)
	fill(evr, f)
	if dF_dU then fill(dF_dU, alpha, gamma_xx, f) end
end

function ADM1D3to5Var:fluxMatrixTransform(solver, m, v)
	local alpha, gamma_xx, f = table.unpack(m)
	local v1, v2, v3, v4, v5 = table.unpack(v)
	local alpha_over_sqrt_gamma_xx = alpha / math.sqrt(gamma_xx)
	return {
		0,
		0,
		v5 * f * alpha_over_sqrt_gamma_xx,
		v5 * 2 * alpha_over_sqrt_gamma_xx,
		v3 * alpha_over_sqrt_gamma_xx
	}
end

function ADM1D3to5Var:eigenLeftTransform(solver, m, v)
	local v1, v2, v3, v4, v5 = table.unpack(v)
	local f = table.unpack(m)
	local sqrt_f = math.sqrt(f)
	return {
		v3 / (2 * f) - v5 / (2 * sqrt_f),	-- first column so it lines up with the min eigenvalue
		v1,
		v2,
		-2*v3/f + v4,
		v3 / (2 * f) + v5 / (2 * sqrt_f)
	}
end

function ADM1D3to5Var:eigenRightTransform(solver, m, v)
	local v1, v2, v3, v4, v5 = table.unpack(v)
	local f = table.unpack(m)
	return {
		v2,
		v3,
		(v1 + v5) * f,
		2 * v1 + v4 + 2 * v5,
		math.sqrt(f) * (-v1 + v5)
	}
end

function ADM1D3to5Var:calcMaxEigenvalue(alpha, gamma_xx)
	local f = self.calc.f(alpha)
	local lambda = alpha * math.sqrt(f / gamma_xx)		
	return lambda
end

function ADM1D3to5Var:calcEigenvalues(alpha, gamma_xx)
	local lambda = self:calcMaxEigenvalue(alpha, gamma_xx)
	return -lambda, 0, 0, 0, lambda
end

function ADM1D3to5Var:calcEigenvaluesFromCons(alpha, gamma_xx, a_x, D_g, KTilde)
	return self:calcEigenvalues(alpha, gamma_xx)
end

function ADM1D3to5Var:sourceTerm(sim, qs)
	local source = sim:newState()
	for i=1,sim.gridsize do
		local alpha, gamma_xx, a_x, D_g, KTilde = unpack(qs[i])
		local f = self.calc.f(alpha)
		local dalpha_f = self.calc.dalpha_f(alpha)
		local K = KTilde / math.sqrt(gamma_xx)
		
		source[i][1] = -alpha * alpha * f * K
		source[i][2] = -2 * alpha * gamma_xx * K
		source[i][3] = ((1/2 * D_g - a_x) * f - alpha * dalpha_f * a_x) * alpha * K
		source[i][4] = (1/2 * D_g - a_x) * 2 * alpha * K
		source[i][5] = (1/2 * D_g - a_x) * a_x * alpha / math.sqrt(gamma_xx)
	end
	return source
end

return ADM1D3to5Var
