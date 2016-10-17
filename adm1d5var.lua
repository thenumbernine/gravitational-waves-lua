--[[
based on http://arxiv.org/pdf/gr-qc/9609015v2.pdf
which itself doesn't specify the formalism, just the equations.
If you follow the book, it explains how to re-cast those same equations as the proper formalism (as I use in adm1d3to5var.lua)


hyperbolic formalism:

state vector: [alpha gamma_xx A_x D_xxx K_xx]	<- even though alpha and gamma_xx are already represented by A_x and D_xxx ...
fluxes vector: alpha K_xx * [0, 0, f/gamma_xx, 1, A_x/K_xx]
source vector: alpha / gamma_xx * [-alpha f K_xx, -2 K_xx gamma_xx, 0, 0, A_x D_xxx - K_xx^2]

alpha,t + (0),x = -alpha^2 f K_xx / gamma_xx
gamma_xx,t + (0),x = -2 alpha K_xx
A_x,t + (alpha K_xx f / gamma_xx),x = 0
D_xxx,t + (alpha K_xx),x = 0
K_xx,t + (alpha A_x),x = alpha / gamma_xx (A_x D_xxx - K_xx^2)

--]]

local class = require 'ext.class'
local table = require 'ext.table'
local Equation = require 'equation'

local ADM1D5Var = class(Equation)
ADM1D5Var.name = 'ADM 1D 5-Var'

ADM1D5Var.numStates = 5 

function ADM1D5Var:init(args, ...)

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

	local dalpha_f = f:diff(args.f_param):simplify()
	self.calc.dalpha_f = dalpha_f:compile{args.f_param}
end

do
	local q = function(self,i) return self.qs[i] end
	local alpha = q:_(1)
	local gamma_xx = q:_(2)
	local A_x = q:_(3)
	local D_xxx = q:_(4)
	local K_xx = q:_(5)
	local volume = alpha * math.sqrt:o(gamma_xx)
	ADM1D5Var:buildGraphInfos{
		{alpha = alpha},
		{A_x = A_x},
		{gamma_xx = gamma_xx},
		{D_xxx = D_xxx},
		{K_xx = K_xx},
		{volume = volume},
	}
end

function ADM1D5Var:initCell(sim,i)
	local x = sim.xs[i]
	local alpha = self.calc.alpha(x)
	local gamma_xx = self.calc.gamma_xx(x)
	local A_x = self.calc.dx_alpha(x) / self.calc.alpha(x)
	local D_xxx = 1/2 * self.calc.dx_gamma_xx(x)
	local K_xx = self.calc.K_xx(x)
	return {alpha, gamma_xx, A_x, D_xxx, K_xx}
end

function ADM1D5Var:calcRoeValues(qL, qR)
	local alpha, gamma_xx, A_x, D_xxx, K_xx = ADM1D5Var.super.calcRoeValues(self, qL, qR)
	local f = self.calc.f(alpha)
	return alpha, gamma_xx, A_x, D_xxx, K_xx, f
end

-- [[
function ADM1D5Var:fluxMatrixTransform(solver, m, v)
	local alpha, gamma_xx, A_x, D_xxx, K_xx, f = table.unpack(m)
	local v1, v2, v3, v4, v5 = table.unpack(v)
	return {
		0,
		0,
		v1*f*K_xx/gamma_xx - v2*alpha*f*K_xx/gamma_xx^2 + v5*alpha*f/gamma_xx,
		v1*K_xx + v5*alpha,
		v1*A_x + v3*alpha
	}
end
--]]
--[[fixme
function ADM1D5Var:eigenLeftTransform(solver, m, v)
	local alpha, gamma_xx, A_x, D_xxx, K_xx, f = table.unpack(m)
	local v1, v2, v3, v4, v5 = table.unpack(v)
	return {
		v1 * (gamma_xx * A_x / f - K_xx * math.sqrt(gamma_xx / f)) / (2 * alpha) + v3 * gamma_xx / (2 * f) - v5 * .5 * math.sqrt(gamma_xx / f),
		v1 / alpha,
		-v1 * (gamma_xx * A_x) / (alpha * f) - v3 * gamma_xx / f + v4,
		v2,
		v1 * (gamma_xx * A_x / f + K_xx * math.sqrt(gamma_xx / f)) / (2 * alpha) + v3 * gamma_xx / (2 * f) + v5 * .5 * math.sqrt(gamma_xx / f)
	}
end
--]]

function ADM1D5Var:calcMaxEigenvalue(alpha, gamma_xx)
	local f = self.calc.f(alpha)
	local lambda = alpha * math.sqrt(f / gamma_xx)
	return lambda
end

function ADM1D5Var:calcEigenvaluesFromCons(alpha, gamma_xx, A_x, D_xxx, K_xx, f)
	local lambda = self:calcMaxEigenvalue(alpha, gamma_xx)
	return -lambda, 0, 0, 0, lambda
end

function ADM1D5Var:calcEigenBasis(lambdas, evr, evl, dF_dU, alpha, gamma_xx, A_x, D_xxx, K_xx, f)
	local sqrt_f = math.sqrt(f)
	local sqrt_g = math.sqrt(gamma_xx)
	local lambda = alpha * sqrt_f / sqrt_g 
	fill(lambdas, -lambda, 0, 0, 0, lambda)
	
	-- row-major, math-indexed
	if dF_dU then
		fill(dF_dU, alpha, gamma_xx, A_x, D_xxx, K_xx, f)
		--[[
		fill(dF_dU,
			{0,0,0,0,0},
			{0,0,0,0,0},
			{f*K_xx/gamma_xx, -alpha*f*K_xx/gamma_xx^2, 0,0, alpha*f/gamma_xx},
			{K_xx,0,0,0,alpha},
			{A_x,0,alpha,0,0}
		)
		--]]
	end
	--[[ where did I get this from? probably eigenvector decomposition on the flux
	fill(evr,	
		{0,			alpha,	0,	0,	0			},	-- alpha
		{0,			0,		0,	1,	0			},	-- gamma_xx
		{f/gamma_xx,		-A_x,		0,	0,	f/gamma_xx			},	-- A_x
		{1,			0,		1,	0,	1			},	-- D_xxx
		{-sqrt_f/sqrt_g,-K_xx,		0,	0,	sqrt_f/sqrt_g	}	-- K_xx
	)
	fill(evl,
		{(gamma_xx * A_x / f - K_xx * math.sqrt(gamma_xx / f)) / (2 * alpha), 0, gamma_xx / (2 * f), 0, -.5 * math.sqrt(gamma_xx / f)}, 
		{1 / alpha, 0, 0, 0, 0}, 
		{-(gamma_xx * A_x) / (alpha * f), 0, -gamma_xx / f, 1, 0}, 
		{0, 1, 0, 0, 0}, 
		{(gamma_xx * A_x / f + K_xx * math.sqrt(gamma_xx / f)) / (2 * alpha), 0, gamma_xx / (2 * f), 0, .5 * math.sqrt(gamma_xx / f)} 
	)
	-- note that because we have zero eigenvalues that the eigendecomposition cannot reconstruct the flux matrix
	--]]
	-- [[ here's from the left eigenvectors
	fill(evl,
		{0, 0, -1/sqrt_g, 0, sqrt_f/gamma_xx},	-- math.sqrt(f) K_xx / gamma_xx - A_x / math.sqrt(gamma_xx)
		{1, 0, 0, 0, 0},								-- alpha
		{0, 1, 0, 0, 0},								-- gamma_xx
		{0, 0, 1, -f/gamma_xx, 0},						-- A_x - f D_xxx / gamma_xx
		{0, 0, 1/sqrt_g, 0, sqrt_f/gamma_xx}	-- math.sqrt(f) K_xx / gamma_xx + A_x / math.sqrt(gamma_xx)
	)
	fill(evr,	
		{0, 1, 0, 0, 0},
		{0, 0, 1, 0, 0},
		{-.5 * sqrt_g, 0, 0, 0, .5 * sqrt_g},
		{-.5 * sqrt_g * gamma_xx / f, 0, 0, -gamma_xx / f, .5 * sqrt_g * gamma_xx / f},
		{.5 * gamma_xx / sqrt_f, 0, 0, 0, .5 * gamma_xx / sqrt_f}
	)
	--]]
end	

function ADM1D5Var:sourceTerm(sim, qs)
	local source = sim:newState()
	for i=1,sim.gridsize do
		local alpha, gamma_xx, A_x, D_xxx, K_xx = unpack(qs[i])
		local f = self.calc.f(alpha)
		source[i][1] = -alpha * alpha * f * K_xx / gamma_xx
		source[i][2] = -2 * alpha * K_xx
		source[i][5] = alpha / gamma_xx * (A_x * D_xxx - K_xx * K_xx)
	end
	return source
end

return ADM1D5Var
