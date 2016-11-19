--[[

Based on the book "Introduction to 3+1 Numerical Relativity" and on the paper "Introduction to Numerical Relativity", both by Alcubierre

a_x = (ln alpha),x = alpha,x / alpha
d_xxx = 1/2 gamma_xx,x
D_g = (ln gamma_xx),x = gamma_xx,x / gamma_xx = 2 d_xxx / gamma_xx
KTilde_xx = sqrt(gamma_xx) K_xx

alpha,x = a_x alpha
gamma_xx,x = gamma_xx D_g
K_xx = KTilde_xx / sqrt(gamma_xx)

a_x,t + (alpha f K_xx),x = 0
D_g,t + (2 alpha K_xx),x = 0
K_xx,t + (alpha a_x / gamma_xx),x = alpha (K_xx^2 - a_x D_g / (2 gamma_xx))

...rewritten for KTilde_xx...
a_x,t + (alpha f KTilde_xx / sqrt(gamma_xx)),x = 0
D_g,t + (2 alpha KTilde_xx / sqrt(gamma_xx)),x = 0
KTilde_xx,t + (alpha a_x / sqrt(gamma_xx)),x = 0

now we look at eigenvectors ...

/* [wxMaxima: input   start ] */
loag("eigen");
/* [wxMaxima: input   end   ] */
/* [wxMaxima: input   start ] */
assume(alpha>0,f>0);
/* [wxMaxima: input   end   ] */
/* [wxMaxima: input   start ] */
F : matrix(
[0, 0, alpha *f/sqrt(gamma_xx)],
[0, 0, 2*alpha/sqrt(gamma_xx)],
[alpha/sqrt(gamma_xx), 0, 0]
);
/* [wxMaxima: input   end   ] */
/* [wxMaxima: input   start ] */
results:eigenvectors(F);
/* [wxMaxima: input   end   ] */
/* [wxMaxima: input   start ] */
eigenvalues : diag_matrix(results[1][1][1], results[1][1][3], results[1][1][2]);
/* [wxMaxima: input   end   ] */
/* [wxMaxima: input   start ] */
eigenvectors : transpose(matrix(
results[2][1][1] * f,
results[2][3][1],
results[2][2][1] * f
));
/* [wxMaxima: input   end   ] */
/* [wxMaxima: input   start ] */
determinant(eigenvectors);
/* [wxMaxima: input   end   ] */
/* [wxMaxima: input   start ] */
invert(eigenvectors);
/* [wxMaxima: input   end   ] */

--]]

local class = require 'ext.class'
local Equation = require 'equation'

local ADM1D3Var = class(Equation)
ADM1D3Var.name = 'ADM 1D 3-Var'

ADM1D3Var.numStates = 5
ADM1D3Var.numWaves = 3

-- initial conditions
function ADM1D3Var:init(args, ...)

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
	exprs.dx_gamma_xx = exprs.gamma_xx:diff(x):simplify()
	exprs.dx_alpha = exprs.alpha:diff(x):simplify()

	-- convert from symbolic functions to Lua functions
	self.calc = exprs:map(function(expr, name)
		return expr:compile{x}, name
	end)

	-- parameters that are symbolic functions -- based on f's alpha
	-- NOTICE: assign these *after* mapping exprs to calc
	exprs.f = makesym'f'
	local f_param = assert(args.f_param)
	self.calc.f = exprs.f:compile{f_param}
	
	local dalpha_f = exprs.f:diff(f_param):simplify()
	self.calc.dalpha_f = dalpha_f:compile{f_param}
end

do
	local q = function(self,i) return self.qs[i] end
	local alpha = q:_(1)
	local gamma_xx = q:_(2)
	local a_x = q:_(3)
	local D_g = q:_(4)
	local d_xxx = D_g * gamma_xx / 2
	local KTilde_xx = q:_(5)
	local K_xx = KTilde_xx / math.sqrt:o(gamma_xx)
	local K = K_xx / gamma_xx
	local volume = alpha * math.sqrt:o(gamma_xx)
	ADM1D3Var:buildGraphInfos{
		{alpha = alpha},
		{a_x = a_x},
		{gamma_xx = gamma_xx},
		{d_xxx = d_xxx},
		{D_g = D_g},
		{K_xx = K_xx},
		{KTilde_xx = KTilde_xx},
		{K = K},
		{volume = volume},
	}
end

function ADM1D3Var:initCell(sim,i)
	local x = sim.xs[i]
	local alpha = self.calc.alpha(x)
	local gamma_xx = self.calc.gamma_xx(x)
	local a_x = self.calc.dx_alpha(x) / self.calc.alpha(x)
	local D_g = 1/2 * self.calc.dx_gamma_xx(x)
	local K_xx = self.calc.K_xx(x) 
	local KTilde_xx = K_xx / math.sqrt(gamma_xx)
	return {alpha, gamma_xx, a_x, D_g, KTilde_xx}
end

function ADM1D3Var:calcEigenvalues(alpha, gamma_xx, f)
	local lambda = alpha * math.sqrt(f / gamma_xx)
	return -lambda, 0, lambda
end

-- arithmetic
function ADM1D3Var:calcRoeValues(qL, qR)
	local alpha = (qL[1] + qR[1]) / 2
	local gamma_xx = (qL[2] + qR[2]) / 2
	local f = self.calc.f(alpha)
	return alpha, gamma_xx, f
end

function ADM1D3Var:calcEigenBasis(lambda, evr, evl, dF_dU, alpha, gamma_xx, f)
	-- store eigenvalues
	fill(lambda, self:calcEigenvalues(alpha, gamma_xx, f))
	-- store information needed to build left and right eigenvectors
	-- this is why I need an 'eigenbasis' object - to hold this information once
	fill(evl, f)
	fill(evr, f)
	if dF_dU then fill(dF_dU, alpha, gamma_xx, f) end
end

-- how can the flux depend on gamma_xx and alpha but not the eigenvectors?
function ADM1D3Var:fluxMatrixTransform(solver, m, v)
	local alpha, gamma_xx, f = table.unpack(m)
	local _, _, v1, v2, v3 = table.unpack(v)
	return {
		0,
		0,
		v3 * alpha * f / math.sqrt(gamma_xx),
		v3 * 2 * alpha / math.sqrt(gamma_xx),
		v1 * alpha / math.sqrt(gamma_xx)
	}
end

function ADM1D3Var:eigenLeftTransform(solver, m, v)
	local f = table.unpack(m)
	local _, _, v1, v2, v3 = table.unpack(v)
	return {
		v1 / (2 * f) - v3 / (2 * math.sqrt(f)),
		-2*v1/f + v2,
		v1 / (2 * f) + v3 / (2 * math.sqrt(f))
	}
end

function ADM1D3Var:eigenRightTransform(solver, m, v)
	local f = table.unpack(m)
	local v1, v2, v3 = table.unpack(v)
	return {
		0,
		0,
		(v1 + v3) * f,
		2 * v1 + v2 + 2 * v3,
		math.sqrt(f) * (-v1 + v3)
	}
end


function ADM1D3Var:calcCellMinMaxEigenvalues(sim, i)
	local alpha, gamma_xx = table.unpack(sim.qs[i])
	local f = self.calc.f(alpha)
	return firstAndLast(self:calcEigenvalues(alpha, gamma_xx, f))
end

function ADM1D3Var:sourceTerm(sim, qs)
	local source = sim:newState()
	for i=1,sim.gridsize do
		local alpha, gamma_xx, a_x, D_g, KTilde_xx = table.unpack(qs[i])
		local f = self.calc.f(alpha)
		local dalpha_f = self.calc.dalpha_f(alpha)
		
		source[i][1] = -alpha * alpha * f * KTilde_xx / (gamma_xx * math.sqrt(gamma_xx))
		source[i][2] = -2 * alpha * KTilde_xx / math.sqrt(gamma_xx)
	end
	return source
end

return ADM1D3Var
