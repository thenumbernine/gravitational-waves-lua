--[[
based on the book "Introduction to 3+1 Numerical Relativity" and on the paper "Introduction to Numerical Relativity", both by Alcubierre

A_x = (ln alpha),x = alpha,x / alpha
D_xxx = (ln gamma_xx),x = gamma_xx,x / gamma_xx
KTilde_xx = sqrt(gamma_xx) K_xx

alpha,x = A_x alpha
gamma_xx,x = gamma_xx D_xxx
K_xx = KTilde_xx / sqrt(gamma_xx)

A_x,t + (alpha f K_xx),x = 0
D_xxx,t + (2 alpha K_xx),x = 0
K_xx,t + (alpha A_x / gamma_xx),x = alpha (K_xx^2 - A_x D_xxx / (2 gamma_xx))

...rewritten for KTilde_xx...
A_x,t + (alpha f KTilde_xx / sqrt(gamma_xx)),x = 0
D_xxx,t + (2 alpha KTilde_xx / sqrt(gamma_xx)),x = 0
KTilde_xx,t + (alpha A_x / sqrt(gamma_xx)),x = 0

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


alright this didn't work.
last bet is to use the eigenspace variables and use identity for the eigenvalues
then use the nonlinear transforms to recover state variables

--]]

local class = require 'ext.class'
local Equation = require 'equation'

local ADM1D3Var = class(Equation)
ADM1D3Var.name = 'ADM 1D 3-Var'

ADM1D3Var.numStates = 3

local State = class(Equation.State)

function State:init(...)
	State.super.init(self, ...)
	for i=1,#self do
		self[i].alpha = 0
		self[i].gamma_xx = 0
	end
end

function State:clone()
	local dst = State.super.clone(self)
	for i=1,#self do
		dst[i].alpha = self[i].alpha
		dst[i].gamma_xx = self[i].gamma_xx
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
		c[i].gamma_xx = a[i].gamma_xx + b[i].gamma_xx
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
		c[i].gamma_xx = (type(a) == 'number' and a or a[i].gamma_xx) * (type(b) == 'number' and b or b[i].gamma_xx)
	end
	return c
end

ADM1D3Var.State = State


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
	local alpha = q:_'alpha'
	local gamma_xx = q:_'gamma_xx'
	local A_x = q:_(1)
	local D_xxx = q:_(2)
	local KTilde_xx = q:_(3)
	local K_xx = KTilde_xx / math.sqrt:o(gamma_xx)
	local volume = alpha * math.sqrt:o(gamma_xx)
	ADM1D3Var:buildGraphInfos{
		{alpha = alpha},
		{A_x = A_x},
		{gamma_xx = gamma_xx},
		{D_xxx = D_xxx},
		{K_xx = K_xx},
		{KTildee_xx = KTilde_xx},
		{volume = volume},
	}
end

local function buildField(call)
	return function(self, sim, i, v)
		local v1, v2, v3 = table.unpack(v)
		
		local avgQ = {}
		for j=1,sim.numStates do 
			avgQ[j] = (sim.qs[i-1][j] + sim.qs[i][j]) / 2
		end
		avgQ.alpha = (sim.qs[i-1].alpha + sim.qs[i].alpha) / 2
		avgQ.gamma_xx = (sim.qs[i-1].gamma_xx + sim.qs[i].gamma_xx) / 2
		
		local A_x, D_xxx, KTilde_xx = table.unpack(avgQ)
		local x = sim.ixs[i]
		local alpha = avgQ.alpha
		local gamma_xx = avgQ.gamma_xx
		local f = self.calc.f(alpha)

		return {call(alpha, f, gamma_xx, A_x, D_xxx, KTilde_xx, v1, v2, v3)}
	end
end

ADM1D3Var.fluxTransform = buildField(function(alpha, f, gamma_xx, A_x, D_xxx, KTilde_xx, v1, v2, v3)
	return
		v3 * alpha * f / math.sqrt(gamma_xx),
		v3 * 2 * alpha / math.sqrt(gamma_xx),
		v1 * alpha / math.sqrt(gamma_xx)
end)

ADM1D3Var.applyLeftEigenvectors = buildField(function(alpha, f, gamma_xx, A_x, D_xxx, KTilde_xx, v1, v2, v3)
	return
		v1 / (2 * f) - v3 / (2 * math.sqrt(f)),
		-2*v1/f + v2,
		v1 / (2 * f) + v3 / (2 * math.sqrt(f))
end)

ADM1D3Var.applyRightEigenvectors = buildField(function(alpha, f, gamma_xx, A_x, D_xxx, KTilde_xx, v1, v2, v3)
	return
		(v1 + v3) * f,
		2 * v1 + v2 + 2 * v3,
		math.sqrt(f) * (-v1 + v3)
end)

function ADM1D3Var:initCell(sim,i)
	local x = sim.xs[i]
	local alpha = self.calc.alpha(x)
	local gamma_xx = self.calc.gamma_xx(x)
	local A_x = self.calc.dx_alpha(x) / self.calc.alpha(x)
	local D_xxx = 1/2 * self.calc.dx_gamma_xx(x)
	local K_xx = self.calc.K_xx(x) 
	local KTilde_xx = K_xx / math.sqrt(gamma_xx)
	return {alpha=alpha, gamma_xx=gamma_xx, A_x, D_xxx, KTilde_xx}
end

function ADM1D3Var:calcEigenvalues(alpha, gamma_xx, f)
	local lambda = alpha * math.sqrt(f / gamma_xx)
	return -lambda, 0, lambda
end

function ADM1D3Var:calcInterfaceEigenBasis(sim,i,qL,qR)
	local alpha = (qL.alpha + qR.alpha) / 2
	local gamma_xx = (qL.gamma_xx + qR.gamma_xx) / 2
	local f = self.calc.f(alpha)
	sim.eigenvalues[i] = {self:calcEigenvalues(alpha, gamma_xx, f)}
end

function ADM1D3Var:calcCellMinMaxEigenvalues(sim, i)
	local q = sim.qs[i]
	local alpha = q.alpha
	local gamma_xx = q.gamma_xx
	local f = self.calc.f(alpha)
	return firstAndLast(self:calcEigenvalues(alpha, gamma_xx, f))
end

function ADM1D3Var:sourceTerm(sim, qs)
	local source = sim:newState()
	for i=1,sim.gridsize do
		local A_x, D_xxx, KTilde_xx = table.unpack(qs[i])
		local alpha = qs[i].alpha
		local gamma_xx = qs[i].gamma_xx
		local f = self.calc.f(alpha)
		local dalpha_f = self.calc.dalpha_f(alpha)
		
		source[i].alpha = -alpha * alpha * f * KTilde_xx / (gamma_xx * math.sqrt(gamma_xx))
		source[i].gamma_xx = -2 * alpha * KTilde_xx / math.sqrt(gamma_xx)
	end
	return source
end

return ADM1D3Var
