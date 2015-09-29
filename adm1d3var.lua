--[[
based on the book "Introduction to 3+1 Numerical Relativity" and on the paper "Introduction to Numerical Relativity", both by Alcubierre

A_x = (ln alpha),x = alpha,x / alpha
D_xxx = (ln g_xx),x = g_xx,x / g_xx
KTilde_xx = sqrt(g_xx) K_xx

alpha,x = A_x alpha
g_xx,x = g_xx D_xxx
K_xx = KTilde_xx / sqrt(g_xx)

A_x,t + (alpha f K_xx),x = 0
D_xxx,t + (2 alpha K_xx),x = 0
K_xx,t + (alpha A_x / g_xx),x = alpha (K_xx^2 - A_x D_xxx / (2 g_xx))

...rewritten for KTilde_xx...
A_x,t + (alpha f KTilde_xx / sqrt(g_xx)),x = 0
D_xxx,t + (2 alpha KTilde_xx / sqrt(g_xx)),x = 0
KTilde_xx,t + (alpha A_x / sqrt(g_xx)),x = 0

now we look at eigenvectors ...

/* [wxMaxima: input   start ] */
loag("eigen");
/* [wxMaxima: input   end   ] */
/* [wxMaxima: input   start ] */
assume(alpha>0,f>0);
/* [wxMaxima: input   end   ] */
/* [wxMaxima: input   start ] */
F : matrix(
[0, 0, alpha *f/sqrt(g_xx)],
[0, 0, 2*alpha/sqrt(g_xx)],
[alpha/sqrt(g_xx), 0, 0]
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
ADM1D3Var.name = 'ADM1D3Var'

ADM1D3Var.numStates = 3

-- TODO move this to Equation? 
-- it's intrinsically bound to the Equation: graphInfos, sourceTerm, eigenvectors, and all
local State = class(Equation.State)

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
	local exprs = table{'alpha', 'g_xx', 'K_xx'}:map(function(name)
		return makesym(name), name
	end)

	-- derived functions
	exprs.dx_g_xx = exprs.g_xx:diff(x):simplify()
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

ADM1D3Var.graphInfos = table{
	{viewport={0/3, 0/3, 1/3, 1/3}, getter=function(self,i) return self.qs[i].alpha end, name='alpha', color={1,0,1}},
	{viewport={0/3, 1/3, 1/3, 1/3}, getter=function(self,i) return self.qs[i][1] end, name='A_x', color={0,1,0}},
	{viewport={1/3, 0/3, 1/3, 1/3}, getter=function(self,i) return self.qs[i].g_xx end, name='g_xx', color={.5,.5,1}},
	{viewport={1/3, 1/3, 1/3, 1/3}, getter=function(self,i) return self.qs[i][2] end, name='D_xxx', color={1,1,0}},
	{viewport={2/3, 0/3, 1/3, 1/3}, getter=function(self,i) return self.qs[i][3] / math.sqrt(self.qs[i].g_xx) end, name='K_xx', color={0,1,1}},
	{viewport={2/3, 1/3, 1/3, 1/3}, getter=function(self,i) return self.qs[i].alpha * math.sqrt(self.qs[i].g_xx) end, name='volume', color={0,1,1}},
	{viewport={0/3, 2/3, 1/3, 1/3}, getter=function(self,i) return math.log(self.eigenbasisErrors[i]) end, name='log eigenbasis error', color={1,0,0}, range={-30, 30}},
	{viewport={1/3, 2/3, 1/3, 1/3}, getter=function(self,i) return math.log(self.fluxMatrixErrors[i]) end, name='log reconstuction error', color={1,0,0}, range={-30, 30}},
}
ADM1D3Var.graphInfoForNames = ADM1D3Var.graphInfos:map(function(info,i)
	return info, info.name
end)

local function buildField(call)
	return function(self, sim, i, v)
		local v1, v2, v3 = unpack(v)
		
		local avgQ = {}
		for j=1,sim.numStates do 
			avgQ[j] = (sim.qs[i-1][j] + sim.qs[i][j]) / 2
		end
		avgQ.alpha = (sim.qs[i-1].alpha + sim.qs[i].alpha) / 2
		avgQ.g_xx = (sim.qs[i-1].g_xx + sim.qs[i].g_xx) / 2
		
		local A_x, D_xxx, KTilde_xx = unpack(avgQ)
		local x = sim.ixs[i]
		local alpha = avgQ.alpha
		local g_xx = avgQ.g_xx
		local f = self.calc.f(alpha)

		return {call(alpha, f, g_xx, A_x, D_xxx, KTilde_xx, v1, v2, v3)}
	end
end

ADM1D3Var.fluxTransform = buildField(function(alpha, f, g_xx, A_x, D_xxx, KTilde_xx, v1, v2, v3)
	return
		v3 * alpha * f / sqrt(g_xx),
		v3 * 2 * alpha / sqrt(g_xx),
		v1 * alpha / sqrt(g_xx)
end)

ADM1D3Var.eigenfields = buildField(function(alpha, f, g_xx, A_x, D_xxx, KTilde_xx, v1, v2, v3)
	return
		v1 / (2 * f) - v3 / (2 * sqrt(f)),
		-2*v1/f + v2,
		v1 / (2 * f) + v3 / (2 * sqrt(f))
end)

ADM1D3Var.eigenfieldsInverse = buildField(function(alpha, f, g_xx, A_x, D_xxx, KTilde_xx, v1, v2, v3)
	return
		(v1 + v3) * f,
		2 * v1 + v2 + 2 * v3,
		sqrt(f) * (-v1 + v3)
end)

function ADM1D3Var:initCell(sim,i)
	local x = sim.xs[i]
	local alpha = self.calc.alpha(x)
	local g_xx = self.calc.g_xx(x)
	local A_x = self.calc.dx_alpha(x) / self.calc.alpha(x)
	local D_xxx = 1/2 * self.calc.dx_g_xx(x)
	local K_xx = self.calc.K_xx(x) 
	local KTilde_xx = K_xx / sqrt(g_xx)
	return {alpha=alpha, g_xx=g_xx, A_x, D_xxx, KTilde_xx}
end

function ADM1D3Var:calcInterfaceEigenBasis(sim,i,qL,qR)
	local alpha = (qL.alpha + qR.alpha) / 2
	local g_xx = (qL.g_xx + qR.g_xx) / 2
	local f = self.calc.f(alpha)
	local lambda = alpha * sqrt(f / g_xx)
	sim.eigenvalues[i] = {-lambda, 0, lambda}
end

function ADM1D3Var:sourceTerm(sim)
	local source = sim:newState()
	for i=1,sim.gridsize do
		local A_x, D_xxx, KTilde_xx = unpack(sim.qs[i])
		local alpha = sim.qs[i].alpha
		local g_xx = sim.qs[i].g_xx
		local f = self.calc.f(alpha)
		local dalpha_f = self.calc.dalpha_f(alpha)
		
		source[i].alpha = -alpha * alpha * f * KTilde_xx / (g_xx * sqrt(g_xx))
		source[i].g_xx = -2 * alpha * KTilde_xx / sqrt(g_xx)
	end
	return source
end

return ADM1D3Var
