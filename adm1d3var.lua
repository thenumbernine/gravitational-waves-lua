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

require 'ext'
local Simulation = require 'simulation'

local ADM1D3VarSim = class(Simulation)

ADM1D3VarSim.numStates = 3

-- initial conditions
function ADM1D3VarSim:init(args, ...)
	ADM1D3VarSim.super.init(self, args, ...)

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

	local get_state = index:bind(self.qs)
	local get_alpha = get_state:index'alpha'
	local get_g_xx = get_state:index'g_xx'
	local get_A_x = get_state:index(1)
	local get_D_xxx = get_state:index(2)
	local get_KTilde_xx = get_state:index(3)
	local get_K_xx = get_KTilde_xx / sqrt:compose(get_g_xx)
	self.graphInfos = {
		{viewport={0/3, 0/3, 1/3, 1/3}, getter=get_alpha, name='alpha', color={1,0,1}},
		{viewport={0/3, 1/3, 1/3, 1/3}, getter=get_A_x, name='A_x', color={0,1,0}},
		{viewport={1/3, 0/3, 1/3, 1/3}, getter=get_g_xx, name='g_xx', color={.5,.5,1}},
		{viewport={1/3, 1/3, 1/3, 1/3}, getter=get_D_xxx, name='D_xxx', color={1,1,0}},
		{viewport={2/3, 0/3, 1/3, 1/3}, getter=get_K_xx, name='K_xx', color={0,1,1}},
		{viewport={2/3, 1/3, 1/3, 1/3}, getter=get_alpha * sqrt:compose(get_g_xx), name='volume', color={0,1,1}},
		{viewport={0/3, 2/3, 1/3, 1/3}, getter=log:compose(index:bind(self.eigenbasisErrors)), name='log eigenbasis error', color={1,0,0}, range={-30, 30}},
		{viewport={1/3, 2/3, 1/3, 1/3}, getter=log:compose(index:bind(self.fluxMatrixErrors)), name='log reconstuction error', color={1,0,0}, range={-30, 30}},
	}

end

local function buildField(call)
	return function(self, i, v)
		local v1, v2, v3 = unpack(v)
		
		local avgQ = {}
		for j=1,self.numStates do 
			avgQ[j] = (self.qs[i-1][j] + self.qs[i][j]) / 2
		end
		avgQ.alpha = (self.qs[i-1].alpha + self.qs[i].alpha) / 2
		avgQ.g_xx = (self.qs[i-1].g_xx + self.qs[i].g_xx) / 2
		
		local A_x, D_xxx, KTilde_xx = unpack(avgQ)
		local x = self.ixs[i]
		local alpha = avgQ.alpha
		local g_xx = avgQ.g_xx
		local f = self.calc.f(alpha)

		return {call(alpha, f, g_xx, A_x, D_xxx, KTilde_xx, v1, v2, v3)}
	end
end

ADM1D3VarSim.fluxTransform = buildField(function(alpha, f, g_xx, A_x, D_xxx, KTilde_xx, v1, v2, v3)
	return
		v3 * alpha * f / sqrt(g_xx),
		v3 * 2 * alpha / sqrt(g_xx),
		v1 * alpha / sqrt(g_xx)
end)

ADM1D3VarSim.eigenfields = buildField(function(alpha, f, g_xx, A_x, D_xxx, KTilde_xx, v1, v2, v3)
	return
		v1 / (2 * f) - v3 / (2 * sqrt(f)),
		-2*v1/f + v2,
		v1 / (2 * f) + v3 / (2 * sqrt(f))
end)

ADM1D3VarSim.eigenfieldsInverse = buildField(function(alpha, f, g_xx, A_x, D_xxx, KTilde_xx, v1, v2, v3)
	return
		(v1 + v3) * f,
		2 * v1 + v2 + 2 * v3,
		sqrt(f) * (-v1 + v3)
end)

function ADM1D3VarSim:initCell(i)
	local x = self.xs[i]
	local alpha = self.calc.alpha(x)
	local g_xx = self.calc.g_xx(x)
	local A_x = self.calc.dx_alpha(x) / self.calc.alpha(x)
	local D_xxx = 1/2 * self.calc.dx_g_xx(x)
	local K_xx = self.calc.K_xx(x) 
	local KTilde_xx = K_xx / sqrt(g_xx)
	return {alpha=alpha, g_xx=g_xx, A_x, D_xxx, KTilde_xx}
end

function ADM1D3VarSim:calcInterfaceEigenBasis(i)
	local alpha = (self.qs[i-1].alpha + self.qs[i].alpha) / 2
	local g_xx = (self.qs[i-1].g_xx + self.qs[i].g_xx) / 2
	local f = self.calc.f(alpha)
	local lambda = alpha * sqrt(f / g_xx)		
	self.eigenvalues[i] = {-lambda, 0, lambda}
end

function ADM1D3VarSim:zeroDeriv(dq_dts)
	ADM1D3VarSim.super.zeroDeriv(self, dq_dts)
	-- zero deriv
	for i=1,self.gridsize do
		dq_dts[i].alpha = 0
		dq_dts[i].g_xx = 0
	end
end

function ADM1D3VarSim:addSourceToDerivCell(dq_dts, i)
	local A_x, D_xxx, KTilde_xx = unpack(self.qs[i])
	local alpha = self.qs[i].alpha
	local g_xx = self.qs[i].g_xx
	local f = self.calc.f(alpha)
	local dalpha_f = self.calc.dalpha_f(alpha)
	
	dq_dts[i].alpha = dq_dts[i].alpha - alpha * alpha * f * KTilde_xx / (g_xx * sqrt(g_xx))
	dq_dts[i].g_xx = dq_dts[i].g_xx - 2 * alpha * KTilde_xx / sqrt(g_xx)
end

function ADM1D3VarSim:integrateDeriv(dq_dts, dt)
	ADM1D3VarSim.super.integrateDeriv(self, dq_dts, dt)
	for i=1,self.gridsize do
		self.qs[i].alpha = self.qs[i].alpha + dt * dq_dts[i].alpha
		self.qs[i].g_xx = self.qs[i].g_xx + dt * dq_dts[i].g_xx
	end
	self.t = self.t + dt
end

return ADM1D3VarSim

