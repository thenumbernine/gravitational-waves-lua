--[[
based on http://arxiv.org/pdf/gr-qc/9609015v2.pdf
which itself doesn't specify the formalism, just the equations.
If you follow the book, it explains how to re-cast those same equations as the proper formalism (as I use in adm1d3to5var.lua)


hyperbolic formalism:

state vector: [alpha g_xx A_x D_xxx K_xx]	<- even though alpha and g_xx are already represented by A_x and D_xxx ...
fluxes vector: alpha K_xx * [0, 0, f/g_xx, 1, A_x/K_xx]
source vector: alpha / g_xx * [-alpha f K_xx, -2 K_xx g_xx, 0, 0, A_x D_xxx - K_xx^2]

alpha,t + (0),x = -alpha^2 f K_xx / g_xx
g_xx,t + (0),x = -2 alpha K_xx
A_x,t + (alpha K_xx f / g_xx),x = 0
D_xxx,t + (alpha K_xx),x = 0
K_xx,t + (alpha A_x),x = alpha / g_xx (A_x D_xxx - K_xx^2)

from here on I put the df/dalpha term in the source

eigenvalues of A_x:
-lambda * (
	-lambda * (
		-lambda * det(
			-lambda		alpha f / g_xx
			alpha		-lamda
		)
	)
) = 0
<=>
lambda = 0 has multiplicity 3,

lambda^2 - alpha^2 f / g_xx = 0 
lambda = +-alpha sqrt(f / g_xx)

for lambda = -alpha sqrt(f/g_xx) the eigenvector is 
	[0, 0, 1, g_xx/f, -sqrt(g_xx/f)] dot state = A_x + D_xxx g_xx/f - K_xx sqrt(g_xx/f)
	= A_x + D_xxx f / g_xx - K_xx * sqrt(f / g_xx)
for lambda = +alpha sqrt(f/g_xx) the eigenvector is 
	[0, 0, 1, g_xx/f, +sqrt(g_xx/f)] dot state = A_x + D_xxx g_xx/f + K_xx sqrt(g_xx/f)
for lambda = 0 the eigenvectors are
[alpha,0,-A_x,0,-K_xx]
[0,0,0,1,0]
[0,1,0,0,0]

eigenvector matrix is :
[	0		  		alpha	0	0	0				]
[	0		  		0		0	1	0				]
[	f/g_xx			-A_x	0	0	f/g_xx			]
[	1				0		1	0	1				]
[	-sqrt(f/g_xx)	-K_xx	0	0	sqrt(f/g_xx)	]
inverse:
[(sqrt(f)*g_xx^(3/2)*A_x-f*g_xx*K_xx)/(2*alpha*f^(3/2)*sqrt(g_xx)),0,g_xx/(2*f),0,-sqrt(g_xx)/(2*sqrt(f))]
[1/alpha,0,0,0,0]
[-(g_xx*A_x)/(alpha*f),0,-g_xx/f,1,0]
[0,1,0,0,0]
[(sqrt(f)*sqrt(g_xx)*K_xx+g_xx*A_x)/(2*alpha*f),0,g_xx/(2*f),0,sqrt(g_xx)/(2*sqrt(f))]

/* [wxMaxima: input   start ] */
fluxMatrix : matrix(
[0,0,0,0,0],
[0,0,0,0,0],
[f*K_xx/g_xx, -alpha*f*K_xx/g_xx^2, 0,0, alpha*f/g_xx],
[K_xx,0,0,0,alpha],
[A_x,0,alpha,0,0]);
/* [wxMaxima: input   end   ] */
/* [wxMaxima: input   start ] */
load("eigen");
/* [wxMaxima: input   end   ] */
/* [wxMaxima: input   start ] */
assume(f>0,g_xx>0,alpha>0);
/* [wxMaxima: input   end   ] */
/* [wxMaxima: input   start ] */
result:ev(similaritytransform(fluxMatrix), hermitianmatrix=true);
/* [wxMaxima: input   end   ] */
/* [wxMaxima: input   start ] */
eigenvalueMatrix : diag_matrix(
result[1][1][1],
result[1][1][3], 
result[1][1][3],
result[1][1][3],
result[1][1][2]);
/* [wxMaxima: input   end   ] */
/* [wxMaxima: input   start ] */
transpose(matrix(
result[2][1][1] * sqrt(f^2 + f*g_xx + g_xx^2)/g_xx,
result[2][3][1] * sqrt(K_xx^2 + A_x^2 + alpha^2),
result[2][3][2],
[0,1,0,0,0],
result[2][2][1] * sqrt(f^2 + f*g_xx + g_xx^2)/g_xx))$
ratsimp(%)$
eigenvectorMatrix : %;
/* [wxMaxima: input   end   ] */
/* [wxMaxima: input   start ] */
invert(eigenvectorMatrix)$
ratsimp(%)$
eigenvectorInverseMatrix : %;
/* [wxMaxima: input   end   ] */
/* [wxMaxima: input   start ] */
eigenvectorMatrix . eigenvectorInverseMatrix$
ratsimp(%);
/* [wxMaxima: input   end   ] */
/* [wxMaxima: input   start ] */
eigenvectorMatrix . eigenvalueMatrix . eigenvectorInverseMatrix $
ratsimp(%);
is(% = fluxMatrix);
fluxMatrix;
/* [wxMaxima: input   end   ] */

--]]

local class = require 'ext.class'
local table = require 'ext.table'

local ADM1D5Var = class()
ADM1D5Var.name = 'ADM1D5Var'

ADM1D5Var.numStates = 5 

function ADM1D5Var:init(args, ...)

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
	exprs.dx_alpha = exprs.alpha:diff(x):simplify()
	exprs.dx_g_xx = exprs.g_xx:diff(x):simplify()

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

ADM1D5Var.graphInfos = table{
	{viewport={0/3, 0/3, 1/3, 1/3}, getter=function(self,i) return self.qs[i][1] end, name='alpha', color={1,0,1}},
	{viewport={0/3, 1/3, 1/3, 1/3}, getter=function(self,i) return self.qs[i][3] end, name='A_x', color={0,1,0}},
	{viewport={1/3, 0/3, 1/3, 1/3}, getter=function(self,i) return self.qs[i][2] end, name='g_xx', color={.5,.5,1}},
	{viewport={1/3, 1/3, 1/3, 1/3}, getter=function(self,i) return self.qs[i][4] end, name='D_xxx', color={1,1,0}},
	{viewport={2/3, 0/3, 1/3, 1/3}, getter=function(self,i) return self.qs[i][5] end, name='K_xx', color={0,1,1}},
	{viewport={2/3, 1/3, 1/3, 1/3}, getter=function(self,i) return self.qs[i][1] * sqrt(self.qs[i][2]) end, name='volume', color={0,1,1}},
	{viewport={0/3, 2/3, 1/3, 1/3}, getter=function(self,i) return math.log(self.eigenbasisErrors[i]) end, name='log eigenbasis error', color={1,0,0}, range={-30, 30}},
	{viewport={1/3, 2/3, 1/3, 1/3}, getter=function(self,i) return math.log(self.fluxMatrixErrors[i]) end, name='log reconstruction error', color={1,0,0}, range={-30, 30}},
}
ADM1D5Var.graphInfoForNames = ADM1D5Var.graphInfos:map(function(info,i)
	return info, info.name
end)

function ADM1D5Var:initCell(sim,i)
	local x = sim.xs[i]
	local alpha = self.calc.alpha(x)
	local g_xx = self.calc.g_xx(x)
	local A_x = self.calc.dx_alpha(x) / self.calc.alpha(x)
	local D_xxx = 1/2 * self.calc.dx_g_xx(x)
	local K_xx = self.calc.K_xx(x)
	return {alpha, g_xx, A_x, D_xxx, K_xx}
end

function ADM1D5Var:calcInterfaceEigenBasis(sim,i,qL,qR)
	local avgQ = {}
	for j=1,self.numStates do
		avgQ[j] = (qL[j] + qR[j]) / 2
	end
	
	local alpha, g_xx, A_x, D_xxx, K_xx = unpack(avgQ)
	local x = sim.ixs[i]
	local f = self.calc.f(alpha)
	local lambda = alpha * sqrt(f / g_xx)		
	sim.eigenvalues[i] = {-lambda, 0, 0, 0, lambda}
	-- row-major, math-indexed
	sim.fluxMatrix[i] = {
		{0,0,0,0,0},
		{0,0,0,0,0},
		{f*K_xx/g_xx, -alpha*f*K_xx/g_xx^2, 0,0, alpha*f/g_xx},
		{K_xx,0,0,0,alpha},
		{A_x,0,alpha,0,0},
	}
	sim.eigenvectors[i] = {
		{0,			alpha,	0,	0,	0			},	-- alpha
		{0,			0,		0,	1,	0			},	-- g_xx
		{f/g_xx,		-A_x,		0,	0,	f/g_xx			},	-- A_x
		{1,			0,		1,	0,	1			},	-- D_xxx
		{-sqrt(f/g_xx),-K_xx,		0,	0,	sqrt(f/g_xx)	},	-- K_xx
	}
	sim.eigenvectorsInverse[i] = {
		{(g_xx * A_x / f - K_xx * sqrt(g_xx / f)) / (2 * alpha), 0, g_xx / (2 * f), 0, -.5 * sqrt(g_xx / f)}, 
		{1 / alpha, 0, 0, 0, 0}, 
		{-(g_xx * A_x) / (alpha * f), 0, -g_xx / f, 1, 0}, 
		{0, 1, 0, 0, 0}, 
		{(g_xx * A_x / f + K_xx * sqrt(g_xx / f)) / (2 * alpha), 0, g_xx / (2 * f), 0, .5 * sqrt(g_xx / f)}, 
	}
	-- note that because we have zero eigenvalues that the eigendecomposition cannot reconstruct the flux matrix
end	
	
return ADM1D5Var

