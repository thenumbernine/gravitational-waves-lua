--[[
based on Alcubierre 1997 "The appearance of coordinate shocks in hyperbolic formalisms in General Relativity" (http://arxiv.org/pdf/gr-qc/9609015v2.pdf)
which itself doesn't specify the eigenvectors, just the equations and eigenfields.
If you follow the book, it explains how to re-cast those same equations as the proper formalism (as I use in adm1d3to5var.lua)

a_x = ln(alpha)_,x
d_xxx = 1/2 gamma_xx,x
D_g = (ln gamma_xx),x = gamma_xx,x / gamma_xx = 2 d_xxx / gamma_xx

hyperbolic formalism:

state vector: [alpha gamma_xx a_x d_xxx K_xx]	<- even though alpha and gamma_xx are already represented by a_x and d_xxx ...
fluxes vector: alpha K_xx * [0, 0, f/gamma_xx, 1, a_x/K_xx]
source vector: alpha / gamma_xx * [-alpha f K_xx, -2 K_xx gamma_xx, 0, 0, a_x d_xxx - K_xx^2]

alpha,t = -alpha^2 f K_xx / gamma_xx
gamma_xx,t = -2 alpha K_xx
a_x,t + (alpha f K_xx / gamma_xx),x = 0
d_xxx,t + (alpha K_xx),x = 0
K_xx,t + (alpha a_x),x = alpha / gamma_xx (a_x d_xxx - K_xx^2)

verify K_xx,t from the K,t equation in "Toy 1+1" chapter of Alcubierre's book
K,t + (alpha a_x / gamma_xx),x = alpha (K^2 - a_x D_g / (2 gamma_xx))
(K_xx / gamma_xx),t + (alpha a_x / gamma_xx),x = alpha (K^2 - a_x D_g / (2 gamma_xx))
K_xx,t / gamma_xx - K_xx gamma_xx,t / gamma_xx^2 + (alpha a_x / gamma_xx),x = alpha (K^2 - a_x D_g / (2 gamma_xx))
K_xx,t + gamma_xx (alpha a_x / gamma_xx),x = alpha (-gamma_xx K^2 - 1/2 a_x D_g)
K_xx,t + gamma_xx (alpha a_x / gamma_xx),x + gamma_xx,x alpha a_x / gamma_xx
	= alpha (-gamma_xx K^2 - 1/2 a_x D_g) + gamma_xx,x alpha a_x / gamma_xx
K_xx,t + (alpha a_x),x = alpha / gamma_xx (a_x d_xxx - K_xx^2)
check.

...linearized

a_x,t + alpha,x f K_xx / gamma_xx
	+ alpha f' alpha,x K_xx / gamma_xx
	+ alpha f K_xx,x / gamma_xx
	- alpha f K_xx gamma_xx,x / gamma_xx^2 = 0
a_x,t + alpha f K_xx,x / gamma_xx
	= ((2 d_xxx / gamma_xx - a_x) f - alpha f' a_x) alpha K_xx / gamma_xx 

d_xxx,t + alpha,x K_xx + alpha K_xx,x = 0
d_xxx,t + alpha K_xx,x = -alpha a_x K_xx 

K_xx,t + alpha,x a_x + alpha a_x,x = alpha / gamma_xx (a_x d_xxx - K_xx^2)
K_xx,t + alpha a_x,x = alpha ((d_xxx / gamma_xx - a_x) a_x - K_xx^2 / gamma_xx)

... as a matrix:
[ a_x ]   [ 0,    0, alpha f / gamma_xx ][ a_x ]   [((2 d_xxx / gamma_xx - a_x) f - alpha f' a_x) alpha K_xx / gamma_xx 
[d_xxx] + [ 0,    0,    alpha           ][d_xxx] = [-alpha a_x K_xx
[ K_xx]   [alpha, 0,      0             ][ K_xx],x [alpha ((d_xxx / gamma_xx - a_x) a_x - K_xx^2 / gamma_xx)


	Maxima to verify linearization and source terms: 
depends([alpha, gamma_xx, a_x, d_xxx, K_xx], x);
depends(f, alpha);
F : transpose(matrix([
	alpha * f * K_xx / gamma_xx,
	alpha * K_xx,
	alpha * a_x
]));
/* here is the flux, differentiated, to start linearization */
diff(F, x)$ ratsimp(%);
/* here is the derivatives substituted for the first-order variables */
subst([diff(alpha,x)=alpha*a_x, diff(gamma_xx,x)=2*d_xxx], %)$ ratsimp(%);
/* here are the derivatives removed, leaving the negative source */
subst([diff(a_x,x)=0,diff(K_xx,x)=0], %)$ ratsimp(%);
/* here I'm subtracing them from the flux source to get the linearized source */
S - %$ ratsimp(%);

	then to find the eigenvectors:


load("eigen");
assume(alpha>0, f>0);
A : matrix([0, 0, alpha*f/gamma_xx], [0, 0, alpha], [alpha, 0, 0]);
results : eigenvectors(A);
lambdas : diag_matrix(results[1][1][1], results[1][1][3], results[1][1][2]);
	matrix(
		[-(alpha*sqrt(f))/sqrt(gamma_xx),0,0],
		[0,0,0],
		[0,0,(alpha*sqrt(f))/sqrt(gamma_xx)]
	)
evR : transpose(matrix(                                                 
     results[2][1][1] * f,
     results[2][3][1],
     results[2][2][1] * f
));
	matrix(
		[f,0,f],
		[gamma_xx,1,gamma_xx],
		[-sqrt(f)*sqrt(gamma_xx),0,sqrt(f)*sqrt(gamma_xx)]
	)
invert(evR)$ ratsimp(%)$ evL : %;
	matrix(
		[1/(2*f),0,-1/(2*sqrt(f)*sqrt(gamma_xx))],
		[-gamma_xx/f,1,0],
		[1/(2*f),0,1/(2*sqrt(f)*sqrt(gamma_xx))]
	)
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
	local a_x = q:_(3)
	local d_xxx = q:_(4)
	local D_g = 2 * d_xxx / gamma_xx
	local K_xx = q:_(5)
	local KTilde = K_xx / math.sqrt:o(gamma_xx)
	local K = K_xx / gamma_xx
	local volume = alpha * math.sqrt:o(gamma_xx)
	ADM1D5Var:buildGraphInfos{
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

function ADM1D5Var:initCell(sim,i)
	local x = sim.xs[i]
	local alpha = self.calc.alpha(x)
	local gamma_xx = self.calc.gamma_xx(x)
	local a_x = self.calc.dx_alpha(x) / self.calc.alpha(x)
	local d_xxx = 1/2 * self.calc.dx_gamma_xx(x)
	local K_xx = self.calc.K_xx(x)
	return {alpha, gamma_xx, a_x, d_xxx, K_xx}
end

function ADM1D5Var:calcRoeValues(qL, qR)
	local alpha, gamma_xx, a_x, d_xxx, K_xx = ADM1D5Var.super.calcRoeValues(self, qL, qR)
	local f = self.calc.f(alpha)
	return alpha, gamma_xx, a_x, d_xxx, K_xx, f
end

function ADM1D5Var:fluxMatrixTransform(solver, m, v)
	local alpha, gamma_xx, a_x, d_xxx, K_xx, f = table.unpack(m)
	local v1, v2, v3, v4, v5 = table.unpack(v)
	-- the alpha and gamma_xx terms are neglected when reconstructing the flux
	-- ... because the eigenvalue is zero?
	return {
	--[[ i think these were favoring derivatives over source terms:
		v1*f*K_xx/gamma_xx - v2*alpha*f*K_xx/gamma_xx^2 + v5*alpha*f/gamma_xx,
		v1*K_xx + v5*alpha,
		v1*a_x + v3*alpha
	--]]
		0,
		0,
		v5*alpha*f/gamma_xx,
		v5*alpha,
		v3*alpha
	}
end

--[[ fixme
function ADM1D5Var:eigenLeftTransform(solver, m, v)
	local alpha, gamma_xx, a_x, d_xxx, K_xx, f = table.unpack(m)
	local v1, v2, v3, v4, v5 = table.unpack(v)
	return {
		v1 * (gamma_xx * a_x / f - K_xx * math.sqrt(gamma_xx / f)) / (2 * alpha) + v3 * gamma_xx / (2 * f) - v5 * .5 * math.sqrt(gamma_xx / f),
		v1 / alpha,
		-v1 * (gamma_xx * a_x) / (alpha * f) - v3 * gamma_xx / f + v4,
		v2,
		v1 * (gamma_xx * a_x / f + K_xx * math.sqrt(gamma_xx / f)) / (2 * alpha) + v3 * gamma_xx / (2 * f) + v5 * .5 * math.sqrt(gamma_xx / f)
	}
end
--]]

function ADM1D5Var:calcMaxEigenvalue(alpha, gamma_xx)
	local f = self.calc.f(alpha)
	local lambda = alpha * math.sqrt(f / gamma_xx)
	return lambda
end

function ADM1D5Var:calcEigenvaluesFromCons(alpha, gamma_xx, a_x, d_xxx, K_xx, f)
	local lambda = self:calcMaxEigenvalue(alpha, gamma_xx)
	return -lambda, 0, 0, 0, lambda
end

function ADM1D5Var:calcEigenBasis(lambdas, evr, evl, dF_dU, alpha, gamma_xx, a_x, d_xxx, K_xx, f)
	local sqrt_f = math.sqrt(f)
	local sqrt_g = math.sqrt(gamma_xx)
	local lambda = alpha * sqrt_f / sqrt_g 
	fill(lambdas, -lambda, 0, 0, 0, lambda)
	
	-- row-major, math-indexed
	if dF_dU then
		fill(dF_dU, alpha, gamma_xx, a_x, d_xxx, K_xx, f)
		--[[
		fill(dF_dU,
			{0,0,0,0,0},
			{0,0,0,0,0},
			{f*K_xx/gamma_xx, -alpha*f*K_xx/gamma_xx^2, 0,0, alpha*f/gamma_xx},
			{K_xx,0,0,0,alpha},
			{a_x,0,alpha,0,0}
		)
		--]]
	end
	--[[ where did I get this from? probably eigenvector decomposition on the flux
	fill(evr,	
		{0,			alpha,	0,	0,	0			},	-- alpha
		{0,			0,		0,	1,	0			},	-- gamma_xx
		{f/gamma_xx,		-a_x,		0,	0,	f/gamma_xx			},	-- a_x
		{1,			0,		1,	0,	1			},	-- d_xxx
		{-sqrt_f/sqrt_g,-K_xx,		0,	0,	sqrt_f/sqrt_g	}	-- K_xx
	)
	fill(evl,
		{(gamma_xx * a_x / f - K_xx * math.sqrt(gamma_xx / f)) / (2 * alpha), 0, gamma_xx / (2 * f), 0, -.5 * math.sqrt(gamma_xx / f)}, 
		{1 / alpha, 0, 0, 0, 0}, 
		{-(gamma_xx * a_x) / (alpha * f), 0, -gamma_xx / f, 1, 0}, 
		{0, 1, 0, 0, 0}, 
		{(gamma_xx * a_x / f + K_xx * math.sqrt(gamma_xx / f)) / (2 * alpha), 0, gamma_xx / (2 * f), 0, .5 * math.sqrt(gamma_xx / f)} 
	)
	-- note that because we have zero eigenvalues that the eigendecomposition cannot reconstruct the flux matrix
	--]]
	-- [[ here's from the left eigenvectors
	fill(evl,
		{0, 0, -1/sqrt_g, 0, sqrt_f/gamma_xx},	-- math.sqrt(f) K_xx / gamma_xx - a_x / math.sqrt(gamma_xx)
		{1, 0, 0, 0, 0},								-- alpha
		{0, 1, 0, 0, 0},								-- gamma_xx
		{0, 0, 1, -f/gamma_xx, 0},						-- a_x - f d_xxx / gamma_xx
		{0, 0, 1/sqrt_g, 0, sqrt_f/gamma_xx}	-- math.sqrt(f) K_xx / gamma_xx + a_x / math.sqrt(gamma_xx)
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
		local alpha, gamma_xx, a_x, d_xxx, K_xx = unpack(qs[i])
		local f = self.calc.f(alpha)
		local dalpha_f = self.calc.dalpha_f(alpha)
		local K = K_xx / gamma_xx

		source[i][1] = -alpha * alpha * f * K
		source[i][2] = -2 * alpha * K_xx
		source[i][5] = alpha / gamma_xx * (a_x * d_xxx - K_xx * K_xx)

--[[
		source[i][3] = ((2 * d_xxx / gamma_xx - a_x) * f - alpha * dalpha_f * a_x) * alpha * K
		source[i][4] = -alpha * a_x * K_xx
		source[i][5] = alpha * ((d_xxx / gamma_xx - a_x) * a_x - K_xx * K_xx / gamma_xx)
--]]	
	end
	return source
end

return ADM1D5Var
