--[[
same as Z4 except using the 1D ADM's 
	D_g = (ln gamma_xx),x = gamma_xx,x / gamma_xx = 2 d_xxx / gamma_xx
	and KTilde = K_xx / sqrt(gamma_xx) = sqrt(gamma_xx) K

alpha,t = -f alpha^2 (K - m Theta)
	= -f alpha^2 (K_xx / gamma_xx - m Theta)

gamma_xx,t = -2 alpha K_xx

a_x,t = - f alpha K_xx,x / gamma_xx 
		+ f alpha m Theta,x 
		- alpha a_x (f' alpha + f) (K_xx / gamma_xx - m Theta)
		+ 2 f alpha K_xx d_xxx / gamma_xx^2

d_xxx,t = -alpha K_xx,x - alpha a_x K_xx

K_xx,t = - alpha a_x,x 
		+ 2 alpha Z_x,x 
		- alpha a_x^2 
		+ alpha d_xxx a_x / gamma_xx
		- 2 alpha d_xxx Z_x / gamma_xx
		- alpha K_xx^2 / gamma_xx
		- 2 alpha Theta K_xx
		- alpha S_xx
		+ 1/2 alpha (S_xx / gamma_xx - tau) gamma_xx

Theta,t = alpha Z_x,x / gamma_xx
		- alpha Z_x d_xxx / gamma_xx^2
		- alpha Theta K_xx / gamma_xx
		- alpha tau
		- alpha a_x Z_x / gamma_xx

Z_x,t = alpha Theta,x
		- 2 alpha Z_x K_xx / gamma_xx
		- alpha S_x
		- alpha Theta a_x

substitue vars...

alpha,t = -f alpha^2 (K - m Theta)
gamma_xx,t = -2 alpha K_xx
a_x,t + f alpha KTilde,x / sqrt(gamma_xx) 
	- f alpha m Theta,x 
	= alpha (1/2 f K D_g - a_x (f' alpha + f) (K - m Theta))
D_g,t = (ln gamma_xx),xt
	= (gamma_xx,x / gamma_xx),t
	= gamma_xx,tx / gamma_xx - gamma_xx,x gamma_xx,t / gamma_xx^2
	= (-2 alpha K_xx),x / gamma_xx - gamma_xx,x (-2 alpha K_xx) / gamma_xx^2
	= -2 alpha (KTilde sqrt(gamma_xx)),x / gamma_xx + 2 alpha K D_g - 2 alpha a_x K 
	= -2 alpha KTilde,x / sqrt(gamma_xx) + alpha K (D_g - 2 a_x)
D_g,t + 2 alpha KTilde,x / sqrt(gamma_xx) = alpha K (D_g - 2 a_x)
KTilde,t = (K_xx / sqrt(gamma_xx)),t
	= K_xx,t / sqrt(gamma_xx) - 1/2 K_xx gamma_xx,t / gamma_xx^(3/2)
	= (- alpha a_x,x 
		+ 2 alpha Z_x,x 
		- alpha a_x^2 
		+ alpha d_xxx a_x / gamma_xx
		- 2 alpha d_xxx Z_x / gamma_xx
		- alpha K_xx^2 / gamma_xx
		- 2 alpha Theta K_xx
		- alpha S_xx
		+ 1/2 alpha (S_xx / gamma_xx - tau) gamma_xx) / sqrt(gamma_xx)
		- 1/2 K_xx (-2 alpha K_xx) / gamma_xx^(3/2)
KTilde,t + alpha a_x,x / sqrt(gamma_xx)
	- 2 alpha Z_x,x / sqrt(gamma_xx)
	= alpha (
		(	1/2 D_g a_x
			- a_x^2
			- D_g Z_x
			- 1/2 S_xx
			- 1/2 tau gamma_xx
		) / sqrt(gamma_xx)
		- 2 Theta KTilde)

Theta,t - alpha Z_x,x / gamma_xx = -alpha ((1/2 D_g + a_x) Z_x / gamma_xx + Theta K + tau)

Z_x,t - alpha Theta,x = -alpha (2 Z_x K + S_x + Theta a_x)



All together as a matrix:

[alpha		]	[	0	0	0						0	0							0			0							][alpha		]	[ -f alpha^2 (K - m Theta)
[gamma_xx	]	[	0	0	0						0	0							0			0							][gamma_xx	]	[ -2 alpha K_xx
[a_x		]	[	0	0	0						0	f alpha / sqrt(gamma_xx)	-f alpha m	0							][a_x		]	[ alpha (1/2 f K D_g - a_x (f' alpha + f) (K - m Theta))
[D_g		] + [	0	0	0						0	2 alpha / sqrt(gamma_xx)	0			0							][D_g		] = [ alpha K (D_g - 2 a_x)
[KTilde		]	[	0	0	alpha / sqrt(gamma_xx)	0	0							0			-2 alpha / sqrt(gamma_xx)	][KTilde	]	[ alpha ( (1/2 D_g - a_x) a_x - D_g Z_x - 1/2 S_xx - 1/2 tau gamma_xx - 2 Theta K_xx) / sqrt(gamma_xx)
[Theta		]	[	0	0	0						0	0							0			-alpha / gamma_xx			][Theta		]	[ -alpha ((1/2 D_g + a_x) Z_x / gamma_xx + Theta K + tau)
[Z_x		]	[	0	0	0						0	0							-alpha		0							][Z_x		],x	[ -alpha (2 Z_x K + S_x + Theta a_x)


load("eigen");
assume(alpha>0, f>0, gamma_xx>0);
A : matrix([0, 0, f * alpha / sqrt(gamma_xx), - f * alpha * m, 0],
[0, 0, 2 * alpha / sqrt(gamma_xx), 0, 0],
[alpha / sqrt(gamma_xx), 0, 0, 0, -2 * alpha / sqrt(gamma_xx)],
[0, 0, 0, 0, -alpha / gamma_xx],
[0, 0, 0, -alpha, 0]);
results : eigenvectors(A);
lambdas : diag_matrix(results[1][1][1], results[1][1][3], results[1][1][5], results[1][1][4], results[1][1][2]);
transpose(matrix(
        results[2][1][1]*f,
        results[2][3][1]*f*(m-2),
        results[2][5][1],
        results[2][4][1]*f*(m-2),
        results[2][2][1]*f))$
ratsimp(%)$
evR : %;
evR$ invert(%)$ ratsimp(%)$ evL:%;
evR . lambdas . evL - A$ ratsimp(%);

eigenvalues:
matrix(
	[-(alpha*sqrt(f))/sqrt(gamma_xx),0,0,0,0],
	[0,-alpha/sqrt(gamma_xx),0,0,0],
	[0,0,0,0,0],
	[0,0,0,alpha/sqrt(gamma_xx),0],
	[0,0,0,0,(alpha*sqrt(f))/sqrt(gamma_xx)])

right eigenvectors:
matrix(
	[f,f*m-2*f,0,f*m-2*f,f],
	[2,2*f*m-4,1,2*f*m-4,2],
	[-sqrt(f),2-f*m,0,f*m-2,sqrt(f)],
	[0,-(f-1)/sqrt(gamma_xx),0,(f-1)/sqrt(gamma_xx),0],
	[0,1-f,0,1-f,0])

left eigenvectors:
matrix(
	[1/(2*f),0,-1/(2*sqrt(f)),(sqrt(gamma_xx)*(f*m-2))/(sqrt(f)*(2*f-2)),(m-2)/(2*f-2)],
	[0,0,0,-sqrt(gamma_xx)/(2*f-2),-1/(2*f-2)],
	[-2/f,1,0,0,2*m],
	[0,0,0,sqrt(gamma_xx)/(2*f-2),-1/(2*f-2)],
	[1/(2*f),0,1/(2*sqrt(f)),-(sqrt(gamma_xx)*(f*m-2))/(sqrt(f)*(2*f-2)),(m-2)/(2*f-2)])

--]]


local class = require 'ext.class'
local Equation = require 'equation'

local Z41Dv2 = class(Equation)
Z41Dv2.name = 'Z4-1D v2'

Z41Dv2.numStates = 7
Z41Dv2.numWaves = 5	-- no waves for alpha and gamma_xx, which are purely source-driven

local m = 2	-- "for f=1 one must have m=2"
local tau = 0
local S_x = 0
local S_xx = 0

-- initial conditions
function Z41Dv2:init(args, ...)

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
	local d_xxx = gamma_xx * D_g / 2
	local KTilde = q:_(5)
	local K_xx = KTilde * math.sqrt:o(gamma_xx)
	local K = K_xx / gamma_xx
	local volume = alpha * math.sqrt:o(gamma_xx)
	local Theta = q:_(6)
	local Z_x = q:_(7)
	Z41Dv2:buildGraphInfos{
		{alpha = alpha},
		{a_x = a_x},
		{gamma_xx = gamma_xx},
		{d_xxx = d_xxx},
		{D_g = D_g},
		{K_xx = K_xx},
		{KTilde = KTilde},
		{K = K},
		{Theta = Theta},
		{Z_x = Z_x},
		{volume = volume},
	}
end

function Z41Dv2:initCell(sim,i)
	local x = sim.xs[i]
	local alpha = self.calc.alpha(x)
	local gamma_xx = self.calc.gamma_xx(x)
	local a_x = self.calc.dx_alpha(x) / self.calc.alpha(x)
	local D_g = self.calc.dx_gamma_xx(x) / self.calc.gamma_xx(x)
	local K_xx = self.calc.K_xx(x) 
	local KTilde = K_xx / math.sqrt(gamma_xx)
	-- what is Theta and Z_i initialized to?
	local Theta = 0	
	local Z_x = 0
	return {alpha, gamma_xx, a_x, D_g, KTilde, Theta, Z_x}
end

function Z41Dv2:calcEigenvalues(alpha, gamma_xx, f)
	local f = self.calc.f(alpha)
	local _1_sqrt_gamma_xx = 1 / math.sqrt(gamma_xx)
	local lambda1 = alpha * math.sqrt(f) * _1_sqrt_gamma_xx
	local lambda2 = alpha * _1_sqrt_gamma_xx
	return -lambda1, -lambda2, 0, lambda2, lambda1
end

-- returns averaging of variables used for interface eigenbasis
function Z41Dv2:calcRoeValues(qL, qR)
	local alpha = (qL[1] + qR[1]) / 2
	local gamma_xx = (qL[2] + qR[2]) / 2
	local f = self.calc.f(alpha)
	return alpha, gamma_xx, f	
end

function Z41Dv2:calcEigenBasis(lambda, evr, evl, dF_dU, alpha, gamma_xx, f)
	fill(lambda, self:calcEigenvalues(alpha, gamma_xx, f))
	fill(evl, f, gamma_xx)
	fill(evr, f, gamma_xx)
	if dF_dU then fill(dF_dU, alpha, gamma_xx, f) end
end

function Z41Dv2:fluxMatrixTransform(solver, A, v)
	local alpha, gamma_xx, f = table.unpack(A)
	local _, _, v1, v2, v3, v4, v5 = table.unpack(v)	-- skip alpha, gamma_xx, and transform the rest: a_x, d_xxx, K_xx, Theta, Z_x
	return {
		0,
		0,
		(alpha*f*v3)/gamma_xx^(3/2)-alpha*f*m*v4,
		2*alpha*gamma_xx^(3/2)*v3,
		alpha*math.sqrt(gamma_xx)*v1-2*alpha*math.sqrt(gamma_xx)*v5,
		-(alpha*v5)/gamma_xx,
		-alpha*v4,
	}
end

function Z41Dv2:eigenLeftTransform(solver, evL, v)
	local f, gamma_xx = table.unpack(evL)
	local _, _, v1, v2, v3, v4, v5 = table.unpack(v)	-- skip alpha, gamma_xx, and transform the rest: a_x, d_xxx, K_xx, Theta, Z_x
	return {
		((m-2)*v5)/(2*f-2)+(math.sqrt(gamma_xx)*(f*m-2)*v4)/(math.sqrt(f)*(2*f-2))-v3/(2*math.sqrt(f))+v1/(2*f),
		-v5/(2*f-2)-(math.sqrt(gamma_xx)*v4)/(2*f-2),
		2*m*v5+v2-(2*v1)/f,
		(math.sqrt(gamma_xx)*v4)/(2*f-2)-v5/(2*f-2),
		((m-2)*v5)/(2*f-2)-(math.sqrt(gamma_xx)*(f*m-2)*v4)/(math.sqrt(f)*(2*f-2))+v3/(2*math.sqrt(f))+v1/(2*f)
	}
end

function Z41Dv2:eigenRightTransform(solver, evR, v)
	local f, gamma_xx = table.unpack(evR)
	local v1, v2, v3, v4, v5 = table.unpack(v)
	return {
		0,
		0,
		f*v5+(f*m-2*f)*v4+(f*m-2*f)*v2+f*v1,
		2*v5+(2*f*m-4)*v4+v3+(2*f*m-4)*v2+2*v1,
		math.sqrt(f)*v5+(f*m-2)*v4+(2-f*m)*v2-math.sqrt(f)*v1,
		((f-1)*v4)/math.sqrt(gamma_xx)+((1-f)*v2)/math.sqrt(gamma_xx),
		(1-f)*v4+(1-f)*v2,
	}
end

function Z41Dv2:calcCellMinMaxEigenvalues(sim, i)
	local alpha, gamma_xx = table.unpack(sim.qs[i])
	local f = self.calc.f(alpha)
	return firstAndLast(self:calcEigenvalues(alpha, gamma_xx, f))
end

function Z41Dv2:sourceTerm(sim, qs)
	local source = sim:newState()
	for i=1,sim.gridsize do
		local alpha, gamma_xx, a_x, D_g, KTilde, Theta, Z_x = table.unpack(qs[i])
		local f = self.calc.f(alpha)
		local dalpha_f = self.calc.dalpha_f(alpha)
		local sqrt_gamma_xx = math.sqrt(gamma_xx)
		local K_xx = KTilde * sqrt_gamma_xx 
		local K = K_xx / gamma_xx

		source[i][1] = -f * alpha^2 * (K - m * Theta)
		source[i][2] = -2 * alpha * K_xx
		source[i][3] = alpha * (.5 * f * K * D_g - a_x * (dalpha_f * alpha + f) * (K - m * Theta))
		source[i][4] = alpha * K * (D_g - 2 * a_x)
		source[i][5] = alpha * ((.5 * D_g - a_x) * a_x - D_g * Z_x - .5 * S_xx - .5 * tau * gamma_xx - 2 * Theta * K_xx) / sqrt_gamma_xx
		source[i][6] = -alpha * ((.5 * D_g + a_x) * Z_x / gamma_xx + Theta * K + tau)
		source[i][7] = -alpha * (2 * Z_x * K + S_x + Theta * a_x)
	end
	return source
end

return Z41Dv2
