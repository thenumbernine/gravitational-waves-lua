--[[
same as Z4 except using the 1D ADM's 
	D_g = (ln gamma_xx),x = gamma_xx,x / gamma_xx
	and KTilde_xx = sqrt(gamma_xx) K_xx

alpha,t = -f alpha^2 (tr K - m Theta)
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

so D_g,t = (ln gamma_xx),xt
	= (gamma_xx,x / gamma_xx),t
	= (gamma_xx,t),x / gamma_xx - gamma_xx,x gamma_xx,t / gamma_xx^2
	= -2 (alpha K_xx),x / gamma_xx + 2 D_g alpha KTilde_xx / gamma_xx^1.5
	= -2 (alpha,x K_xx + alpha K_xx,x) / gamma_xx + 2 D_g alpha KTilde_xx / gamma_xx^1.5
	= - 2 alpha (KTilde_xx / sqrt(gamma_xx)),x / gamma_xx 
		- 2 a_x alpha KTilde_xx / gamma_xx^1.5 
		+ 2 D_g alpha KTilde_xx / gamma_xx^1.5
	= - 2 alpha (KTilde_xx,x / sqrt(gamma_xx) - .5 KTilde_xx gamma_xx,x / gamma_xx^1.5) / gamma_xx 
		- 2 a_x alpha KTilde_xx / gamma_xx^1.5 
		+ 2 D_g alpha KTilde_xx / gamma_xx^1.5
	= -2 alpha KTilde_xx,x / gamma_xx^1.5 
		+ (3 D_g - 2 a_x) alpha KTilde_xx / gamma_xx^1.5

and KTilde_xx,t = (sqrt(gamma_xx) K_xx),t
	= .5 gamma_xx,t K_xx / sqrt(gamma_xx) + sqrt(gamma_xx) K_xx,t
	= (-2 alpha K_xx) K_xx / (2 sqrt(gamma_xx)) + sqrt(gamma_xx) (
		- alpha a_x,x 
		+ 2 alpha Z_x,x 
		- alpha a_x^2 
		+ alpha d_xxx a_x / gamma_xx
		- 2 alpha d_xxx Z_x / gamma_xx
		- alpha K_xx^2 / gamma_xx
		- 2 alpha Theta K_xx
		- alpha S_xx
		+ 1/2 alpha (S_xx / gamma_xx - tau) gamma_xx)
	= - alpha sqrt(gamma_xx) a_x,x
		+ 2 alpha sqrt(gamma_xx) Z_x,x
		- 2 alpha KTilde_xx^2 / gamma_xx^1.5
		- alpha a_x^2 sqrt(gamma_xx)
		+ .5 alpha D_g a_x sqrt(gamma_xx)
		- alpha D_g Z_x sqrt(gamma_xx)
		- 2 alpha Theta KTilde_xx / gamma_xx^.5
		- .5 alpha S_xx gamma_xx^.5 
		- .5 alpha tau gamma_xx^1.5

alpha,t = -f alpha^2 (KTilde_xx / gamma_xx^1.5 - m Theta)

gamma_xx,t = -2 alpha KTilde_xx / gamma_xx^.5

a_x,t = - f alpha K_xx,x / gamma_xx
		+ f alpha m Theta,x
		+ f alpha KTilde_xx D_g / gamma_xx^1.5
		- alpha a_x (f' alpha + f) (KTilde_xx / gamma_xx^1.5 - m Theta)
	= -f alpha KTilde_xx,x / gamma_xx^1.5 
		+ f alpha m Theta,x
		+ 1.5 f alpha D_g KTilde_xx / gamma_xx^1.5
		- alpha a_x (f' alpha + f) (KTilde_xx / gamma_xx^1.5 - m Theta)

Theta,t = alpha Z_x,x / gamma_xx
		- 1/2 alpha Z_x D_g / gamma_xx
		- alpha Theta KTilde_xx / gamma_xx^1.5
		- alpha a_x Z_x / gamma_xx
		- alpha tau

Z_x,t = alpha Theta,x
		- 2 alpha Z_x KTilde_xx / gamma_xx^1.5
		- alpha S_x
		- alpha Theta a_x


All together as a matrix:

[alpha		]	[	0	0	0						0	0						0			0						][alpha		]	[ -f alpha^2 (KTilde_xx / gamma_xx^1.5 - m Theta)																							]
[gamma_xx	]	[	0	0	0						0	0						0			0						][gamma_xx	]	[ -2 alpha KTilde_xx / gamma_xx^.5																											]
[a_x		]	[	0	0	0						0	f alpha / gamma_xx^1.5	-f alpha m	0						][a_x		]	[ alpha ((1.5 f D_g - a_x (f' alpha + f)) KTilde_xx / gamma_xx^1.5 + a_x (f' alpha + f) m Theta)											]
[D_g		] + [	0	0	0						0	2 alpha / gamma_xx^1.5	0			0						][D_g		] = [ (3 D_g - 2 a_x) alpha KTilde_xx / gamma_xx^1.5																							]
[KTilde_xx	]	[	0	0	alpha sqrt(gamma_xx)	0	0						0			-2 alpha sqrt(gamma_xx) ][KTilde_xx	]	[ alpha sqrt(gamma_xx) (.5 D_g a_x - a_x^2 - D_g Z_x - 2 KTilde_xx / gamma_xx (KTilde_xx / gamma_xx + Theta) - .5 S_xx - .5 tau gamma_xx)	]
[Theta		]	[	0	0	0						0	0						0			-alpha / gamma_xx		][Theta		]	[ -alpha (Theta KTilde_xx / gamma_xx^1.5 + (a_x + 1/2 D_g) Z_x / gamma_xx + tau)															]
[Z_x		]	[	0	0	0						0	0						-alpha		0						][Z_x		],x	[ -alpha (2 Z_x KTilde_xx / gamma_xx^1.5 + S_x + Theta a_x)																					]

gets results:

eigenvalues:
matrix([-(alpha*sqrt(f))/sqrt(g_xx),0,0,0,0],[0,-alpha/sqrt(g_xx),0,0,0],[0,0,0,0,0],[0,0,0,alpha/sqrt(g_xx),0],[0,0,0,0,(alpha*sqrt(f))/sqrt(g_xx)])

right eigenvectors:
matrix([f,f*m-2*f,0,f*m-2*f,f],[2*g_xx^3,2*f*g_xx^3*m-4*g_xx^3,1,2*f*g_xx^3*m-4*g_xx^3,2*g_xx^3],[-sqrt(f)*g_xx,2*g_xx-f*g_xx*m,0,f*g_xx*m-2*g_xx,sqrt(f)*g_xx],[0,-(f-1)/sqrt(g_xx),0,(f-1)/sqrt(g_xx),0],[0,1-f,0,1-f,0])

left eigenvectors:
matrix([1/(2*f),0,-1/(2*sqrt(f)*g_xx),(f*g_xx*m-2*g_xx)/(sqrt(f)*(2*f-2)*sqrt(g_xx)),(m-2)/(2*f-2)],[0,0,0,-sqrt(g_xx)/(2*f-2),-1/(2*f-2)],[-(2*g_xx^3)/f,1,0,0,2*g_xx^3*m],[0,0,0,sqrt(g_xx)/(2*f-2),-1/(2*f-2)],[1/(2*f),0,1/(2*sqrt(f)*g_xx),-(f*g_xx*m-2*g_xx)/(sqrt(f)*(2*f-2)*sqrt(g_xx)),(m-2)/(2*f-2)])
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
	local KTilde_xx = q:_(5)
	local K_xx = KTilde_xx / math.sqrt:o(gamma_xx)
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
		{KTilde_xx = KTilde_xx},
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
	local KTilde_xx = K_xx / math.sqrt(gamma_xx)
	-- what is Theta and Z_i initialized to?
	local Theta = 0	
	local Z_x = 0
	return {alpha, gamma_xx, a_x, D_g, KTilde_xx, Theta, Z_x}
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
		((m-2)*v5)/(2*f-2)+((f*gamma_xx*m-2*gamma_xx)*v4)/(math.sqrt(f)*(2*f-2)*math.sqrt(gamma_xx))-v3/(2*math.sqrt(f)*gamma_xx)+v1/(2*f),
		-v5/(2*f-2)-(math.sqrt(gamma_xx)*v4)/(2*f-2),
		2*gamma_xx^3*m*v5+v2-(2*gamma_xx^3*v1)/f,
		(math.sqrt(gamma_xx)*v4)/(2*f-2)-v5/(2*f-2),
		((m-2)*v5)/(2*f-2)-((f*gamma_xx*m-2*gamma_xx)*v4)/(math.sqrt(f)*(2*f-2)*math.sqrt(gamma_xx))+v3/(2*math.sqrt(f)*gamma_xx)+v1/(2*f),
	}
end

function Z41Dv2:eigenRightTransform(solver, evR, v)
	local f, gamma_xx = table.unpack(evR)
	local v1, v2, v3, v4, v5 = table.unpack(v)
	return {
		0,
		0,
		f*v5+(f*m-2*f)*v4+(f*m-2*f)*v2+f*v1,
		2*gamma_xx^3*v5+(2*f*gamma_xx^3*m-4*gamma_xx^3)*v4+v3+(2*f*gamma_xx^3*m-4*gamma_xx^3)*v2+2*gamma_xx^3*v1,
		math.sqrt(f)*(gamma_xx*v5-gamma_xx*v1)+(f*gamma_xx*m-2*gamma_xx)*v4+(2*gamma_xx-f*gamma_xx*m)*v2,
		((f-1)*v4+(1-f)*v2)/math.sqrt(gamma_xx),
		(1-f)*v4+(1-f)*v2
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
		local alpha, gamma_xx, a_x, d_xxx, K_xx, Theta, Z_x = table.unpack(qs[i])
		local f = self.calc.f(alpha)
		local dalpha_f = self.calc.dalpha_f(alpha)
		
		source[i][1] = -f * alpha * alpha * (K_xx / gamma_xx - m * Theta)
		source[i][2] = -2 * alpha * K_xx
		--source[i][3] = -alpha * a_x * (dalpha_f * alpha + f) * (K_xx / gamma_xx - m * Theta) + 2 * f * alpha * K_xx * d_xxx / (gamma_xx * gamma_xx)
		--source[i][4] = -alpha * a_x * K_xx
		source[i][5] = alpha * (-a_x * a_x + d_xxx * (a_x - 2 * Z_x) / gamma_xx - K_xx * (K_xx / gamma_xx + 2 * Theta) - 1/2 * (S_xx + tau * gamma_xx))
		--source[i][6] = -alpha * ((Z_x * (a_x + d_xxx / gamma_xx) + Theta * K_xx) / gamma_xx + tau)
		--source[i][7] = -alpha * (2 * Z_x * K_xx / gamma_xx + Theta * a_x + S_x)
	end
	return source
end

return Z41Dv2
