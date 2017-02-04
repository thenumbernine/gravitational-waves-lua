--[[
http://arxiv.org/pdf/1106.2254v2.pdf

Shift-less (no beta's or b's)
Is this system weakly hyperbolic?
Should I be using D_g = (ln gamma_xx),x and KTilde_xx = sqrt(gamma_xx) K_xx as with ADM?

here's our variables:
alpha
gamma_xx
a_x = (ln alpha),x
d_xxx = 1/2 gamma_xx,x
K_xx
Theta
Z_x

derived values:
tr K = K_ij gamma^ij = K_xx / gamma_xx

here's our evolution:
alpha,t = -f alpha^2 (tr K - m Theta)
	= -f alpha^2 (K_xx / gamma_xx - m Theta)
gamma_xx,t = -2 alpha K_xx
a_x,t = (ln alpha),xt = (ln alpha),tx
	= (alpha,t / alpha),x
	= (-f alpha (K_xx / gamma_xx - m Theta)),x
	= -alpha,x (f' alpha + f) (K_xx / gamma_xx - m Theta)
		- f alpha (K_xx,x / gamma_xx - K_xx gamma_xx,x / gamma_xx^2 - m Theta,x)
	= 
		- f alpha K_xx,x / gamma_xx 
		+ f alpha m Theta,x 
		- alpha a_x (f' alpha + f) (K_xx / gamma_xx - m Theta)
		+ 2 f alpha K_xx d_xxx / gamma_xx^2
d_xxx,t = -alpha K_xx,x - alpha a_x K_xx
K_xx,t = -alpha a_x,x + 2 alpha Z_x,x - alpha a_x a_x + alpha d_xxx a_x / gamma_xx
		- alpha d_xxx d_xxx / gamma_xx
		- alpha d_xxx d_xxx / gamma_xx 
		+ alpha d_xxx d_xxx / gamma_xx
		+ alpha d_xxx d_xxx / gamma_xx
		- 2 alpha d_xxx Z_x / gamma_xx
		- 2 alpha K_xx K_xx / gamma_xx
		+ alpha (K_xx / gamma_xx - 2 Theta) K_xx
		- alpha S_xx
		+ 1/2 alpha (S_xx / gamma_xx - tau) gamma_xx
	= 
		- alpha a_x,x 
		+ 2 alpha Z_x,x 
		- alpha a_x^2 
		+ alpha d_xxx a_x / gamma_xx
		- 2 alpha d_xxx Z_x / gamma_xx
		- alpha K_xx^2 / gamma_xx
		- 2 alpha Theta K_xx
		- alpha S_xx
		+ 1/2 alpha (S_xx / gamma_xx - tau) gamma_xx
Theta,t = 1/2 alpha Z_x gamma_xx,x / gamma_xx^2
		+ alpha (gamma^xx Z_x,x - Z_x gamma_xx,x / gamma_xx^2)
		+ 1/2 alpha (K_xx / gamma_xx - 2 Theta) K_xx / gamma_xx
		- 1/2 alpha K_xx^2 / gamma_xx^2
		- alpha tau
		- alpha a_x Z_x / gamma_xx
	= -1/2 alpha Z_x gamma_xx,x / gamma_xx^2
		+ alpha Z_x,x / gamma_xx
		- alpha Theta K_xx / gamma_xx
		- alpha tau
		- alpha a_x Z_x / gamma_xx
Z_x,t = alpha Theta,x
		- 2 alpha Z_x K_xx / gamma_xx
		- alpha S_x
		- alpha Theta a_x

summary:

alpha,t = -f alpha^2 (K_xx / gamma_xx - m Theta)
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


as a matrix:

[ alpha		]	[	0,	0,	0,		0,	0,					0,			0					] [ alpha		]		[ -f alpha^2 (K_xx / gamma_xx - m Theta ]
[ gamma_xx	]	[	0,	0,	0,		0,	0,					0,			0					] [ gamma_xx	] 		[ -2 alpha K_xx ]
[ a_x		]	[	0,	0,	0,		0,	f alpha / gamma_xx,	-f alpha m,	0					] [ a_x			]		[ -alpha a_x (f' alpha + f) (K_xx / gamma_xx - m Theta) + 2 f alpha K_xx d_xxx / gamma_xx^2	]
[ d_xxx		] =	[	0,	0,	0,		0,	alpha,				0,			0					] [ d_xxx		]	=	[ -alpha a_x K_xx	]
[ K_xx		]	[	0,	0,	alpha,	0,	0,					0,			-2 alpha			] [ K_xx		]		[ alpha (- a_x^2 + d_xxx (a_x - 2 Z_x) / gamma_xx - K_xx (K_xx / gamma_xx + 2 Theta) - 1/2 (S_xx + tau gamma_xx) ]
[ Theta		]	[	0,	0,	0,		0,	0,					0,			-alpha / gamma_xx	] [ Theta		]		[ -alpha ((Z_x (a_x + d_xxx / gamma_xx) + Theta K_xx) / gamma_xx + tau)	]
[ Z_x		]	[	0,	0,	0,		0,	0,					-alpha,		0					] [ Z_x			],x 	[ -alpha (2 Z_x K_xx / gamma_xx + Theta a_x + S_x) ]

V_x = -Z_x in 1D

eigenvalue = 0:
eigenfields: alpha, gamma_ij, a_p, d_pij, a_k - f d_k + f m V_k
1D: alpha, gamma_xx, a_x - f d_xxx / gamma_xx - f m Z_x

eigenvalue = +- alpha
eigenfields: (K_ij - n_i n_j K) +- (lambda^n_ij - n_i n_j tr lambda^n)
		Theta +- V^n
1D: 	K_xx - K_xx / gamma_xx +- lambda^x_xx -+ tr lambda^x
		Theta -+ Z_x / gamma_xx

eigenvalue = +- alpha sqrt(f)
eigenfields: sqrt(f) (tr K - mu Theta) +- (a^n + (2 - mu) V^n)
1D: sqrt(f) (K_xx / gamma_xx - mu Theta) +- (a_x / gamma_xx - (2 - mu) Z_x / gamma_xx)

hmm, there's 9 here, which means 2 of them don't apply to the 1D case ...

let's see what the eigenvalues of the linear system give:

eigenvalues:
matrix([-(alpha*sqrt(f))/sqrt(g_xx),0,0,0,0],[0,-alpha*sqrt(g_xx),0,0,0],[0,0,0,0,0],[0,0,0,alpha*sqrt(g_xx),0],[0,0,0,0,(alpha*sqrt(f))/sqrt(g_xx)])

right eigenvectors:
matrix([f,f*g_xx^2*m-2*f,0,f*g_xx^2*m-2*f,f],[g_xx,f*g_xx*m-2*g_xx,1,f*g_xx*m-2*g_xx,g_xx],[-sqrt(f)*sqrt(g_xx),sqrt(g_xx)*(2*g_xx-f*g_xx*m),0,sqrt(g_xx)*(f*g_xx*m-2*g_xx),sqrt(f)*sqrt(g_xx)],[0,sqrt(g_xx)*(g_xx^2-f),0,sqrt(g_xx)*(f-g_xx^2),0],[0,g_xx^2-f,0,g_xx^2-f,0])

left eigenvectors:
matrix([1/(2*f),0,-1/(2*sqrt(f)*sqrt(g_xx)),-(sqrt(g_xx)*(f*m-2))/(sqrt(f)*(2*g_xx^2-2*f)),-(g_xx^2*m-2)/(2*g_xx^2-2*f)],[0,0,0,1/(sqrt(g_xx)*(2*g_xx^2-2*f)),1/(2*g_xx^2-2*f)],[-g_xx/f,1,0,0,g_xx*m],[0,0,0,-1/(sqrt(g_xx)*(2*g_xx^2-2*f)),1/(2*g_xx^2-2*f)],[1/(2*f),0,1/(2*sqrt(f)*sqrt(g_xx)),(sqrt(g_xx)*(f*m-2))/(sqrt(f)*(2*g_xx^2-2*f)),-(g_xx^2*m-2)/(2*g_xx^2-2*f)])

eigenfields (left eigenvectors times state):
matrix([-(sqrt(g_xx)*(Z_x*f*g_xx^2*m-a_x*g_xx^2+(a_x-2*Z_x)*f)+sqrt(f)*(Theta*f*g_xx*m+K_xx*g_xx^2-2*Theta*g_xx-K_xx*f))/sqrt(g_xx)],[(Z_x*f*sqrt(g_xx)+Theta*f)/sqrt(g_xx)],[Z_x*f*g_xx*m-a_x*g_xx+d_xxx*f],[(Z_x*f*sqrt(g_xx)-Theta*f)/sqrt(g_xx)],[-(sqrt(g_xx)*(Z_x*f*g_xx^2*m-a_x*g_xx^2+(a_x-2*Z_x)*f)+sqrt(f)*(-Theta*f*g_xx*m-K_xx*g_xx^2+2*Theta*g_xx+K_xx*f))/sqrt(g_xx)])

hmm, this isn't getting the light-cone (+- alpha) eigenvalues at all, and is instead getting +- alpha sqrt(f / gamma_xx) and the expected +- alpha sqrt(gamma_xx)
hmm, Alcubierre has no eigenvalues for Z4, but for Bona-Masso ADM it has gauge as alpha sqrt(f / gamma_xx) and light as alpha / sqrt(gamma_xx)
	... while the Z4 paper says alpha sqrt(f) is gauge waves and alpha is light cone 

--]]

local class = require 'ext.class'
local Equation = require 'equation'

local Z41D = class(Equation)
Z41D.name = 'Z4-1D'

Z41D.numStates = 7
Z41D.numWaves = 5	-- no waves for alpha and gamma_xx, which are purely source-driven

local m = 2	-- "for f=1 one must have m=2"
local tau = 0
local S_x = 0
local S_xx = 0

-- initial conditions
function Z41D:init(args, ...)

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
	local d_xxx = q:_(4)
	local D_g = 2 * d_xxx / gamma_xx
	local K_xx = q:_(5)
	local KTilde_xx = K_xx * math.sqrt:o(gamma_xx)
	local K = K_xx / gamma_xx
	local volume = alpha * math.sqrt:o(gamma_xx)
	local Theta = q:_(6)
	local Z_x = q:_(7)
	Z41D:buildGraphInfos{
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

function Z41D:initCell(sim,i)
	local x = sim.xs[i]
	local alpha = self.calc.alpha(x)
	local gamma_xx = self.calc.gamma_xx(x)
	local a_x = self.calc.dx_alpha(x) / self.calc.alpha(x)
	local d_xxx = 1/2 * self.calc.dx_gamma_xx(x)
	local K_xx = self.calc.K_xx(x) 
	-- what is Theta and Z_i initialized to?
	local Theta = 0	
	local Z_x = 0
	return {alpha, gamma_xx, a_x, d_xxx, K_xx, Theta, Z_x}
end

function Z41D:calcEigenvalues(alpha, gamma_xx, f)
	local f = self.calc.f(alpha)
	local sqrt_gamma_xx = math.sqrt(gamma_xx)
	local lambda1 = alpha * math.sqrt(f) / sqrt_gamma_xx
	local lambda2 = alpha * sqrt_gamma_xx
	return -lambda1, -lambda2, 0, lambda2, lambda1
end

-- returns averaging of variables used for interface eigenbasis
function Z41D:calcRoeValues(qL, qR)
	local alpha = (qL[1] + qR[1]) / 2
	local gamma_xx = (qL[2] + qR[2]) / 2
	local f = self.calc.f(alpha)
	return alpha, gamma_xx, f	
end

function Z41D:calcEigenBasis(lambda, evr, evl, dF_dU, alpha, gamma_xx, f)
	fill(lambda, self:calcEigenvalues(alpha, gamma_xx, f))
	fill(evl, f, gamma_xx)
	fill(evr, f, gamma_xx)
	if dF_dU then fill(dF_dU, alpha, gamma_xx, f) end
end

function Z41D:fluxMatrixTransform(solver, A, v)
	local alpha, gamma_xx, f = table.unpack(A)
	local _, _, v1, v2, v3, v4, v5 = table.unpack(v)	-- skip alpha, gamma_xx, and transform the rest: a_x, d_xxx, K_xx, Theta, Z_x
	return {
		0,
		0,
		(alpha*f*v3)/gamma_xx-alpha*f*m*v4,
		alpha*v3,
		alpha*v1-2*alpha*v5,
		-alpha*gamma_xx*v5,
		-alpha*v4,
	}
end

function Z41D:eigenLeftTransform(solver, evL, v)
	local f, gamma_xx = table.unpack(evL)
	local _, _, v1, v2, v3, v4, v5 = table.unpack(v)	-- skip alpha, gamma_xx, and transform the rest: a_x, d_xxx, K_xx, Theta, Z_x
	return {
		-((gamma_xx^2*m-2)*v5)/(2*gamma_xx^2-2*f)-(math.sqrt(gamma_xx)*(f*m-2)*v4)/(math.sqrt(f)*(2*gamma_xx^2-2*f))-v3/(2*math.sqrt(f)*math.sqrt(gamma_xx))+v1/(2*f),
		v5/(2*gamma_xx^2-2*f)+v4/(math.sqrt(gamma_xx)*(2*gamma_xx^2-2*f)),
		gamma_xx*m*v5+v2-(gamma_xx*v1)/f,
		v5/(2*gamma_xx^2-2*f)-v4/(math.sqrt(gamma_xx)*(2*gamma_xx^2-2*f)),
		-((gamma_xx^2*m-2)*v5)/(2*gamma_xx^2-2*f)+(math.sqrt(gamma_xx)*(f*m-2)*v4)/(math.sqrt(f)*(2*gamma_xx^2-2*f))+v3/(2*math.sqrt(f)*math.sqrt(gamma_xx))+v1/(2*f),
	}
end

function Z41D:eigenRightTransform(solver, evR, v)
	local f, gamma_xx = table.unpack(evR)
	local v1, v2, v3, v4, v5 = table.unpack(v)
	return {
		0,
		0,
		f*v5+(f*gamma_xx^2*m-2*f)*v4+(f*gamma_xx^2*m-2*f)*v2+f*v1,
		gamma_xx*v5+(f*gamma_xx*m-2*gamma_xx)*v4+v3+(f*gamma_xx*m-2*gamma_xx)*v2+gamma_xx*v1,
		math.sqrt(f)*math.sqrt(gamma_xx)*v5+math.sqrt(gamma_xx)*(f*gamma_xx*m-2*gamma_xx)*v4+math.sqrt(gamma_xx)*(2*gamma_xx-f*gamma_xx*m)*v2-math.sqrt(f)*math.sqrt(gamma_xx)*v1,
		math.sqrt(gamma_xx)*(f-gamma_xx^2)*v4+math.sqrt(gamma_xx)*(gamma_xx^2-f)*v2,
		(gamma_xx^2-f)*v4+(gamma_xx^2-f)*v2,
	}
end

function Z41D:calcCellMinMaxEigenvalues(sim, i)
	local alpha, gamma_xx = table.unpack(sim.qs[i])
	local f = self.calc.f(alpha)
	return firstAndLast(self:calcEigenvalues(alpha, gamma_xx, f))
end

function Z41D:sourceTerm(sim, qs)
	local source = sim:newState()
	for i=1,sim.gridsize do
		local alpha, gamma_xx, a_x, d_xxx, K_xx, Theta, Z_x = table.unpack(qs[i])
		local f = self.calc.f(alpha)
		local dalpha_f = self.calc.dalpha_f(alpha)
		
		source[i][1] = -f * alpha * alpha * (K_xx / gamma_xx - m * Theta)
		source[i][2] = -2 * alpha * K_xx
		source[i][3] = -alpha * a_x * (dalpha_f * alpha + f) * (K_xx / gamma_xx - m * Theta) + 2 * f * alpha * K_xx * d_xxx / (gamma_xx * gamma_xx)
		source[i][4] = -alpha * a_x * K_xx
		source[i][5] = alpha * (-a_x * a_x + d_xxx * (a_x - 2 * Z_x) / gamma_xx - K_xx * (K_xx / gamma_xx + 2 * Theta) - 1/2 * (S_xx + tau * gamma_xx))
		source[i][6] = -alpha * ((Z_x * (a_x + d_xxx / gamma_xx) + Theta * K_xx) / gamma_xx + tau)
		source[i][7] = -alpha * (2 * Z_x * K_xx / gamma_xx + Theta * a_x + S_x)
	end
	return source
end

return Z41D
