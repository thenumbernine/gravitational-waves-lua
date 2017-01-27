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
Theta,t = -1/2 alpha Z_x gamma_xx,x / gamma_xx^2
		+ alpha Z_x,x / gamma_xx
		- alpha Theta K_xx / gamma_xx
		- alpha tau
		- alpha a_x Z_x / gamma_xx
Z_x,t = alpha Theta,x
		- 2 alpha Z_x K_xx / gamma_xx
		- alpha S_x
		- alpha Theta a_x
--]]

local class = require 'ext.class'
local Equation = require 'equation'

local Z41D = class(Equation)
Z41D.name = 'Z4-1D'

Z41D.numStates = 7
Z41D.numWaves = 5	-- no waves for alpha and gamma_xx, which are purely source-driven

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


