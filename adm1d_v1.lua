--[[

Based on the book "Introduction to 3+1 Numerical Relativity" and on the paper "Introduction to Numerical Relativity", both by Alcubierre

In the "Toy 1+1" section the book uses D_alpha = (ln alpha),x
 but the "Hyperbolic Formalisms" sections use a_i = (ln alpha),i ... 
a_x = (ln alpha),x = alpha,x / alpha
so alpha,x = a_x alpha

D_g = (ln gamma_xx),x 
	= gamma_xx,x / gamma_xx
so gamma_xx,x = gamma_xx D_g

note in 1D, gamma = gamma_xx

The d_ijk variable is also used in the "Hyperbolic Formalisms" chapter:
d_xxx = 1/2 gamma_xx,x = 1/2 gamma_xx D_g

KTilde = sqrt(gamma_xx) K
K = KTilde / sqrt(gamma_xx)
K = gamma^ij K_ij = gamma^xx K_xx = K_xx / gamma_xx 
K_xx / gamma_xx = KTilde / sqrt(gamma_xx)
K_xx = sqrt(gamma_xx) KTilde
KTilde = K_xx / sqrt(gamma_xx)

ADM:

alpha,t = -alpha^2 f K

gamma_ij,t = -2 alpha K_ij

gamma,t = gamma gamma^ij gamma_ij,t = -2 alpha gamma gamma^ij K_ij = -2 alpha gamma K

K_ij,t = -D_i D_j alpha + alpha (R_ij + K K_ij - 2 K_ik K^k_j)
K_ij,t = -alpha,ij + Gamma^k_ij alpha,k + alpha (R_ij + K K_ij - 2 K_ik K^k_j)

K,t = (gamma^ij K_ij),t 
	= gamma^ij,t K_ij + gamma^ij K_ij,t
	= -gamma^im gamma_mn,t gamma^nj K_ij + gamma^ij K_ij,t
	= gamma^ij K_ij,t - K^ij gamma_ij,t
	= gamma^ij (-alpha,ij + Gamma^k_ij alpha,k + alpha (R_ij + K K_ij - 2 K_ik K^k_j)) - K^ij (-2 alpha K_ij)
	= -gamma^ij alpha,ij + gamma^ij Gamma^k_ij alpha,k + alpha (R + K^2)
... needs to become ...
... so R must become 0 ...
	= -gamma^ij alpha,ij + gamma^ij gamma^kl Gamma_lij alpha,k + alpha K^2
... for 1D, Gamma_xxx = 1/2 gamma_xx,x so Gamma^x_xx = 1/2 gamma^xx gamma_xx,x
	= -gamma^xx alpha,xx + gamma^xx * gamma^xx * 1/2 gamma_xx,x alpha,x + alpha K^2
	= -alpha,xx / gamma + 1/2 alpha,x gamma,x / gamma^2 + alpha K^2
	= -alpha,xx / gamma + alpha,x gamma,x / gamma^2 + alpha K^2 - alpha,x gamma,x / (2 gamma^2)
	= -(alpha,x / gamma),x + (alpha K^2 - alpha,x gamma,x / (2 gamma^2))
K,t + (alpha a_x / gamma),x = alpha (K^2 - a_x D_g / (2 gamma))

in 1D variables:

a_x,t = (ln alpha),xt 
	= ((ln alpha),t),x 
	= (alpha,t / alpha),x 
	= (-alpha f K),x

gamma_xx,t = -2 alpha K_xx
(ln gamma_xx),t = gamma_xx,t / gamma_xx = -2 alpha K_xx / gamma_xx 
	= -2 alpha gamma^xx K_xx = -2 alpha K

D_g,t = (ln gamma_xx),xt = ((ln gamma_xx),t),x = (-2 alpha K),x

KTilde,t = (sqrt(gamma_xx) K),t 
	= sqrt(gamma_xx),t K + sqrt(gamma_xx) K,t
	= 1/2 K gamma_xx,t / sqrt(gamma_xx) + sqrt(gamma_xx) K,t
	= 1/2 K (-2 alpha gamma_xx K) / sqrt(gamma_xx) - sqrt(gamma_xx) ((alpha,x / gamma),x + (alpha K^2 - alpha,x gamma,x / (2 gamma^2)))
	= -2 alpha sqrt(gamma_xx) K^2 
		- sqrt(gamma_xx) (alpha,x / gamma),x 
		+ 1/2 alpha,x gamma,x / gamma_xx^(3/2)
...needs to become...
KTilde,t + (alpha a_x / sqrt(gamma_xx)),x = 0

1D ADM formalism:
a_x,t + (alpha f K),x = 0
D_g,t + (2 alpha K),x = 0
K_xx,t + (alpha a_x / gamma_xx),x = alpha (K_xx^2 - a_x D_g / (2 gamma_xx))

...rewritten for D_g and KTilde...
a_x,t + (alpha f KTilde / sqrt(gamma_xx)),x = 0
D_g,t + (2 alpha KTilde / sqrt(gamma_xx)),x = 0
KTilde,t + (alpha a_x / sqrt(gamma_xx)),x = 0

...linearized

a_x,t + alpha,x f KTilde / sqrt(gamma_xx) 
	+ alpha f,x KTilde / sqrt(gamma_xx) 
	+ alpha f KTilde,x / sqrt(gamma_xx)
	- 1/2 alpha f KTilde gamma_xx,x / gamma_xx^(3/2) = 0
a_x,t + alpha f KTilde,x / sqrt(gamma_xx) = (1/2 f D_g - (f + alpha f') a_x) alpha KTilde / sqrt(gamma_xx)

D_g,t + 2 alpha,x KTilde / sqrt(gamma_xx)
	+ 2 alpha KTilde,x / sqrt(gamma_xx)
	- alpha KTilde gamma_xx,x / gamma_xx^(3/2) = 0
D_g,t + 2 alpha KTilde,x / sqrt(gamma_xx) = (D_g - 2 a_x) alpha KTilde / sqrt(gamma_xx)

KTilde,t + alpha,x a_x / sqrt(gamma_xx)
	+ alpha a_x,x / sqrt(gamma_xx)
	+ alpha a_x (-1/2 gamma_xx,x / gamma_xx^(3/2))
	= 0
KTilde,t + alpha a_x,x / sqrt(gamma_xx) = (1/2 D_g - a_x) alpha a_x / sqrt(gamma_xx)

... as a matrix

[  a_x ]   [        0               0, alpha f / sqrt(gamma_xx)] [  a_x ]   [ ((1/2 D_g - a_x) f - alpha a_x f' ) KTilde alpha / sqrt(gamma_xx)
[  D_g ] + [        0               0, 2 alpha / sqrt(gamma_xx)] [  D_g ] = [ (1/2 D_g - a_x) 2 KTilde alpha / sqrt(gamma_xx)
[KTilde]   [alpha / sqrt(gamma_xx), 0,           0             ] [KTilde],x [ (1/2 D_g - a_x) a_x alpha / sqrt(gamma_xx)

	here's some Maxima code to verify:

depends([alpha, KTilde, gamma_xx], x);
depends(f, alpha);
F : transpose(matrix([
	alpha * f * KTilde / sqrt(gamma_xx),
	2 * alpha * KTilde / sqrt(gamma_xx),
	alpha * a_x / sqrt(gamma_xx)
]));
diff(F, x);
subst([diff(gamma_xx,x) = gamma_xx * D_g, diff(alpha, x) = alpha * a_x], %);

what exactly is 1/2 D_g - a_x ?

1/2 gamma_xx,x / gamma_xx - alpha,x / alpha
(1/2 ln(gamma_xx) - ln(alpha)),x
ln(sqrt(gamma_xx) / alpha),x

x-derivative of the logarithm of the ratio of the spatial volume to lapse

now we look at eigenvectors ...

	Maxima code:

loag("eigen");
assume(alpha>0,f>0);
F : matrix(
	[0, 0, alpha*f/sqrt(gamma_xx)],
	[0, 0, 2*alpha/sqrt(gamma_xx)],
	[alpha/sqrt(gamma_xx), 0, 0]
);
results:eigenvectors(F);
eigenvalues : diag_matrix(results[1][1][1], results[1][1][3], results[1][1][2]);
eigenvectors : transpose(matrix(
	results[2][1][1] * f,
	results[2][3][1],
	results[2][2][1] * f
));
determinant(eigenvectors);
invert(eigenvectors);

--]]

local class = require 'ext.class'
local Equation = require 'equation'

local ADM1Dv1 = class(Equation)
ADM1Dv1.name = 'ADM 1D v.1'

ADM1Dv1.numStates = 5
ADM1Dv1.numWaves = 3

-- initial conditions
function ADM1Dv1:init(args, ...)

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
	local _2dx = function(self,i) return self.xs[i+1] - self.xs[i-1] end
	local alpha = q:_(1)
	local gamma_xx = q:_(2)
	local a_x = q:_(3)
	local D_g = q:_(4)
	local d_xxx = D_g * gamma_xx / 2
	local KTilde = q:_(5)
	local K_xx = KTilde * math.sqrt:o(gamma_xx)
	local K = K_xx / gamma_xx
	local volume = alpha * math.sqrt:o(gamma_xx)
	local f = function(self, i) return self.equation.calc.f(self.qs[i][1]) end
	local dalpha_f = function(self, i) return self.equation.calc.dalpha_f(self.qs[i][1]) end
	ADM1Dv1:buildGraphInfos{
		{alpha = alpha},
		{a_x = a_x},
		{gamma_xx = gamma_xx},
		{d_xxx = d_xxx},
		{D_g = D_g},
		{K_xx = K_xx},
		{KTilde = KTilde},
		{K = K},
		{f = f},
		{dalpha_f = dalpha_f},
		{volume = volume},
	
		{['log alpha vs a_x'] = function(self, i)
			-- a_x = (ln alpha),x = alpha,x / alpha
			-- alpha a_x = alpha,x
			local dx_alpha = (alpha(self,i+1) - alpha(self,i-1)) / _2dx(self,i)
			return math.log(math.abs(dx_alpha - alpha(self,i) * a_x(self,i)), 10)
		end},
	
		{['log gamma_xx vs D_g'] = function(self, i)
			--D_g = (ln gamma_xx),x  = gamma_xx,x / gamma_xx
			-- gamma_xx,x = gamma_xx D_g
			local dx_gamma_xx = (gamma_xx(self,i+1) - gamma_xx(self,i-1)) / _2dx(self,i)
			return math.log(math.abs(dx_gamma_xx - gamma_xx(self,i) * D_g(self,i)), 10)
		end},
	}
end

function ADM1Dv1:initCell(sim,i)
	local x = sim.xs[i]
	local alpha = self.calc.alpha(x)
	local gamma_xx = self.calc.gamma_xx(x)
	local a_x = self.calc.dx_alpha(x) / self.calc.alpha(x)
	local D_g = self.calc.dx_gamma_xx(x) / self.calc.gamma_xx(x)
	local K_xx = self.calc.K_xx(x) 
	local KTilde = K_xx / math.sqrt(gamma_xx)
	return {alpha, gamma_xx, a_x, D_g, KTilde}
end

function ADM1Dv1:calcEigenvalues(alpha, gamma_xx, f)
	local lambda = alpha * math.sqrt(f / gamma_xx)
	return -lambda, 0, lambda
end

-- arithmetic
function ADM1Dv1:calcRoeValues(qL, qR)
	local alpha = (qL[1] + qR[1]) / 2
	local gamma_xx = (qL[2] + qR[2]) / 2
	local f = self.calc.f(alpha)
	return alpha, gamma_xx, f
end

function ADM1Dv1:calcEigenBasis(lambda, evr, evl, dF_dU, alpha, gamma_xx, f)
	-- store eigenvalues
	fill(lambda, self:calcEigenvalues(alpha, gamma_xx, f))
	-- store information needed to build left and right eigenvectors
	-- this is why I need an 'eigenbasis' object - to hold this information once
	fill(evl, f)
	fill(evr, f)
	if dF_dU then fill(dF_dU, alpha, gamma_xx, f) end
end

-- how can the flux depend on gamma_xx and alpha but not the eigenvectors?
function ADM1Dv1:fluxMatrixTransform(solver, m, v)
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

function ADM1Dv1:eigenLeftTransform(solver, m, v)
	local f = table.unpack(m)
	local _, _, v1, v2, v3 = table.unpack(v)
	return {
		v1 / (2 * f) - v3 / (2 * math.sqrt(f)),
		-2*v1/f + v2,
		v1 / (2 * f) + v3 / (2 * math.sqrt(f))
	}
end

function ADM1Dv1:eigenRightTransform(solver, m, v)
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

function ADM1Dv1:calcCellMinMaxEigenvalues(sim, i)
	local alpha, gamma_xx = table.unpack(sim.qs[i])
	local f = self.calc.f(alpha)
	return firstAndLast(self:calcEigenvalues(alpha, gamma_xx, f))
end

function ADM1Dv1:sourceTerm(sim, qs, dt)
	local source = sim:newState()
	for i=2,sim.gridsize-1 do
		local alpha, gamma_xx, a_x, D_g, KTilde = table.unpack(qs[i])
		local f = self.calc.f(alpha)
		local dalpha_f = self.calc.dalpha_f(alpha)
		local K = KTilde / math.sqrt(gamma_xx)

		source[i][1] = -alpha * alpha * f * K
		source[i][2] = -2 * alpha * gamma_xx * K
		source[i][3] = ((1/2 * D_g - a_x) * f - alpha * dalpha_f * a_x) * alpha * K
		source[i][4] = (1/2 * D_g - a_x) * 2 * alpha * K
		source[i][5] = (1/2 * D_g - a_x) * a_x * alpha / math.sqrt(gamma_xx)

		-- [[ using dissipation for damping.  but why not just directly constrain the variables?
		-- modification from Bona et al "Elements of Numerical Relativity..." 2009 eqns 4.9 & 4.10
		local eta = 1/dt	-- damping term / constraint enforcing of 1st order terms
		local _2dx = sim.xs[i+1] - sim.xs[i-1]
		-- a_x = alpha,x / alpha <=> a_x += eta (alpha,x / alpha - a_x)
		local dx_alpha = (sim.qs[i+1][1] - sim.qs[i-1][1]) / _2dx
		source[i][3] = source[i][3] + eta * (dx_alpha / alpha - a_x)
		
		--D_g = gamma_xx,x / gamma_xx <=> D_g += eta (gamma_xx,x / gamma_xx - D_g)
		local dx_gamma_xx = (sim.qs[i+1][2] - sim.qs[i-1][2]) / _2dx
		source[i][4] = source[i][4] + eta * (dx_gamma_xx / gamma_xx - D_g)
		--]]
	end
	return source
end

--[[ directly constraining the first-order variable constraints
--[=[ nullspace of a rectangular matrix of finite difference of alpha,x; then -alpha a_x
[h/2 0 -h/2   ... | -a_x_1 0 ... ] [ a_x_1 ]   [ 0 ]
[0 h/2 0 -h/2 ... | 0 -a_x_2 ... ] [ a_x_2 ] = [ 0 ]
 ...                               [  ...  ]
                                   [alpha_1]
					   			   [alpha_2]
                                   [  ...  ]
--]=]
function ADM1Dv1:postIterate(sim, qs)
	for i=2,sim.gridsize-1 do
	end
end
--]]

return ADM1Dv1
