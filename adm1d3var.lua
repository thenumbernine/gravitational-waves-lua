--[[
based on the book "Introduction to 3+1 Numerical Relativity" and on the paper "Introduction to Numerical Relativity", both by Alcubierre

A = (ln alpha),x = alpha,x / alpha
D = (ln g),x = g,x / g
KTilde = sqrt(g) K

alpha,x = A alpha
g,x = g D
K = KTilde / sqrt(g)

A,t + (alpha f K),x = 0
D,t + (2 alpha K),x = 0
K,t + (alpha A / g),x = alpha (K^2 - A D / (2 g))

...rewritten for KTilde...
A,t + (alpha f KTilde / sqrt(g)),x = 0
D,t + (2 alpha KTilde / sqrt(g)),x = 0
KTilde,t + (alpha A / sqrt(g)),x = 0

A,t + (f KTilde/sqrt(g) + alpha KTilde/sqrt(g) f,alpha) alpha,x + alpha f / sqrt(g) KTilde,x - 1/2 alpha f KTilde / g^(3/2) g,x = 0
D,t + 2 KTilde/sqrt(g) alpha,x + 2 alpha/sqrt(g) KTilde,x - alpha KTilde / g^(3/2) g,x = 0
KTilde,t + A / sqrt(g) alpha,x + alpha / sqrt(g) A,x - 1/2 alpha A / g^(3/2) g,x = 0

A,t + alpha f / sqrt(g) KTilde,x = alpha KTilde / sqrt(g) (f (1/2 D - A) - A alpha f,alpha)
D,t + 2 alpha / sqrt(g) KTilde,x = 2 alpha KTilde / sqrt(g) (1/2 D - A)
KTilde,t + alpha / sqrt(g) A,x = A alpha / sqrt(g) (1/2 D - A)

[  A   ]     [0,               0, alpha f / sqrt(g)] [  A   ]     [alpha KTilde / sqrt(g) (f (1/2 D - A) - A alpha f,alpha)]
[  D   ],t + [0,               0, 2 alpha / sqrt(g)] [  D   ],x = [2 alpha KTilde / sqrt(g) (1/2 D - A)                    ]
[KTilde]     [alpha / sqrt(g), 0, 0                ] [KTilde]     [A alpha / sqrt(g) (1/2 D - A)                           ]

...and this is the matrix he has.
now we look at eigenvectors ...

/* [wxMaxima: input   start ] */
loag("eigen");
/* [wxMaxima: input   end   ] */
/* [wxMaxima: input   start ] */
assume(alpha>0,f>0);
/* [wxMaxima: input   end   ] */
/* [wxMaxima: input   start ] */
F : matrix(
[0, 0, alpha *f/sqrt(g)],
[0, 0, 2*alpha/sqrt(g)],
[alpha/sqrt(g), 0, 0]
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

	local x = assert(args.x)

	local h = symmath.clone(assert(args.h)):simplify()
	self.calc_h = h:compile{x}
	
	local dx_h = h:diff(x):simplify()
	self.calc_dx_h = dx_h:compile{x}
	
	local d2x_h = dx_h:diff(x):simplify()
	self.calc_d2x_h = d2x_h:compile{x}

	local g = symmath.clone(assert(args.g)):simplify()
	self.calc_g = g:compile{x}

	local dx_g = g:diff(x):simplify()
	self.calc_dx_g = dx_g:compile{x}

	local alpha = symmath.clone(assert(args.alpha)):simplify()
	self.calc_alpha = alpha:compile{x}

	local dx_alpha = alpha:diff(x):simplify()
	self.calc_dx_alpha = dx_alpha:compile{x}

	local f_param = assert(args.f_param)

	local f = symmath.clone(assert(args.alpha)):simplify()
	self.calc_f = f:compile{f_param}

	local dalpha_f = f:diff(f_param):simplify()
	self.calc_dalpha_f = dalpha_f:compile{f_param}

	local get_state = index:bind(self.qs)
	local get_alpha = get_state:index'alpha'
	local get_g = get_state:index'g'
	local get_A = get_state:index(1)
	local get_D = get_state:index(2)
	local get_KTilde = get_state:index(3)
	self.graphInfos = {
		{viewport={0/3, 0/3, 1/3, 1/3}, getter=get_alpha, name='alpha', color={1,0,1}},
		{viewport={0/3, 1/3, 1/3, 1/3}, getter=get_A, name='A', color={0,1,0}},
		{viewport={1/3, 0/3, 1/3, 1/3}, getter=get_g, name='g', color={.5,.5,1}},
		{viewport={1/3, 1/3, 1/3, 1/3}, getter=get_D, name='D', color={1,1,0}},
		{viewport={2/3, 0/3, 1/3, 1/3}, getter=get_KTilde, name='KTilde', color={0,1,1}},
		{viewport={2/3, 1/3, 1/3, 1/3}, getter=get_alpha * sqrt:compose(get_g), name='volume', color={0,1,1}},
		{viewport={0/3, 2/3, 1/3, 1/3}, getter=log:compose(index:bind(self.eigenbasisErrors)), name='log eigenbasis error', color={1,0,0}, range={-30, 30}},
		{viewport={1/3, 2/3, 1/3, 1/3}, getter=log:compose(index:bind(self.fluxMatrixErrors)), name='log reconstuction error', color={1,0,0}, range={-30, 30}},
	}
end

function ADM1D3VarSim:initCell(i)
	local x = self.xs[i]
	local alpha = self.calc_alpha(x)
	local g = self.calc_g(x)
	local A = self.calc_dx_alpha(x) / self.calc_alpha(x)
	local D = 1/2 * self.calc_dx_g(x)
	local K = -self.calc_d2x_h(x) / sqrt(self.calc_g(x))
	local KTilde = K / sqrt(g)
	return {A, D, KTilde, alpha=alpha, g=g}
end

function ADM1D3VarSim:calcInterfaceEigenBasis(i)
	local alpha = (self.qs[i-1].alpha + self.qs[i].alpha) / 2
	local g = (self.qs[i-1].g + self.qs[i].g) / 2
	local f = self.calc_f(alpha)
	local lambda = alpha * sqrt(f / g)		
	self.eigenvalues[i] = {-lambda, 0, lambda}

	local function buildField(call)
		return function(i, ...)
			local v1, v2, v3 = ...
			
			local avgQ = {}
			for j=1,self.numStates do 
				avgQ[j] = (self.qs[i-1][j] + self.qs[i][j]) / 2
			end
			avgQ.alpha = (self.qs[i-1].alpha + self.qs[i].alpha) / 2
			avgQ.g = (self.qs[i-1].g + self.qs[i].g) / 2
			
			local A, D, KTilde = unpack(avgQ)
			local x = self.ixs[i]
			local alpha = avgQ.alpha
			local g = avgQ.g
			local f = self.calc_f(alpha)

			return call(alpha, f, g, A, D, KTilde, v1, v2, v3)
		end
	end

	self:buildFields{
		fluxTransform = buildField(function(alpha, f, g, A, D, KTilde, v1, v2, v3)
			return
				v3 * alpha * f / sqrt(g),
				v3 * 2 * alpha / sqrt(g),
				v1 * alpha / sqrt(g)
		end),
		eigenfields = buildField(function(alpha, f, g, A, D, KTilde, v1, v2, v3)
			return
				v1 / (2 * f) - v3 / (2 * sqrt(f)),
				-2*v1/f + v2,
				v1 / (2 * f) + v3 / (2 * sqrt(f))
		end),
		eigenfieldsInverse = buildField(function(alpha, f, g, A, D, KTilde, v1, v2, v3)
			return
				(v1 + v3) * f,
				2 * v1 + v2 + 2 * v3,
				sqrt(f) * (-v1 + v3)
		end),
	} 
end

function ADM1D3VarSim:zeroDeriv(dq_dts)
	ADM1D3VarSim.super.zeroDeriv(self, dq_dts)
	-- zero deriv
	for i=1,self.gridsize do
		dq_dts[i].alpha = 0
		dq_dts[i].g = 0
	end
end

function ADM1D3VarSim:addSourceToDerivCell(dq_dts, i)
	local A, D, KTilde = unpack(self.qs[i])
	local alpha = self.qs[i].alpha
	local g = self.qs[i].g
	local f = self.calc_f(alpha)
	local dalpha_f = self.calc_dalpha_f(alpha)
	
	dq_dts[i].alpha = dq_dts[i].alpha - alpha * alpha * f * KTilde / (g * sqrt(g))
	dq_dts[i].g = dq_dts[i].g - 2 * alpha * KTilde / sqrt(g)
	
	local tmp1 = alpha / sqrt(g)
	local tmp2 = .5 * D - A
	dq_dts[i][1] = dq_dts[i][1] + KTilde * tmp1 * (f * tmp2 - A * alpha * dalpha_f)
	dq_dts[i][2] = dq_dts[i][2] + 2 * KTilde * tmp1 * tmp2
	dq_dts[i][3] = dq_dts[i][3] + A * tmp1 * tmp2
end

function ADM1D3VarSim:integrateDeriv(dq_dts, dt)
	ADM1D3VarSim.super.integrateDeriv(self, dq_dts, dt)
	for i=1,self.gridsize do
		self.qs[i].alpha = self.qs[i].alpha + dt * dq_dts[i].alpha
		self.qs[i].g = self.qs[i].g + dt * dq_dts[i].g
	end
	self.t = self.t + dt
end

return ADM1D3VarSim

