--[[
based on the 5-var system, but with the hyperbolic formalism only based on A D K while alpha and g are updated separately

A = (ln alpha),x = alpha,x / alpha
alpha,x = alpha A

D = 1/2 g,x
g,x = 2 D

prim evolution:
d_t alpha = -alpha^2 f K / g
d_t g = -2 alpha K
state evolution:
d_t A + alpha K / g f,x + alpha f / g K,x = alpha K f / g (2 D / g - A)
d_t D + alpha K,x = -alpha K A
d_t K + alpha A,x = alpha ((A D - K^2) / g - A^2)

linear system:
       [0     0  alpha f / g]     [A]   [ alpha K f / g (2 D / g - A) ]
d_t  + [0     0    alpha    ] d_x [D] = [          -alpha K A         ]
       [alpha 0      0      ]     [K]   [alpha ((A D - K^2) / g - A^2)]
--]]

require 'ext'
local Simulation = require 'simulation'

local ADM1D3VarSim = class(Simulation)

ADM1D3VarSim.numStates = 3

-- initial conditions
function ADM1D3VarSim:init(args, ...)
	Simulation.init(self, args, ...)

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

	local f = symmath.clone(assert(args.alpha)):simplify()
	self.calc_f = f:compile{assert(args.f_var)}

	self.graphInfos = {
		{viewport={0/3, 0/3, 1/3, 1/3}, getter=index:bind(self.qs):index'alpha', name='alpha', color={1,0,1}},
		{viewport={0/3, 1/3, 1/3, 1/3}, getter=index:bind(self.qs):index(1), name='A', color={0,1,0}},
		{viewport={1/3, 0/3, 1/3, 1/3}, getter=index:bind(self.qs):index'g', name='g', color={.5,.5,1}},
		{viewport={1/3, 1/3, 1/3, 1/3}, getter=index:bind(self.qs):index(2), name='D', color={1,1,0}},
		{viewport={2/3, 0/3, 1/3, 1/3}, getter=index:bind(self.qs):index(3), name='K', color={0,1,1}},
		{viewport={0/3, 2/3, 1/3, 1/3}, getter=log:compose(index:bind(self.eigenbasisErrors)), name='error', color={1,0,0}, range={-30, 30}},
		{viewport={1/3, 2/3, 1/3, 1/3}, getter=log:compose(index:bind(self.fluxMatrixErrors)), name='error', color={1,0,0}, range={-30, 30}},
	}
end

function ADM1D3VarSim:initCell(i)
	local x = self.xs[i]
	local alpha = self.calc_alpha(x)
	local g = self.calc_g(x)
	local A = self.calc_dx_alpha(x) / self.calc_alpha(x)
	local D = 1/2 * self.calc_dx_g(x)
	local K = -self.calc_d2x_h(x) / sqrt(self.calc_g(x))
	return {A, D, K, alpha=alpha, g=g}
end

function ADM1D3VarSim:calcInterfaceEigenBasis(i)
	local avgQ = {}
	for j=1,self.numStates do 
		avgQ[j] = (self.qs[i-1][j] + self.qs[i][j]) / 2
	end
	avgQ.alpha = (self.qs[i-1].alpha + self.qs[i].alpha) / 2
	avgQ.g = (self.qs[i-1].g + self.qs[i].g) / 2

	local A, D, K = unpack(avgQ)
	local x = self.ixs[i]
	local alpha = avgQ.alpha
	local g = avgQ.g
	local f = self.calc_f(alpha)
	local lambda = alpha * sqrt(f / g)		
	self.eigenvalues[i] = {-lambda, 0, 0, 0, lambda}
	-- row-major, math-indexed
	self.fluxMatrix[i] = {
		{0,0, alpha*f/g},
		{0,0,alpha},
		{alpha,0,0},
	}
	self.eigenvectors[i] = {
		{f,			 0,	f			},
		{g,			 1,	g			},
		{-sqrt(f*g), 0,	sqrt(f*g)	},
	}
	self.eigenvectorsInverse[i] = {
		{1/(2*f), 	0, -1/(2*sqrt(f*g))	},
		{-g/f, 		1, 0				},
		{1/(2*f), 	0, 1/(2*sqrt(f*g))	},
	}
end

function ADM1D3VarSim:zeroDeriv()
	ADM1D3VarSim.super.zeroDeriv(self)
	-- zero deriv
	for i=1,self.gridsize do
		self.dq_dts[i].alpha = 0
		self.dq_dts[i].g = 0
	end
end

function ADM1D3VarSim:addSourceToDerivCell(i)
	local A, D, K = unpack(self.qs[i])
	local alpha = self.qs[i].alpha
	local g = self.qs[i].g
	local f = self.calc_f(alpha)
	self.dq_dts[i].alpha = self.dq_dts[i].alpha - alpha * alpha * f * K / g
	self.dq_dts[i].g = self.dq_dts[i].g - 2 * alpha * K
	self.dq_dts[i][1] = self.dq_dts[i][1] + alpha * K * f / g * (2 * D / g - A)
	self.dq_dts[i][2] = self.dq_dts[i][2] - alpha * K * A
	self.dq_dts[i][3] = self.dq_dts[i][3] + alpha * ((A * D - K * K) / g - A * A)
end

function ADM1D3VarSim:integrateDeriv(dt)
	ADM1D3VarSim.super.integrateDeriv(self, dt)
	for i=1,self.gridsize do
		self.qs[i].alpha = self.qs[i].alpha + dt * self.dq_dts[i].alpha
		self.qs[i].g = self.qs[i].g + dt * self.dq_dts[i].g
	end
	self.t = self.t + dt
end

return ADM1D3VarSim

