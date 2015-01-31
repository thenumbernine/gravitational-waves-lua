require 'ext'
local Simulation = require 'adm1d.simulation'

local ADM1D3VarSim = class(Simulation)

ADM1D3VarSim.numStates = 3

-- initial conditions
function ADM1D3VarSim:init(...)
	Simulation.init(self, ...)

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

do
	local sigma = 10
	local xc = 150
	local H = 5
	-- aux var for init
	local function calc_h(x) return H * exp(-(x - xc)^2 / sigma^2) end
	local function d_calc_h(x) return -2 * (x - xc) / sigma^2 * calc_h(x) end
	local function d2_calc_h(x) return (-2 / sigma^2 + 4 * (x - xc)^2 / sigma^4) * calc_h(x) end
	local function calc_g(x) return 1 - d_calc_h(x)^2 end
	local function d_calc_g(x) return -2 * d_calc_h(x) * d2_calc_h(x) end
	local function calc_alpha(x) return 1 end
	local function d_calc_alpha(x) return 0 end

	function ADM1D3VarSim:initCell(i)
		local x = self.xs[i]
		-- state variables:
		local alpha = calc_alpha(x)
		local g = calc_g(x)
		local A = d_calc_alpha(x) / alpha
		local D = .5 * d_calc_g(x)
		local K = -d2_calc_h(x) / sqrt(calc_g(x))
		return {A, D, K, alpha=alpha, g=g}
	end
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
		{0,0, alpha*f/sqrt(g)},
		{0,0,2*alpha/sqrt(g)},
		{alpha/sqrt(g),0,0},
	}
	self.eigenvectors[i] = {
		{f,			0,	f		},
		{2,			1,	2		},
		{-sqrt(f),	0,	sqrt(f)	},
	}
	self.eigenvectorsInverse[i] = {
		{1/(2*f), 	0, -1/(2*sqrt(f))	},
		{-2/f, 		1, 0				},
		{1/(2*f), 	0, 1/(2*sqrt(f))	},
	}
	-- doesn't work completely well
	-- the 5D flux matrix has zeroes on the top two alpha & g rows, so that equates to this just fine
	-- but on the first two alpha & g columns it has nonzero values, 
	-- which means in the 3D system we are neglecting alpha & g's explicit contribution to the eigenbasis
	-- that *should* be made up for through the reformulation, which shows us a different flux matrix for the lower 3x3 portion
	-- however something is still being lost
	-- it looks like alpha graph is growing on its rhs due to whatever influences 
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
	self.dq_dts[i][3] = self.dq_dts[i][3] + alpha * (A * D - K * K) / g
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

