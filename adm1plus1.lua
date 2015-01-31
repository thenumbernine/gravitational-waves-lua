#!/usr/bin/env luajit

require 'ext'

-- utility functions:

-- usage:
-- 	lookup(something)
-- 	lookup(something, where)
-- returns a list of locations of where to find that something
local function lookup(what, where, prefix, checked, found)
	assert(what ~= nil, "have to look for something")
	where = where or _G
	checked = checked or {}
	found = found or {}
	for k,v in pairs(where) do
		if v == what then	-- test without __eq? rawequals? 
			table.insert(found, prefix and prefix..'.'..k or k)
		end
	end
	for k,v in pairs(where) do
		if type(v) == 'table' then
			if not checked[v] then
				checked[v] = true	-- prevent recursive loops
				lookup(what, where, prefix..'.'..k, checked)
			end
		end
	end
	return unpack(found)
end



-- meta __index operation -- equivalent of operator[]() in C++
local function index(k,v) return k[v] end


-- makes math easier!
for k,v in pairs(math) do _G[k] = v end


-- update loop

local function class(...)
	local classobj = {}
	for _,parent in ipairs{...} do
		for k,v in pairs(parent) do
			classobj[k] = v
		end
	end
	classobj.__index = classobj
	classobj.super = ...
	setmetatable(classobj, {
		__call = function(classobj, ...)
			local obj = setmetatable({}, classobj)
			if obj.init then return obj, obj:init(...) end
			return obj
		end,
	})
	return classobj
end



local slopeLimiters = {
	donorCell = function(r) return 0 end,
	LaxWendroff = function(r) return 1 end,
	MC = function(r) return max(0, min(2, .5 * (1 + r), 2 * r)) end,
	superbee = function(r) return max(0, min(1, 2*r), min(2,r)) end,
	BeamWarming = function(r) return r end,
	Fromm = function(r) return .5 * (1 + r) end,
	vanLeer = function(r) return (r + abs(r)) / (1 + abs(r)) end,
}

local boundaryMethods = {}

boundaryMethods.mirror = function(self)
	for i=1,self.numStates do
		self.qs[1][i] = self.qs[4][i]
		self.qs[2][i] = self.qs[3][i]
		self.qs[self.gridsize-1][i] = self.qs[self.gridsize-2][i]
		self.qs[self.gridsize][i] = self.qs[self.gridsize-3][i]
	end
	-- ... and negative the wavespeed variable ... which varies per-equation 
	self.qs[1][2] = -self.qs[4][2]
	self.qs[2][2] = -self.qs[3][2]
	self.qs[self.gridsize-1][2] = -self.qs[self.gridsize-2][2]
	self.qs[self.gridsize][2] = -self.qs[self.gridsize-3][2]
end

boundaryMethods.freeFlow = function(self)
	for i=1,self.numStates do
		self.qs[1][i] = self.qs[3][i]
		self.qs[2][i] = self.qs[3][i]
		self.qs[self.gridsize][i] = self.qs[self.gridsize-2][i]
		self.qs[self.gridsize-1][i] = self.qs[self.gridsize-2][i]
	end
end

local calcFluxSchemes = {}

calcFluxSchemes.Roe = function(self)
	-- Roe solver:
	-- 1) calculate eigenbasis at interface with Roe-weighted average of states
	for i=2,self.gridsize do
		self:calcInterfaceEigenBasis(i)

		-- collect error for reporting
		local eigenbasisError = 0
		local fluxMatrixError = 0
		for j=1,self.numStates do
			for k=1,self.numStates do
				local eigenbasisSum = 0
				local fluxMatrixSum = 0
				for l=1,self.numStates do
					eigenbasisSum = eigenbasisSum + self.eigenvectors[i][j][l] * self.eigenvectorsInverse[i][l][k]
					fluxMatrixSum = fluxMatrixSum + self.eigenvectors[i][j][l] * self.eigenvectorsInverse[i][l][k] * self.eigenvalues[i][l]
				end
				eigenbasisError = eigenbasisError + abs(eigenbasisSum - (j == k and 1 or 0))
				fluxMatrixError = fluxMatrixError + abs(fluxMatrixSum - self.fluxMatrix[i][j][k])
			end
		end
		self.eigenbasisErrors[i] = eigenbasisError
		self.fluxMatrixErrors[i] = fluxMatrixError
	end

	-- determine timestep based on eigenvalue interfaces
	local dt
	if self.fixed_dt then
		dt = self.fixed_dt
	else
		local result = huge
		for i=1,self.gridsize do
			local eigenvaluesL = self.eigenvalues[i]
			local eigenvaluesR = self.eigenvalues[i+1]
			local maxLambda = max(0, unpack(eigenvaluesL))
			local minLambda = min(0, unpack(eigenvaluesR))
			local dx = self.ixs[i+1] - self.ixs[i]
			local dum = dx / (abs(maxLambda - minLambda) + 1e-9)
			result = min(result, dum)
		end
		dt = result * self.cfl
	end

	-- 2) calculate interface state difference in eigenbasis coordinates
	for i=2,self.gridsize do
		for j=1,self.numStates do
			local s = 0
			for k=1,self.numStates do
				s = s + self.eigenvectorsInverse[i][j][k] * (self.qs[i][k] - self.qs[i-1][k])
			end
			self.deltaQTildes[i][j] = s
		end
	end
	
	-- 3) slope limit on interface difference
	-- 4) transform back
	for i=2,self.gridsize do
		local fluxTilde = {}
		for j=1,self.numStates do
			local rTilde
			if self.deltaQTildes[i][j] == 0 then
				rTilde = 0
			else
				if self.eigenvalues[i][j] >= 0 then
					rTilde = self.deltaQTildes[i-1][j] / self.deltaQTildes[i][j]
				else
					rTilde = self.deltaQTildes[i+1][j] / self.deltaQTildes[i][j]
				end
			end
			local phi = self.slopeLimiter(rTilde)
			local theta = self.eigenvalues[i][j] >= 0 and 1 or -1
			local dx = self.xs[i] - self.xs[i-1]
			local epsilon = self.eigenvalues[i][j] * dt / dx
			local deltaFluxTilde = self.eigenvalues[i][j] * self.deltaQTildes[i][j]
			fluxTilde[j] = -.5 * deltaFluxTilde * (theta + phi * (epsilon - theta))
		end
		for j=1,self.numStates do
			local s = 0
			for k=1,self.numStates do
				s = s + self.eigenvectors[i][j][k] * fluxTilde[k]
				-- using the flux matrix itself allows for complete reconstruction even in the presence of zero-self.eigenvalues
				s = s + self.fluxMatrix[i][j][k] * (self.qs[i-1][k] + self.qs[i][k]) * .5
			end
			self.fluxes[i][j] = s
		end
	end

	return dt
end




-- base functions
local Simulation = class()

function Simulation:init(args)
	args = args or {}
	self.gridsize = assert(args.gridsize)
	self.domain = assert(args.domain)
	self.boundaryMethod = assert(args.boundaryMethod)
	self.slopeLimiter = assert(args.slopeLimiter)

	self.calcFlux = calcFluxSchemes.Roe 
	self.t = 0
	self.cfl = .5
	self.xs = {}
	self.ixs = {}
	self.qs = {}
	self.deltaQTildes = {}
	self.fluxes = {}
	self.dq_dts = {}
	-- used by Roe
	self.fluxMatrix = {}
	self.eigenvalues = {}
	self.eigenvectors = {}
	self.eigenvectorsInverse = {}
	self.eigenbasisErrors = {}
	self.fluxMatrixErrors = {}
end

function Simulation:reset()
	for i=1,self.gridsize do
		self.xs[i] = (i-.5)/self.gridsize*(self.domain.xmax - self.domain.xmin) + self.domain.xmin
	end
	-- interface center of coordinate system
	for i=1,self.gridsize+1 do
		self.ixs[i] = (i-1)/self.gridsize*(self.domain.xmax - self.domain.xmin) + self.domain.xmin
	end

	-- state at cell centers
	for i=1,self.gridsize do
		self.qs[i] = self:initCell(i)
	end

	-- state interfaces
	for i=1,self.gridsize+1 do
		self.deltaQTildes[i] = {}
		self.fluxes[i] = {}
		for j=1,self.numStates do
			self.deltaQTildes[i][j] = 0
			self.fluxes[i][j] = 0
		end
	end

	-- integration vars
	for i=1,self.gridsize do
		self.dq_dts[i] = {}
		for j=1,self.numStates do
			self.dq_dts[i][j] = 0
		end
	end
	
	for i=1,self.gridsize+1 do
		self.fluxMatrix[i] = {}
		self.eigenvalues[i] = {}
		self.eigenvectors[i] = {}
		self.eigenvectorsInverse[i] = {}
		for j=1,self.numStates do
			self.fluxMatrix[i][j] = {}
			self.eigenvectors[i][j] = {}
			self.eigenvectorsInverse[i][j] = {}
		end
		self.eigenbasisErrors[i] = 0
		self.fluxMatrixErrors[i] = 0
	end
end

function Simulation:zeroDeriv()
	-- zero deriv
	for i=1,self.gridsize do
		for j=1,self.numStates do
			self.dq_dts[i][j] = 0
		end
	end
end

function Simulation:addSourceToDeriv()
	for i=1,self.gridsize do
		self:addSourceToDerivCell(i)
	end
end

function Simulation:addFluxToDeriv()
	-- 5) integrate
	for i=1,self.gridsize do
		local dx = self.ixs[i+1] - self.ixs[i]
		for j=1,self.numStates do
			self.dq_dts[i][j] = self.dq_dts[i][j] - (self.fluxes[i+1][j] - self.fluxes[i][j]) / dx
		end
	end
end

function Simulation:integrateDeriv(dt)
	for i=1,self.gridsize do
		for j=1,self.numStates do
			self.qs[i][j] = self.qs[i][j] + dt * self.dq_dts[i][j]	
		end
	end
	self.t = self.t + dt
end

function Simulation:iterate()
	
	self:boundaryMethod()
	
	local dt = self:calcFlux()

	self:zeroDeriv()	
	self:addSourceToDeriv()
	self:addFluxToDeriv()
	self:integrateDeriv(dt)
end




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




--[[
hyperbolic formalism:

state vector: [alpha g A D K]	<- even though alpha and g are already represented by A and D ...
fluxes vector: alpha K * [0, 0, f/g, 1, A/K]
source vector: alpha / g * [-alhpa f K, -2 K g, 0, 0, A D - K^2]

alpha,t + (0),x = -alpha^2 f K / g
g,t + (0),x = -2 alpha K
A,t + (alpha K f / g),x = 0
D,t + (alpha K),x = 0
K,t + (alpha A),x = alpha / g (A D - K^2)

d_t alpha + (0),x = -alpha^2 f K / g
d_t g + (0),x = -2 alpha K
d_t A + K f / g alpha,x + alpha K / g f,x + alpha f / g K,x - alpha K f / g^2 g,x = 0 
d_t D + K alpha,x + alpha K,x = 0
d_t K + A alpha,x + alpha A,x = alpha / g (A D - K^2)

    [alpha]   [   0,            0,          0,   0,      0     ]     [alpha]   [  -alpha^2 f K / g   ]
    [  g  ]   [   0,            0,          0,   0,      0     ]     [  g  ]   [     -2 alpha K      ]
d_t [  A  ] + [f K / g, -alpha f K / g^2,   0,   0, alpha f / g] d_x [  A  ] = [          0          ]
    [  D  ]   [   K,            0,          0,   0,    alpha   ]     [  D  ]   [          0          ]
    [  K  ]   [   A,            0,        alpha, 0,      0     ]     [  K  ]   [alpha / g (A D - K^2)]

eigenvalues of A:
-lambda * (
	-lambda * (
		-lambda * det(
			-lambda		alpha f / g
			alpha		-lamda
		)
	)
) = 0
<=>
lambda = 0 has multiplicity 3,

lambda^2 - alpha^2 f / g = 0 
lambda = +-alpha sqrt(f / g)

for lambda = -alpha sqrt(f/g) the eigenvector is 
	[0, 0, 1, g/f, -sqrt(g/f)] dot state = A + D g/f - K sqrt(g/f)
	= A + D f / g - K * sqrt(f / g)
for lambda = +alpha sqrt(f/g) the eigenvector is 
	[0, 0, 1, g/f, +sqrt(g/f)] dot state = A + D g/f + K sqrt(g/f)
for lambda = 0 the eigenvectors are
[alpha,0,-A,0,-K]
[0,0,0,1,0]
[0,1,0,0,0]

eigenvector matrix is :
[	0			alpha	0	0	0			]
[	0			0		0	1	0			]
[	f/g			-A		0	0	f/g			]
[	1			0		1	0	1			]
[	-sqrt(f/g)	-K		0	0	sqrt(f/g)	]
inverse:
[(sqrt(f)*g^(3/2)*A-f*g*K)/(2*alpha*f^(3/2)*sqrt(g)),0,g/(2*f),0,-sqrt(g)/(2*sqrt(f))]
[1/alpha,0,0,0,0]
[-(g*A)/(alpha*f),0,-g/f,1,0]
[0,1,0,0,0]
[(sqrt(f)*sqrt(g)*K+g*A)/(2*alpha*f),0,g/(2*f),0,sqrt(g)/(2*sqrt(f))]

/* [wxMaxima: input   start ] */
fluxMatrix : matrix(
[0,0,0,0,0],
[0,0,0,0,0],
[f*K/g, -alpha*f*K/g^2, 0,0, alpha*f/g],
[K,0,0,0,alpha],
[A,0,alpha,0,0]);
/* [wxMaxima: input   end   ] */
/* [wxMaxima: input   start ] */
load("eigen");
/* [wxMaxima: input   end   ] */
/* [wxMaxima: input   start ] */
assume(f>0,g>0,alpha>0);
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
result[2][1][1] * sqrt(f^2 + f*g + g^2)/g,
result[2][3][1] * sqrt(K^2 + A^2 + alpha^2),
result[2][3][2],
[0,1,0,0,0],
result[2][2][1] * sqrt(f^2 + f*g + g^2)/g))$
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


local ADM1D5VarSim = class(Simulation)
	
ADM1D5VarSim.numStates = 5 

function ADM1D5VarSim:init(...)
	Simulation.init(self, ...)

	self.graphInfos = {
		{viewport={0/3, 0/3, 1/3, 1/3}, getter=index:bind(self.qs):index(1), name='alpha', color={1,0,1}},
		{viewport={0/3, 1/3, 1/3, 1/3}, getter=index:bind(self.qs):index(3), name='A', color={0,1,0}},
		{viewport={1/3, 0/3, 1/3, 1/3}, getter=index:bind(self.qs):index(2), name='g', color={.5,.5,1}},
		{viewport={1/3, 1/3, 1/3, 1/3}, getter=index:bind(self.qs):index(4), name='D', color={1,1,0}},
		{viewport={2/3, 0/3, 1/3, 1/3}, getter=index:bind(self.qs):index(5), name='K', color={0,1,1}},
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

	function ADM1D5VarSim:initCell(i)
		local x = self.xs[i]
		-- primitives:
		local alpha = calc_alpha(x)
		local g = calc_g(x)
		-- state variables:
		local A = d_calc_alpha(x) / calc_alpha(x)
		local D = 1/2 * d_calc_g(x)
		local K = -d2_calc_h(x) / sqrt(calc_g(x))
		return {alpha, g, A, D, K}
	end
end

function ADM1D5VarSim:calcInterfaceEigenBasis(i)
	local avgQ = {}
	for j=1,self.numStates do
		avgQ[j] = (self.qs[i-1][j] + self.qs[i][j]) / 2
	end
	
	local alpha, g, A, D, K = unpack(avgQ)
	local x = self.ixs[i]
	local f = self.calc_f(alpha)
	local lambda = alpha * sqrt(f / g)		
	self.eigenvalues[i] = {-lambda, 0, 0, 0, lambda}
	-- row-major, math-indexed
	self.fluxMatrix[i] = {
		{0,0,0,0,0},
		{0,0,0,0,0},
		{f*K/g, -alpha*f*K/g^2, 0,0, alpha*f/g},
		{K,0,0,0,alpha},
		{A,0,alpha,0,0},
	}
	self.eigenvectors[i] = {
		{0,			alpha,	0,	0,	0			},	-- alpha
		{0,			0,		0,	1,	0			},	-- g
		{f/g,		-A,		0,	0,	f/g			},	-- A
		{1,			0,		1,	0,	1			},	-- D
		{-sqrt(f/g),-K,		0,	0,	sqrt(f/g)	},	-- K
	}
	self.eigenvectorsInverse[i] = {
		{(g * A / f - K * sqrt(g / f)) / (2 * alpha), 0, g / (2 * f), 0, -.5 * sqrt(g / f)}, 
		{1 / alpha, 0, 0, 0, 0}, 
		{-(g * A) / (alpha * f), 0, -g / f, 1, 0}, 
		{0, 1, 0, 0, 0}, 
		{(g * A / f + K * sqrt(g / f)) / (2 * alpha), 0, g / (2 * f), 0, .5 * sqrt(g / f)}, 
	}
	-- note that because we have zero eigenvalues that the eigendecomposition cannot reconstruct the flux matrix
end

function ADM1D5VarSim:addSourceToDerivCell(i)
	local alpha, g, A, D, K = unpack(self.qs[i])
	local f = self.calc_f(alpha)
	self.dq_dts[i][1] = self.dq_dts[i][1] - alpha * alpha * f * K / g
	self.dq_dts[i][2] = self.dq_dts[i][2] - 2 * alpha * K
	self.dq_dts[i][5] = self.dq_dts[i][5] + alpha * (A * D - K * K) / g
end

--[[
variables:
alpha, g_rr, g_hh, A_r, D_rrr, D_rhh, K_rr, K_hh, V_r

definitions:
A_r = (ln alpha),r = alpha,r / alpha
D_kij = g_ij,k/2
V_k = D_km^m - D^m_mk

equations:
alpha,t = -alpha^2 f tr K
g_rr,t = -2 alpha K_rr
g_hh,t = -2 alpha K_hh

A_r,t + (alpha f tr K),r = 0
D_rrr,t + (alpha K_rr),r = 0
D_rhh,t + (alpha K_hh),r = 0
K_rr,t + (alpha lambda^r_rr),r = alpha S_rr
K_hh,t + (alpha lambda^r_hh),r = alpha S_hh
V_r,t = alpha P_r

for:
tr K = g^rr K_rr g^hh K_hh = K_rr / g_rr + K_hh / g_hh
lambda^r_rr = A_r + 2 V_r - 2 D_rhh / g_hh
lambda^r_hh = D_rhh / g_rr
S_rr = K_rr(2 K_hh / g_hh - K_rr/g_rr) + A_r(D_rrr/g_rr - 2 D_rhh/g_hh) + 2 D_rhh/g_hh (D_rrr/g_rr - D_rhh/g_hh) + 2 A_r V_r
S_hh = K_rr K_hh/g_rr - D_rrr D_rhh / g_rr^2 + 1
P_r = -2 / g_hh (A_r K_hh - D_rhh (K_hh / g_hh - K_rr / g_rr))

constraint:
V_r = 2 D_rhh / g_hh

sphrical metric: ds^2 = -dt^2 + dr^2 + r^2 (dh^2 + sin(h)^2 dp^2)
g_uv = diag(-1, 1, r^2, r^2 sin(h)^2)
g^uv = diag(-1, 1, r^-2, r^-2 sin(h)^-2)
... h = theta, p = phi
g_tt is fixed?  or is derived?  is alpha^2?
g_pp is fixed?
either way, none of the skew elements exist so the inverse metric is the diagonal of inverse elements 

... evolution equations become ...
alpha,t = -alpha^2 f tr K
g_rr,t = -2 alpha K_rr
g_hh,t = -2 alpha K_hh
A_r,t 
	+ alpha,r (f tr K) 
	+ f,r (alpha tr K) 	... which equals ... alpha,r (df/dalpha alpha tr K)
	+ g_rr,r (-alpha f K_rr / g_rr^2) 
	+ g_hh,r (-alpha f K_hh / g_hh^2) 
	+ K_rr,r (alpha f / g_rr)
	+ K_hh,r (alpha f / g_hh)
	= 0
D_rrr,t + alpha,r K_rr + K_rr,r alpha = 0
D_rhh,t + alpha,r K_hh + K_hh,r alpha = 0
K_rr,t 
	+ alpha,r lambda^r_rr 
	+ g_hh,r (2 alpha D_rhh / g_hh^2)
	+ A_r,r alpha
	+ D_rhh,r (-2 alpha / g_hh)
	+ V_r,r (2 alpha)
	= alpha S_rr
K_hh,t 
	+ alpha,r lambda^r_hh 
	+ g_rr,r (-alpha / g_rr^2)
	+ D_rhh,r (alpha / g_rr)
	= alpha S_hh
V_r,t = alpha P_r

alpha, g_rr, g_hh, A_r, D_rrr, D_rhh, K_rr, K_hh, V_r
flux matrix:
[0,	0,	0,	0,	0,	0,	0,	0,	0]
[0,	0,	0,	0,	0,	0,	0,	0,	0]
[0,	0,	0,	0,	0,	0,	0,	0,	0]
[(f + df/dalpha alpha) tr K, -alpha f K_rr/g_rr^2, -alpha f K_hh/g_hh^2, 0, 0, 0, alpha f/g_rr, alpha f/g_hh, 0]
[K_rr, 0, 0, 0, 0, 0, alpha, 0, 0]
[K_hh, 0, 0, 0, 0, 0, 0, alpha, 0]
[lambda^r_rr, 0, 2 alpha D_rhh/g_hh^2, alpha, 0, -2 alpha/g_hh, 0, 0, 2 alpha]
[lambda^r_hh, -alpha/g_rr^2, 0, 0, 0, alpha/g_rr, 0, 0, 0]
[0,	0,	0,	0,	0,	0,	0,	0,	0]

eigenvalues:
-alpha*sqrt(f/g_rr), -alpha/sqrt(g_rr), 0 x5, alpha/sqrt(g_rr), alpha*sqrt(f/g_rr), 

eigenvector columns:
lambda = -alpha*sqrt(f/g_rr)
v = [0,0,0,f,g_rr,0,-sqrt(f g_rr),0,0]
v is equal to sqrt(g_rr) K_rr - D_rrr g_rr / sqrt(f) - A_r sqrt(f)
Alcubierre says sqrt(g_rr) K_hh - D_rhh == [0,0,0,0,0,-1,0,sqrt(g_rr),0]
so does D_rrr g_rr / sqrt(f) + A_r sqrt(f) == D_rhh ?
	does g_rr,r g_rr / sqrt(f) + 2 alpha,r sqrt(f) / alpha == g_hh,r
	does g_rr,r (g_rr / sqrt(f) - 1) + 2 alpha,r sqrt(f) / alpha == 0?
also which of many forms should Alcubierre's eqn be reformulated into a vector?

lambda = -alpha/sqrt(g_rr)
v = [0,0,0,f,-(f-2)*g_rr,(f-1)*g_hh,(f-2)*sqrt(g_rr),-((f-1)*g_hh)/sqrt(g_rr),0]
Alcubierre says sqrt(f g_rr) tr K - (A_r + 2 V_r) = [0, 0, 0, -1, 0, 0, sqrt(f/g_rr), sqrt(f g_rr) / g_hh, -2]

lambda = 0
Alcubierre says alpha, g_rr, g_hh, V_r, A_r - f D_r^m_m

lambda = alpha/sqrt(g_rr)
v = [0,0,0,f,-(f-2)*g_rr,(f-1)*g_hh,-(f-2)*sqrt(g_rr),((f-1)*g_hh)/sqrt(g_rr),0]
Alcubierre says says sqrt(f g_rr) tr K + (A_r + 2 V_r)

lambda = alpha*sqrt(f/g_rr)
v = [0,0,0,f,g_rr,0,sqrt(f g_rr),0,0]
... which amounts to A_r sqrt(f) + D_rrr g_rr / sqrt(f) + K_rr sqrt(g_rr)
Alcubierre says sqrt(g_rr) K_hh + D_rhh == [0,0,0,0,0,1,0,sqrt(g_rr),0]



--]]
local ADM2DSim = class(Simulation)
	
ADM2DSim.numStates = 9

-- initial conditions
function ADM2DSim:init(...)
	Simulation.init(self, ...)

	self.graphInfos = {
		{viewport={0/3, 0/3, 1/3, 1/3}, getter=index:bind(self.qs):index(1), name='alpha', color={1,0,1}},
		{viewport={0/3, 1/3, 1/3, 1/3}, getter=index:bind(self.qs):index(3), name='A', color={0,1,0}},
		{viewport={1/3, 0/3, 1/3, 1/3}, getter=index:bind(self.qs):index(2), name='g', color={.5,.5,1}},
		{viewport={1/3, 1/3, 1/3, 1/3}, getter=index:bind(self.qs):index(4), name='D', color={1,1,0}},
		{viewport={2/3, 0/3, 1/3, 1/3}, getter=index:bind(self.qs):index(5), name='K', color={0,1,1}},
		{viewport={0/3, 2/3, 1/3, 1/3}, getter=log:compose(index:bind(self.eigenbasisErrors)), name='error', color={1,0,0}, range={-30, 30}},
		{viewport={1/3, 2/3, 1/3, 1/3}, getter=log:compose(index:bind(self.fluxMatrixErrors)), name='error', color={1,0,0}, range={-30, 30}},
	}
end

do
	local sigma = 10
	local xc = 150
	local H = 5
	-- aux var for init
	local function calc_h(r) return H * exp(-(r - xc)^2 / sigma^2) end
	local function d_calc_h(r) return -2 * (r - xc) / sigma^2 * calc_h(r) end
	local function d2_calc_h(r) return (-2 / sigma^2 + 4 * (r - xc)^2 / sigma^4) * calc_h(r) end
	local function calc_g_rr(r) return 1 - d_calc_h(r)^2 end
	local function dr_calc_g_rr(r) return -2 * d_calc_h(r) * d2_calc_h(r) end
	local function calc_g_hh(r) return r^2 end
	local function dr_calc_g_hh(r) return 2*r end
	local function calc_alpha(r) return 1 end
	local function dr_calc_alpha(r) return 0 end

	function ADM2DSim:initCell(i)
		local r = self.xs[i]
		local alpha = calc_alpha(r)
		local g_rr = calc_g_rr(r)
		local g_hh = calc_g_hh(r)
		local A_r = dr_calc_alpha(r)
		local D_rrr = dr_calc_g_rr(r)/2
		local D_rhh = dr_calc_g_hh(r)/2
		local K_rr = -d2_calc_h(r) / sqrt(g_rr)
		local K_hh = -r * d_calc_h(r) / sqrt(g_rr)
		local V_r = D_rhh / g_hh
		return {alpha, g_rr, g_hh, A_r, D_rrr, D_rhh, K_rr, K_hh, V_r}
	end
end

function ADM2DSim:calcInterfaceEigenBasis(i)
	local avgQ = {}
	for j=1,self.numStates do
		avgQ[j] = (self.qs[i-1][j] + self.qs[i][j]) / 2
	end
	
	local alpha, g_rr, g_hh, A_r, D_rrr, D_rhh, K_rr, K_hh, V_r = unpack(avgQ)
	local x = self.ixs[i]
	local f = self.calc_f(alpha)
	
	local l1 = alpha / sqrt(g_rr)
	local l2 = l1 * (f)
	
	local lambdaUr_rr = A_r + 2 * V_r - 2 * D_rhh / g_hh
	local lambdaUr_hh = D_rhh / g_rr
	
	self.eigenvalues[i] = {-l2, -l1, 0, 0, 0, 0, 0, l1, l2}
	-- row-major, math-indexed
	self.fluxMatrix[i] = {
		{0,0,0,0,0,0,0,0,0},
		{0,0,0,0,0,0,0,0,0},
		{0,0,0,0,0,0,0,0,0},
		{f * tr_K, -alpha * f * K_rr/g_rr^2, -alpha * f * K_hh/g_hh^2, 0, 0, 0, alpha * f / g_rr, alpha * f / g_hh, 0},
		{K_rr, 0, 0, 0, 0, 0, alpha, 0, 0},
		{K_hh, 0, 0, 0, 0, 0, 0, alpha, 0},
		{lambdaUr_rr, 0, 2 * alpha * D_rhh/g_hh^2, alpha, 0, -2 * alpha/g_hh, 0, 0, 2 * alpha},
		{lambdaUr_hh, -alpha/g_rr^2, 0, 0, 0, alpha/g_rr, 0, 0, 0},
		{0,0,0,0,0,0,0,0,0},
	}
	-- TODO figure these out
	self.eigenvectors[i] = {
		{0,			0,					1,	0,	0,	0,	0,	0,					0,			},
		{0,			0,					0,	1,	0,	0,	0,	0,					0,			},
		{0,			0,					0,	0,	1,	0,	0,	0,					0,			},
		{0,			-1,					0,	0,	0,	0,	0,	1,					0,			},
		{0,			0,					0,	0,	0,	1,	0,	0,					0,			},
		{-1,		0,					0,	0,	0,	0,	0,	0,					1,			},
		{0,			sqrt(f/g_rr),		0,	0,	0,	0,	0,	sqrt(f/g_rr),		0,			},
		{sqrt(g_rr),sqrt(f*g_rr)/g_hh,	0,	0,	0,	0,	0,	sqrt(f*g_rr)/g_hh,	sqrt(g_rr),	},
		{0,			-2,					0,	0,	0,	0,	0,	2,					0,			},
	}
	self.eigenvectorsInverse[i] = {
		{(g * A / f - K * sqrt(g / f)) / (2 * alpha), 0, g / (2 * f), 0, -.5 * sqrt(g / f)}, 
		{1 / alpha, 0, 0, 0, 0}, 
		{-(g * A) / (alpha * f), 0, -g / f, 1, 0}, 
		{0, 1, 0, 0, 0}, 
		{(g * A / f + K * sqrt(g / f)) / (2 * alpha), 0, g / (2 * f), 0, .5 * sqrt(g / f)}, 
	}
	-- note that because we have zero eigenvalues that the eigendecomposition cannot reconstruct the flux matrix
end

function ADM2DSim:addSourceToDerivCell(i)
	local alpha, g, A, D, K = unpack(self.qs[i])
	local f = self.calc_f(alpha)
	self.dq_dts[i][1] = self.dq_dts[i][1] - alpha * alpha * f * K / g
	self.dq_dts[i][2] = self.dq_dts[i][2] - 2 * alpha * K
	self.dq_dts[i][5] = self.dq_dts[i][5] + alpha * (A * D - K * K) / g
end




local EulerSim = class(Simulation)

EulerSim.numStates = 3
EulerSim.gamma = 5/3	

function EulerSim:init(...)
	Simulation.init(self, ...)
	
	--index:bind(qs) => f(k) = qs[k] takes 1 arg and returns an array of 3 elements
	--index:bind(qs)[1] => f(k) = qs[k][1] takes 1 arg and returns the 1st element of the 3
	self.graphInfos = {
		{viewport={0/3, 0/3, 1/3, 1/3}, getter=index:bind(self.qs):index(1), name='rho', color={1,0,1}},
		{viewport={1/3, 0/3, 1/3, 1/3}, getter=index:bind(self.qs):index(2) / index:bind(self.qs):index(1), name='u', color={0,1,0}},
		{viewport={2/3, 0/3, 1/3, 1/3}, getter=index:bind(self.qs):index(3) / index:bind(self.qs):index(1), name='E', color={.5,.5,1}},
		{viewport={0/3, 1/3, 1/3, 1/3}, getter=log:compose(index:bind(self.eigenbasisErrors)), name='error', color={1,0,0}, range={-30, 30}},
		{viewport={1/3, 1/3, 1/3, 1/3}, getter=log:compose(index:bind(self.fluxMatrixErrors)), name='error', color={1,0,0}, range={-30, 30}},
	}
end

function EulerSim:initCell(i)
	local rho = self.xs[i] < 0 and .1 or 1
	local u = 0
	local E = 1	+ .5 * u * u	-- internal + kinetic
	return {rho, rho * u, rho * E}
end

function EulerSim:calcInterfaceEigenBasis(i)
	local rhoL = self.qs[i-1][1]
	local uL = self.qs[i-1][2] / rhoL 
	local EL = self.qs[i-1][3] / rhoL
	local eIntL = EL - .5 * uL^2
	local PL = (self.gamma - 1) * rhoL * eIntL
	local HL = EL + PL / rhoL
	local weightL = sqrt(rhoL)

	local rhoR = self.qs[i-1][1]
	local uR = self.qs[i-1][2] / rhoR 
	local ER = self.qs[i-1][3] / rhoR
	local eIntR = ER - .5 * uR^2
	local PR = (self.gamma - 1) * rhoR * eIntR
	local HR = ER + PR / rhoR
	local weightR = sqrt(rhoR)

	local u = (weightL * uL + weightR * uR) / (weightL + weightR)
	local H = (weightL * HL + weightR * HR) / (weightL + weightR)
	local E = (weightL * EL + weightR * ER) / (weightL + weightR)
	
	local Cs = sqrt((self.gamma - 1) * (H - .5 * u^2))

	self.fluxMatrix[i] = {
		{0, 1, 0},
		{(self.gamma-3)/2*u^2, (3-self.gamma)*u, self.gamma-1},
		{-u*(self.gamma*E + (1-self.gamma)*u^2), H + (1-self.gamma) * u^2, self.gamma * u},
	}
	self.eigenvalues[i] = {u - Cs, u, u + Cs}
	self.eigenvectors[i] = {
		{1,				1,			1,			},
		{u - Cs,		u,			u + Cs,		},
		{H - Cs * u,	.5 * u^2,	H + Cs * u,	},
	}
	-- [[ symbolically
	self.eigenvectorsInverse[i] = {
		{
			(.5 * (self.gamma - 1) * u^2 + Cs * u) / (2 * Cs^2),
			-(Cs + (self.gamma - 1) * u) / (2 * Cs^2),
			(self.gamma - 1) / (2 * Cs^2),
		}, {
			1 - (self.gamma - 1) * u^2 / (2 * Cs^2),
			(self.gamma - 1) * u / Cs^2,
			-(self.gamma - 1) / Cs^2,
		}, {
			(.5 * (self.gamma - 1) * u^2 - Cs * u) / (2 * Cs^2),
			(Cs - (self.gamma - 1) * u) / (2 * Cs^2),
			(self.gamma - 1) / (2 * Cs^2),		
		}
	}
	--]]
	--[[ numerically via cramers rule
	local det = eigenvectors[i][1][1] * eigenvectors[i][2][2] * eigenvectors[i][3][3]
			+ eigenvectors[i][2][1] * eigenvectors[i][3][2] * eigenvectors[i][1][3]
			+ eigenvectors[i][3][1] * eigenvectors[i][1][2] * eigenvectors[i][2][3]
			- eigenvectors[i][3][1] * eigenvectors[i][2][2] * eigenvectors[i][1][3]
			- eigenvectors[i][2][1] * eigenvectors[i][1][2] * eigenvectors[i][3][3]
			- eigenvectors[i][1][1] * eigenvectors[i][3][2] * eigenvectors[i][2][3];
	if det == 0 then
		for j=1,3 do
			for k=1,3 do
				console.log('A('+i+','+j+') = '+eigenvectors[i][j][k]);
			end
		end
		error 'singular!'
	end
	local invDet = 1 / det
	for j=1,3 do
		local j1 = j % 3 + 1 
		local j2 = j1 % 3 + 1
		for k=1,3 do
			local k1 = k % 3 + 1
			local k2 = k1 % 3 + 1
			eigenvectorsInverse[i][k][j] = invDet * (eigenvectors[i][j1][k1] * eigenvectors[i][j2][k2] - eigenvectors[i][j1][k2] * eigenvectors[i][j2][k1])
		end
	end
	--]]
end

function EulerSim:addSourceToDerivCell(i) end




-- setup
-- [[
--local sim = ADM1D3VarSim{
local sim = ADM1D5VarSim{
	gridsize = 1200,
	domain = {xmin=0, xmax=300},
	boundaryMethod = boundaryMethods.freeFlow,	-- still reflecting despite freeflow ...
	slopeLimiter = slopeLimiters.donorCell,
}
sim.fixed_dt = nil	--.125 
sim.tmax = 70
	-- Bona-Masso slicing conditions:
--sim.calc_f = function(alpha) return 1 end 
--sim.calc_f = function(alpha) return 1.69 end 
--sim.calc_f = function(alpha) return .49 end 
sim.calc_f = function(alpha) return 1 + 1/(alpha*alpha) end 
--]]

--[[
local sim = EulerSim{
	gridsize = 200,
	domain = {xmin=-1, xmax=1},
	boundaryMethod = boundaryMethods.mirror,
	slopeLimiter = slopeLimiters.superbee,
}
--]]

sim:reset()


--[[ text
for iter=1,100 do
	sim:iterate()
	for j=1,sim.numStates do
		io.write(({'alpha','g','A','D','K'})[j])
		for i=1,sim.gridsize do
			io.write('\t', qs[i][j])
		end
		print()
	end
end
--]]

-- [=[ graphics
local GLApp = require 'glapp'
local gl = require 'ffi.OpenGL'
local sdl = require 'ffi.sdl'
local GLApp = require 'glapp.glapp'

local TestApp = class(GLApp)

TestApp.width = 640
TestApp.height = 640
TestApp.showFPS = false

function TestApp:init()
	TestApp.super.init(self)
	self.doIteration = false
end

function TestApp:event(event)
	if event.type == sdl.SDL_KEYDOWN then
		local callback = self.keyDownCallbacks[event.key.keysym.sym]
		if callback then
			callback(self)
		end
	end
end

TestApp.keyDownCallbacks = {
	[sdl.SDLK_r] = function(self)
		sim:reset()
	end,
	[sdl.SDLK_e] = function(self)
		self.reportError = not self.reportError	
	end,
	[sdl.SDLK_SPACE] = function(self)
		self.doIteration = not self.doIteration
	end,
	[sdl.SDLK_u] = function(self)
		self.doIteration = 'once'
	end,
}

function TestApp:update(...)
	gl.glClear(gl.GL_COLOR_BUFFER_BIT)

	if sim.tmax then
		local t = sim.t
		if self.oldt then
			if t >= sim.tmax and self.oldt < sim.tmax then
				self.doIteration = false
			end
		end
		self.oldt = t
	end

	if self.doIteration then
		if self.doIteration == 'once' then
			self.doIteration = false
		end
		sim:iterate()
	end
	
	local xs = sim.xs
	local w, h = self:size()
	for infoIndex,info in ipairs(sim.graphInfos) do
		local ys = {}
		local ymin, ymax
		for i=1,sim.gridsize do
			local y = info.getter(i)
			ys[i] = y
			if y == y and abs(y) < huge then
				if not ymin or y < ymin then ymin = y end
				if not ymax or y > ymax then ymax = y end
			end
		end
		if self.reportError then
			print(info.name, 'min', ymin, 'max', ymax)
		end
	
		if info.range then
			ymin, ymax = unpack(info.range)
		else
			if not ymin or not ymax or ymin ~= ymin or ymax ~= ymax then
				ymin = -1
				ymax = 1
			--elseif abs(ymin) == huge or abs(ymax) == huge then
			else
				local base = 10	-- round to nearest base-10
				local scale = 10 -- ...with increments of 10
				ymin = (ymin<0 and -1 or 1)*(abs(ymin)==huge and 1e+100 or base^((log(abs(1.1*ymin-.1*ymax),base)*scale)/scale))
				ymax = (ymax<0 and -1 or 1)*(abs(ymax)==huge and 1e+100 or base^((log(abs(1.1*ymax-.1*ymin),base)*scale)/scale))
			end
		end
		
		local xmin, xmax = sim.domain.xmin, sim.domain.xmax
		xmin, xmax = 1.1 * xmin - .1 * xmax, 1.1 * xmax - .1 * xmin	
		
		-- display
		-- TODO viewports per variable and maybe ticks too
		gl.glViewport(
			info.viewport[1] * w,
			(1 - info.viewport[2] - info.viewport[4]) * h,
			info.viewport[3] * w,
			info.viewport[4] * h)
		gl.glMatrixMode(gl.GL_PROJECTION)
		gl.glLoadIdentity()
		gl.glOrtho(xmin, xmax, ymin, ymax, -1, 1)
		gl.glMatrixMode(gl.GL_MODELVIEW)
		gl.glLoadIdentity()

-- [[
		gl.glColor3f(.25, .25, .25)
		local xrange = xmax - xmin
		local xstep = 10^floor(log(xrange, 10) - .5)
		gl.glBegin(gl.GL_LINES)
		for x=floor(xmin/xstep)*xstep,ceil(xmax/xstep)*xstep,xstep do
			gl.glVertex2f(x,ymin)
			gl.glVertex2f(x,ymax)
		end
		gl.glEnd()
		local yrange = ymax - ymin
		local ystep = 10^floor(log(yrange, 10) - .5)
		gl.glBegin(gl.GL_LINES)
		for y=floor(ymin/ystep)*ystep,ceil(ymax/ystep)*ystep,ystep do
			gl.glVertex2f(xmin,y)
			gl.glVertex2f(xmax,y)
		end
		gl.glEnd()
		gl.glColor3f(.5, .5, .5)
		gl.glBegin(gl.GL_LINES)
		gl.glVertex2f(xmin, 0)
		gl.glVertex2f(xmax, 0)
		gl.glVertex2f(0, ymin)
		gl.glVertex2f(0, ymax)
		gl.glEnd()
--]]

		gl.glColor3f(unpack(info.color))
		gl.glPointSize(2)
		for _,mode in ipairs{gl.GL_POINTS, gl.GL_LINE_STRIP} do
			gl.glBegin(mode)
			for i=1,sim.gridsize do
				gl.glVertex2f(xs[i], ys[i])
			end
			gl.glEnd()
		end
		gl.glPointSize(1)
	end
	if self.reportError then
		self.reportError = false
	end
	gl.glViewport(0,0,w,h)
end
TestApp():run()
--]=]

