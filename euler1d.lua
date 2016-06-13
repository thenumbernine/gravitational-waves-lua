--[[
rho,t = -m,x
m,t = -(m^2/rho + p),x
E,t = -(m/rho*(E+p)),x

-- flux matrix:
[rho],t   [0											1						0			] [rho],x	[ 0 ]
[ m ],t + [(gamma-3)/2*m^2/rho^2						(3-gamma)*m/rho			gamma-1		] [ m ],x = [ 0 ]
[ E ],t   [(-gamma*m*E/rho^2 + (gamma-1)*m^3/rho^3)		H+(1-gamma)*m^2/rho^2	gamma*m/rho	] [ E ],x	[ 0 ]
--]]
local table = require 'ext.table'
local class = require 'ext.class'
local Equation = require 'equation'

local Euler1D = class(Equation)
Euler1D.name = 'Euler 1D'
Euler1D.numStates = 3
Euler1D.gamma = 5/3	

do
	--[[ just to know i could ...
	local I = function(...) return ... end -- identity function
	local q = I._:bind(I:_'qs'):uncurry(2):swap()
	local gamma = I:_'equation':_'gamma'
	--]]
	-- [[
	local q = function(self, i) return self.qs[i] end
	local gamma = function(self) return self.equation.gamma end
	--]]
	local rho = q:_(1)
	local mx = q:_(2)
	local ETotal = q:_(3)
	local vx = mx/rho
	local eKin = .5 * vx * vx
	local EKin = rho * eKin
	local EInt = ETotal - EKin
	local eInt = EInt / rho
	local P = (gamma - 1) * EInt
	local S = P / rho^gamma
	local H = EInt + P 
	local h = H / rho 
	local HTotal = ETotal + P 
	local hTotal = HTotal / rho
	
	Euler1D:buildGraphInfos{
		-- prims
		{rho = rho},
		{vx = vx},
		{P = P},
		-- other vars
		{EInt = EInt},
		{EKin = EKin},
		{H = H},
		{S = S},
		-- conservative
		--{mx = mx},
		{ETotal = ETotal},
	}
end

function Euler1D:initCell(sim,i)
	local x = sim.xs[i]
	--[[ constant
	local rho = 1
	local vx = 0
	local P = 1
	--]]
	--[[ linear
	local rho = 2 + x
	local vx = 0
	local P = 1
	--]]
	--[[ gaussian curve
	local sigma = 1/math.sqrt(10)
	local rho = math.exp(-x^2/sigma^2) + .1
	-- drho/dx = (-2x/sigma^2) exp(-x^2/sigma^2)
	local vx = 0
	local P = 1 + .1 * (math.exp(-x^2/sigma^2) + 1) / ((self.gamma - 1) * rho)
	--]]
	--[[ rarefaction wave
	local delta = .1
	local rho = 1	--x<0 and .2 or .8
	local vx = x<0 and .5-delta or .5+delta
	local P = 1
	--]]
	-- [[ Sod
	local rho = x < 0 and 1 or .125
	local vx = 0
	local P = x < 0 and 1 or .1
	--]]
	--[[ shock wave from numerical srhd paper marti & muller 2008
	local rho = 1
	local vx = x < 0 and .5 or 0
	local P = x < 0 and 1e+3 or 1
	--]]
	--[[ Sedov
	local rho = 1
	local vx = 0
	local P = i == math.floor(sim.gridsize/2) and 1e+3 or 1
	--]]
	

	-- matching the SRHD results
	
	--[[ Sod
	self.gamma = 7/5
	local rho = x < 0 and 1 or .125
	local vx = 0
	local P = x < 0 and 1 or .1
	--]]
	--[[ relativistic blast wave interaction
	self.gamma = 7/5
	local lhs = .9 * sim.domain.xmin + .1 * sim.domain.xmax
	local rhs = .1 * sim.domain.xmin + .9 * sim.domain.xmax
	local rho = 1
	local vx = 0
	local P = x < lhs and 1000 or (x > rhs and 100 or .01)
	--]]
	--[[ relativistic blast wave test problem 1
	self.gamma = 5/3
	local rho = x < 0 and 10 or 1
	local vx = 0
	local P = (self.gamma - 1) * rho * (x < 0 and 2 or 1e-6)
	--]]

	return {self:calcConsFromPrim(rho, vx, P)}
end

-- This is for an ideal gas
-- Euler1D uses rho, vx, P as primitives (because Euler1D is commonly used with ideal gasses, which relate P and eInt.  And P is more intuitive than eInt.)
-- SRHD uses rho, vx, eInt (because P is a function of eInt in general) 
-- maybe I should switch to rho, vx, eInt for generalization's sake ...
function Euler1D:calcConsFromPrim(rho, vx, P)
	return rho, rho * vx, P/(self.gamma-1) + .5 * rho * vx*vx
end

-- rho is th density
-- mx is momentum in x direction (densitized velocity)
-- ETotal is densitized total energy
function Euler1D:calcPrimFromCons(rho, mx, ETotal)
	local vx = mx / rho
	local P = (self.gamma - 1) * (ETotal - .5 * rho * vx * vx)
	return rho, vx, P
end

-- calculates the roe-averaged variables (often primitives) used for calculating interface eigenbasis
function Euler1D:calcRoeValues(qL, qR)
	local gamma = self.gamma

	local ETotalL = qL[3]
	local rhoL, vxL, PL = self:calcPrimFromCons(table.unpack(qL))
	local hTotalL = self:calc_hTotal(rhoL, PL, ETotalL)
	local sqrtRhoL = math.sqrt(rhoL)
	
	local ETotalR = qR[3]
	local rhoR, vxR, PR = self:calcPrimFromCons(table.unpack(qR))
	local hTotalR = self:calc_hTotal(rhoR, PR, ETotalR)
	local sqrtRhoR = math.sqrt(rhoR)
	
	local rho = sqrtRhoL * sqrtRhoR
	local vx = (sqrtRhoL * vxL + sqrtRhoR * vxR) / (sqrtRhoL + sqrtRhoR)
	local hTotal = (sqrtRhoL * hTotalL + sqrtRhoR * hTotalR) / (sqrtRhoL + sqrtRhoR)
	local Cs = self:calcSpeedOfSound(vx, hTotal)

	return rho, vx, hTotal, Cs
end

function Euler1D:calcEigenvalues(vx, Cs)
	return vx - Cs, vx, vx + Cs
end

function Euler1D:calcSpeedOfSound(vx, hTotal)
	return math.sqrt((self.gamma - 1) * (hTotal - .5 * vx * vx))
end

function Euler1D:calc_hTotal(rho, P, ETotal)
	return (ETotal + P) / rho
end

function Euler1D:calcEigenvaluesFromCons(rho, mx, ETotal)
	local rho, vx, P = self:calcPrimFromCons(rho, mx, ETotal)
	local hTotal = self:calc_hTotal(rho, P, ETotal)
	local Cs = self:calcSpeedOfSound(vx, hTotal)
	return self:calcEigenvalues(vx, Cs)
end

function Euler1D:calcConsFromState(...)
	return ...
end

-- used by HLL
function Euler1D:calcFluxForState(q)
	local gamma = self.gamma
	return 
		q[2],
		(gamma - 1) * q[3] + (3 - gamma) / 2 * q[2] * q[2] / q[1],
		gamma * q[2] * q[3] / q[1] + (1 - gamma) / 2 * q[2] * q[2] * q[2] / (q[1] * q[1])
end

function Euler1D:calcEigenBasisWrtPrims(rho, vx, hTotal, Cs, F, lambda, evL, evR)
	Cs = Cs or self:calcSpeedOfSound(vx, hTotal)
	local gamma = self.gamma

	fill(lambda, self:calcEigenvalues(vx, Cs))

	local CsSq = Cs * Cs

	if F then
		F[1][1] = vx
		F[1][2] = rho
		F[1][3] = 0
		F[2][1] = 0
		F[2][2] = vx
		F[2][3] = 1/rho
		F[3][1] = 0
		F[3][2] = rho * CsSq
		F[3][3] = vx
	end

	evL[1][1] = 1
	evL[1][2] = rho / (2 * Cs)
	evL[1][3] = -rho / (2 * Cs)
	evL[2][1] = 0
	evL[2][2] = 1/2
	evL[2][3] = 1/2
	evL[3][1] = 0
	evL[3][2] = rho * Cs / 2
	evL[3][3] = -rho * Cs / 2

	evR[1][1] = 1
	evR[2][1] = 0
	evR[3][1] = 0
	evR[1][2] = 0
	evR[2][2] = 1
	evR[3][2] = 1
	evR[1][3] = -1/CsSq
	evR[2][3] = 1 / (rho * Cs)
	evR[3][3] = -1 / (rho * Cs)
end

function Euler1D:calcEigenBasis(rho, vx, hTotal, Cs, F, lambda, evL, evR)
	Cs = Cs or self:calcSpeedOfSound(vx, hTotal)
	
	local gamma = self.gamma
	local CsSq = Cs * Cs
	local vxSq = vx * vx

	fill(F[1], 0, 1, 0)
	fill(F[2], .5*(gamma-3)*vxSq, (3-gamma)*vx, gamma-1)
	fill(F[3], vx*(.5*(gamma-1)*vxSq - hTotal), hTotal-(gamma-1)*vxSq, gamma*vx)
	
	fill(lambda, self:calcEigenvalues(vx, Cs))

	fill(evR[1], 1, 1, 1)
	fill(evR[2], vx - Cs, vx, vx + Cs)
	fill(evR[3], hTotal - Cs * vx, .5 * vxSq, hTotal + Cs * vx)

	fill(evL[1],
		(.5 * (gamma - 1) * vxSq + Cs * vx) / (2 * CsSq),
		 -(Cs + (gamma - 1) * vx) / (2 * CsSq),
		 (gamma - 1) / (2 * CsSq))
	fill(evL[2], 
		1 - (gamma - 1) * vxSq / (2 * CsSq),
		(gamma - 1) * vx / CsSq,
		-(gamma - 1) / CsSq)
	fill(evL[3],
		(.5 * (gamma - 1) * vxSq - Cs * vx) / (2 * CsSq),
		(Cs - (gamma - 1) * vx) / (2 * CsSq),
		(gamma - 1) / (2 * CsSq))
end

-- functions that use sim:

-- used by HLL
-- TODO how often do we create new tables of this?
function Euler1D:calcInterfaceEigenvalues(sim, i, qL, qR)
	local rho, vx, hTotal, Cs = self:calcRoeValues(qL, qR)
	fill(sim.eigenvalues[i], self:calcEigenvalues(vx, Cs))
end

-- used by Roe
function Euler1D:calcInterfaceEigenBasis(sim,i,qL,qR)
	local rho, vx, hTotal, Cs = self:calcRoeValues(qL, qR)
	local F = sim.fluxMatrix[i] 
	local lambda = sim.eigenvalues[i]
	local evL = sim.eigenvectorsInverse[i]
	local evR = sim.eigenvectors[i]
	self:calcEigenBasis(rho, vx, hTotal, Cs, F, lambda, evL, evR) 
end

return Euler1D
