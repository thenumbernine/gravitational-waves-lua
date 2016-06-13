--[[

dU/dt + F(U) = 0 <= conservative form
dU/dt + A(U) dU/dx = 0
for A(U) = dF/dU
dF/dU = A(U) = RU Lambda LU	<- eigenvector decomposition

dU/dW dW/dt + dF/dU dU/dW dW/dx = 0
dW/dt + dW/dU dF/dU dU/dW dW/dx = 0
dW/dt + A(W) dW/dx = 0	<= quasilinear form
so A(W) = dW/dU dF/dU dU/dW
A(W) = dW/dU A(U) dU/dW
A(U) = dU/dW A(W) dW/dU
A(W) = RW Lambda LW

so RU = dU/dW RW <=> RW = dW/dU RU
and LU = LW dW/dU <=> LW = LU dU/dW

dU/dt - dF/dx = 0
dU/dt = (Fl - Fr)/dx <= integrate-by-parts
for F = A(U) (UR + UL)/2 - RU Lambda (sgn(Lambda) + phi(rTilde) (Lambda * dt/dx - sgn(Lambda)) LU (UR - UL)/2
F = dU/dW [A(W) dW/dU (dU/dW WR + dU/dW WL)/2 - dW/dU RW Lambda (sgn(Lambda) + phi(rTilde) (Lambda * dt/dx - sgn(Lambda)) LW dW/dU (dU/dW WR - dU/dW WL)/2]
F = dU/dW [A(W) (WR + WL)/2 - RW Lambda (sgn(Lambda) + phi(rTilde) (Lambda * dt/dx - sgn(Lambda)) LW (WR - WL)/2]
F = dU/dW FW
for FW = A(W) (WR + WL)/2 - RW Lambda (sgn(Lambda) + phi(rTilde) (Lambda * dt/dx - sgn(Lambda)) LW (WR - WL)/2

dU/dW dW/dt = (Fl - Fr)/dx
dU/dW dW/dt = (dU/dW FWl - dU/dW FWr)/dx
dW/dt = (FWl - FWr)/dx

... except the dU/dW's have to be evaluated at the correct left and right locations ...
the flux dU/dW is at the interface, the left-and-right are at the cell center ...
F_i+1/2 = dU/dW_i+1/2 [A(U_i+1/2) dW/dU_i+1/2 (dU/dW_i+1 W_i+1 + dU/dW_i W_i)/2 - dW/dU_i+1/2 RW_i+1/2 Lambda (sgn(Lambda) + phi(rTilde) (Lambda * dt/dx - sgn(Lambda)) LW_i+1/2 dW/dU_i+1/2 (dU/dW_i+1 W_i+1 - dU/dW_i W_i)/2]

--]]
local table = require 'ext.table'
local class = require 'ext.class'
local Equation = require 'equation'

local Euler1DQuasiLinear = class(Equation)
Euler1DQuasiLinear.name = 'Euler 1D QuasiLinear'
Euler1DQuasiLinear.numStates = 3
Euler1DQuasiLinear.gamma = 5/3	

do
	local q = function(self, i) return self.qs[i] end
	local gamma = function(self) return self.equation.gamma end
	local rho = q:_(1)
	local vx = q:_(2)
	local P = q:_(3)
	local mx = rho*vx
	local eKin = .5 * vx * vx
	local EKin = rho * eKin
	local EInt = P / (gamma - 1)
	local ETotal = EKin + EInt
	local S = P / rho^gamma
	local H = EInt + P 
	local h = H / rho 
	local HTotal = ETotal + P 
	local hTotal = HTotal / rho
	
	Euler1DQuasiLinear:buildGraphInfos{
		-- prims
		{rho = rho},
		{vx = vx},
		{P = P},
		-- other vars
		{EInt = EInt},
		{H = H},
		{S = S},
		-- conservative
		{mx = mx},
		{ETotal = ETotal},
	}
end

function Euler1DQuasiLinear:initCell(sim,i)
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
	--[[ Sod
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
	
	-- [[ Sod
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

	return {rho, vx, P}
end

-- This is for an ideal gas
-- Euler1DQuasiLinear uses rho, vx, P as primitives (because Euler1DQuasiLinear is commonly used with ideal gasses, which relate P and eInt.  And P is more intuitive than eInt.)
-- SRHD uses rho, vx, eInt (because P is a function of eInt in general) 
-- maybe I should switch to rho, vx, eInt for generalization's sake ...
function Euler1DQuasiLinear:calcConsFromPrim(rho, vx, P)
	return rho, rho * vx, P/(self.gamma-1) + .5 * rho * vx*vx
end

-- rho is th density
-- mx is momentum in x direction (densitized velocity)
-- ETotal is densitized total energy
function Euler1DQuasiLinear:calcPrimFromCons(rho, mx, ETotal)
	local vx = mx / rho
	local P = (self.gamma - 1) * (ETotal - .5 * rho * vx * vx)
	return rho, vx, P
end

-- calculates the roe-averaged variables (often primitives) used for calculating interface eigenbasis
function Euler1DQuasiLinear:calcRoeValues(qL, qR)
	local gamma = self.gamma

	local rhoL, vxL, PL = table.unpack(qL)
	local hTotalL = self:calc_hTotal(rhoL, vxL, PL) 
	local sqrtRhoL = math.sqrt(rhoL)

	local rhoR, vxR, PR = table.unpack(qR)
	local hTotalR = self:calc_hTotal(rhoR, vxR, PR)
	local sqrtRhoR = math.sqrt(rhoR)

	local rho = sqrtRhoL * sqrtRhoR
	local vx = (sqrtRhoL * vxL + sqrtRhoR * vxR) / (sqrtRhoL + sqrtRhoR)
	local hTotal = (sqrtRhoL * hTotalL + sqrtRhoR * hTotalR) / (sqrtRhoL + sqrtRhoR)
	
	local a = self:calcSpeedOfSound_from_hTotal(vx, hTotal)
	return rho, vx, a
end

function Euler1DQuasiLinear:calc_hTotal(rho, vx, P)
	local gamma = self.gamma
	return .5 * vx * vx + gamma * P / (rho * (gamma - 1))
end

function Euler1DQuasiLinear:calcSpeedOfSound_from_hTotal(vx, hTotal)
	return math.sqrt((self.gamma - 1) * (hTotal - .5 * vx * vx))
end

function Euler1DQuasiLinear:calcSpeedOfSoundFromPrim(rho, P)
	return math.sqrt(self.gamma * P / rho)
end

function Euler1DQuasiLinear:calcEigenvalues(vx, a)
	return vx - a, vx, vx + a
end

Euler1DQuasiLinear.calcConsFromState = Euler1DQuasiLinear.calcConsFromPrim

function Euler1DQuasiLinear:calcEigenvaluesFromState(rho, vx, P)
	local a = self:calcSpeedOfSoundFromPrim(rho, P)
	return self:calcEigenvalues(vx, a)
end

function Euler1DQuasiLinear:calcMinMaxEigenvaluesFromState(rho, vx, P)
	local a = self:calcSpeedOfSoundFromPrim(rho, P)
	return vx - a, vx + a
end

-- used by HLL
-- what is the quasilinear flux vector?
-- if the quasilinear jacobian is dW/dU.dF/dW, then what does the new F vector become?
function Euler1DQuasiLinear:calcFluxForState(q)
	local gamma = self.gamma
	local U = self:calcConsFromState(table.unpack(q))
	-- conservative form flux:
	return 
		U[2],
		(gamma - 1) * U[3] + (3 - gamma) / 2 * U[2] * U[2] / U[1],
		gamma * U[2] * U[3] / U[1] + (1 - gamma) / 2 * U[2] * U[2] * U[2] / (U[1] * U[1])
end

function Euler1DQuasiLinear:calcEigenBasis(rho, vx, a, F, lambda, evL, evR)
	local gamma = self.gamma
	local aSq = a * a
	local vxSq = vx * vx
	
	fill(F[1], vx, rho, 0)
	fill(F[2], 0, vx, 1/rho)
	fill(F[3], 0, rho*aSq, vx)
	
	fill(lambda, self:calcEigenvalues(vx, a))
	
	fill(evR[1], 1, 1, 1)
	fill(evR[2], -a/rho, 0, a/rho)
	fill(evR[3], aSq, 0, aSq)
	
	fill(evL[1], 0, -.5*rho/a, .5/aSq)
	fill(evL[2], 1, 0, -1/aSq)
	fill(evL[3], 0, .5*rho/a, .5/aSq)
end

-- functions that use sim:

-- used by HLL
-- TODO how often do we create new tables of this?
function Euler1DQuasiLinear:calcInterfaceEigenvalues(sim, i, qL, qR)
	local rho, vx, a = self:calcRoeValues(qL, qR)
	fill(sim.eigenvalues[i], self:calcEigenvalues(vx, a))	
end

-- used by Roe
function Euler1DQuasiLinear:calcInterfaceEigenBasis(sim,i,qL,qR)
	local rho, vx, a = self:calcRoeValues(qL, qR)
	local F = sim.fluxMatrix[i] 
	local lambda = sim.eigenvalues[i]
	local evL = sim.eigenvectorsInverse[i]
	local evR = sim.eigenvectors[i]
	self:calcEigenBasis(rho, vx, a, F, lambda, evL, evR) 
end

return Euler1DQuasiLinear
