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
	local q = function(self,i) return self.qs[i] end
	local gamma = function(self,i) return self.equation.gamma end
	local rho = q:_(1)
	local mom = q:_(2)
	local ETotal = q:_(3)
	local eTotal = ETotal / rho
	local vx = mom/rho
	local EKin = .5 * rho * vx^2
	local EInt = ETotal - EKin
	local P = (gamma - 1) * EInt
	local H = EInt * gamma
	local S = P / rho^gamma

	Euler1D:buildGraphInfos{
		-- prims
		{rho = rho},
		{vx = vx},
		{P = P},
		-- other vars
		{EInt = EInt},
		{H = H},
		{S = S},
		-- conservative
		{mom = mom},
		{ETotal = ETotal},
		-- roe scheme
		{['log eigenbasis error'] = function(self,i) return self.eigenbasisErrors and math.log(self.eigenbasisErrors[i]) end},
		{['log reconstruction error'] = function(self,i) return self.fluxMatrixErrors and math.log(self.fluxMatrixErrors[i]) end},
	}
end

function Euler1D:initCell(sim,i)
	local gamma = self.gamma
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
	local P = 1 + .1 * (math.exp(-x^2/sigma^2) + 1) / ((gamma - 1) * rho)
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
	return {rho, rho * vx, P/(gamma-1) + .5 * rho * vx^2}
end

-- used by HLL
function Euler1D:calcFluxForState(sim, i, q, flux)
	flux = flux or {}
	local gamma = self.gamma
	flux[1] = q[2]
	flux[2] = (gamma - 1) * q[3] + (3 - gamma) / 2 * q[2] * q[2] / q[1]
	flux[3] = gamma * q[2] * q[3] / q[1] + (1 - gamma) / 2 * q[2] * q[2] * q[2] / (q[1] * q[1])
	return flux
end

-- used by HLL
function Euler1D:calcInterfaceEigenvalues(sim, i, qL, qR, S)
	S = S or {}
	
	local gamma = self.gamma
	
	local rhoL = qL[1]
	local vxL = qL[2] / rhoL 
	local eTotalL = qL[3] / rhoL
	local eIntL = eTotalL - .5 * vxL^2
	local PL = (gamma - 1) * rhoL * eIntL
	local hTotalL = eTotalL + PL / rhoL
	local weightL = sqrt(rhoL)
	
	local rhoR = qR[1]
	local vxR = qR[2] / rhoR 
	local eTotalR = qR[3] / rhoR
	local eIntR = eTotalR - .5 * vxR^2
	local PR = (gamma - 1) * rhoR * eIntR
	local hTotalR = eTotalR + PR / rhoR
	local weightR = sqrt(rhoR)
	
	local vx = (weightL * vxL + weightR * vxR) / (weightL + weightR)
	local hTotal = (weightL * hTotalL + weightR * hTotalR) / (weightL + weightR)
	
	local Cs = sqrt((gamma - 1) * (hTotal - .5 * vx^2))
	
	S[1] = vx - Cs
	S[2] = vx
	S[3] = vx + Cs
	
	return S
end

-- used by Roe
-- TODO fluxMatrix and reconstruction error is broken
function Euler1D:calcInterfaceEigenBasis(sim,i,qL,qR)
	local gamma = self.gamma
	
	local rhoL = qL[1]
	local vxL = qL[2] / rhoL 
	local eTotalL = qL[3] / rhoL
	local eIntL = eTotalL - .5 * vxL^2
	local PL = (gamma - 1) * rhoL * eIntL
	local hTotalL = eTotalL + PL / rhoL
	local weightL = sqrt(rhoL)
	
	local rhoR = qR[1]
	local vxR = qR[2] / rhoR 
	local eTotalR = qR[3] / rhoR
	local eIntR = eTotalR - .5 * vxR^2
	local PR = (gamma - 1) * rhoR * eIntR
	local hTotalR = eTotalR + PR / rhoR
	local weightR = sqrt(rhoR)
	
	local rho = sqrt(weightL * weightR)
	local vx = (weightL * vxL + weightR * vxR) / (weightL + weightR)
	local hTotal = (weightL * hTotalL + weightR * hTotalR) / (weightL + weightR)
	
	local Cs = sqrt((gamma - 1) * (hTotal - .5 * vx^2))
	
	local F = sim.fluxMatrix[i]
	F[1][1] = 0
	F[1][2] = 1
	F[1][3] = 0
	F[2][1] = .5 * (gamma - 3) * vx * vx
	F[2][2] = (3 - gamma) * vx
	F[2][3] = gamma - 1
	F[3][1] = vx * (.5 * (gamma - 1) * vx * vx - hTotal)
	F[3][2] = hTotal - (gamma - 1) * vx * vx
	F[3][3] = gamma * vx
	
	local S = sim.eigenvalues[i]
	S[1] = vx - Cs
	S[2] = vx
	S[3] = vx + Cs
	
	local U = sim.eigenvectors[i]
	-- slow
	U[1][1] = 1
	U[2][1] = vx - Cs
	U[3][1] = hTotal - Cs * vx
	-- vel
	U[1][2] = 1
	U[2][2] = vx
	U[3][2] = .5 * vx * vx
	-- fast
	U[1][3] = 1
	U[2][3] = vx + Cs
	U[3][3] = hTotal + Cs * vx
	
	-- [[ symbolically
	local V = sim.eigenvectorsInverse[i]
	V[1][1] = (.5 * (gamma - 1) * vx^2 + Cs * vx) / (2 * Cs^2)
	V[1][2] = -(Cs + (gamma - 1) * vx) / (2 * Cs^2)
	V[1][3] = (gamma - 1) / (2 * Cs^2)
	V[2][1] = 1 - (gamma - 1) * vx^2 / (2 * Cs^2)
	V[2][2] = (gamma - 1) * vx / Cs^2
	V[2][3] = -(gamma - 1) / Cs^2
	V[3][1] = (.5 * (gamma - 1) * vx^2 - Cs * vx) / (2 * Cs^2)
	V[3][2] = (Cs - (gamma - 1) * vx) / (2 * Cs^2)
	V[3][3] = (gamma - 1) / (2 * Cs^2)
	--]]
	--[[ numerically via cramers rule
	local mat33 = require 'mat33'
	sim.eigenvectorsInverse[i] = mat33.inv(U)
	--]]
end

--[[ do something to prove source terms are working ...
function Euler1D:sourceTerm(sim, qs)
	local source = sim:newState()
	for i=1,sim.gridsize do
		local x = sim.xs[i]
		local rho = qs[i][1]
		local v = math.abs(x) < .1 and -.01 or 0
		qs[i][2] = qs[i][2] + rho * v
		qs[i][3] = qs[i][3] + .5 * rho * v * v
	end
	return source
end
--]]

return Euler1D
