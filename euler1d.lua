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

Euler1D.numStates = 3
Euler1D.gamma = 5/3	

Euler1D.graphInfos = table{
	{viewport={0/3, 0/2, 1/3, 1/2}, getter=function(self,i) return self.qs[i][1] end, name='rho', color={1,0,1}},
	{viewport={1/3, 0/2, 1/3, 1/2}, getter=function(self,i) return self.qs[i][2] / self.qs[i][1] end, name='u', color={0,1,0}},
	{viewport={2/3, 0/2, 1/3, 1/2}, getter=function(self,i) return self.qs[i][3] / self.qs[i][1] end, name='E', color={.5,.5,1}},
	{viewport={0/3, 1/2, 1/3, 1/2}, getter=function(self,i) return self.eigenbasisErrors and math.log(self.eigenbasisErrors[i]) end, name='log eigenbasis error', color={1,0,0}, range={-30, 30}},
	{viewport={1/3, 1/2, 1/3, 1/2}, getter=function(self,i) return self.fluxMatrixErrors and math.log(self.fluxMatrixErrors[i]) end, name='log reconstruction error', color={1,0,0}, range={-30, 30}},
}
Euler1D.graphInfoForNames = Euler1D.graphInfos:map(function(info,i)
	return info, info.name
end)

function Euler1D:initCell(sim,i)
	--[[ constant
	local rho = 1
	local u = 0
	local E = 1 + .5 * u * u
	--]]
	--[[ linear
	local rho = 2 + sim.xs[i]
	local u = 0
	local E = 1
	--]]
	--[[ gaussian curve
	local x = sim.xs[i]
	local sigma = 1/math.sqrt(10)
	local rho = math.exp(-x^2/sigma^2) + .1
	-- drho/dx = (-2x/sigma^2) exp(-x^2/sigma^2)
	local u = 0
	local E = 1 + .5 * u * u + .1 * (math.exp(-x^2/sigma^2) + 1) / ((self.gamma - 1) * rho)
	--]]
	-- [[ Sod
	local rho = sim.xs[i] < 0 and 1 or .1
	local u = 0
	local E = 1	+ .5 * u * u	-- internal + kinetic
	--]]
	--[[ Sedov
	local rho = 1
	local u = 0
	local E
	if i == math.floor(sim.gridsize/2) then
		E = 1e+3 / ((self.gamma - 1) * rho)
	else
		E = 1 / ((self.gamma - 1) * rho)
	end
	--]]
	return {rho, rho * u, rho * E}
end

-- used by HLL
function Euler1D:calcFluxForState(q, flux)
	flux = flux or {}
	local gamma = self.gamma
	flux[1] = q[2]
	flux[2] = (gamma - 1) * q[3] + (3 - gamma) / 2 * q[2] * q[2] / q[1]
	flux[3] = gamma * q[2] * q[3] / q[1] + (1 - gamma) / 2 * q[2] * q[2] * q[2] / (q[1] * q[1])
	return flux
end

-- used by HLL
function Euler1D:calcInterfaceEigenvalues(sim, qL, qR, S)
	S = S or {}
	
	local gamma = self.gamma
	
	local rhoL = qL[1]
	local uL = qL[2] / rhoL 
	local EL = qL[3] / rhoL
	local eIntL = EL - .5 * uL^2
	local PL = (gamma - 1) * rhoL * eIntL
	local HL = EL + PL / rhoL
	local weightL = sqrt(rhoL)
	
	local rhoR = qR[1]
	local uR = qR[2] / rhoR 
	local ER = qR[3] / rhoR
	local eIntR = ER - .5 * uR^2
	local PR = (gamma - 1) * rhoR * eIntR
	local HR = ER + PR / rhoR
	local weightR = sqrt(rhoR)
	
	local u = (weightL * uL + weightR * uR) / (weightL + weightR)
	local H = (weightL * HL + weightR * HR) / (weightL + weightR)
	
	local Cs = sqrt((gamma - 1) * (H - .5 * u^2))
	
	S[1] = u - Cs
	S[2] = u
	S[3] = u + Cs
	
	return S
end

-- used by Roe
-- TODO fluxMatrix and reconstruction error is broken
function Euler1D:calcInterfaceEigenBasis(sim,i,qL,qR)
	local gamma = self.gamma
	
	local rhoL = qL[1]
	local uL = qL[2] / rhoL 
	local EL = qL[3] / rhoL
	local eIntL = EL - .5 * uL^2
	local PL = (gamma - 1) * rhoL * eIntL
	local HL = EL + PL / rhoL
	local weightL = sqrt(rhoL)
	
	local rhoR = qR[1]
	local uR = qR[2] / rhoR 
	local ER = qR[3] / rhoR
	local eIntR = ER - .5 * uR^2
	local PR = (gamma - 1) * rhoR * eIntR
	local HR = ER + PR / rhoR
	local weightR = sqrt(rhoR)
	
	local rho = sqrt(weightL * weightR)
	local u = (weightL * uL + weightR * uR) / (weightL + weightR)
	local H = (weightL * HL + weightR * HR) / (weightL + weightR)
	local E = (weightL * EL + weightR * ER) / (weightL + weightR)
	
	local Cs = sqrt((gamma - 1) * (H - .5 * u^2))
	
	local F = sim.fluxMatrix[i]
	F[1][1] = 0
	F[1][2] = 1
	F[1][3] = 0
	F[2][1] = (gamma-3)/2*u*u
	F[2][2] = (3-gamma)*u
	F[2][3] = gamma-1
	F[3][1] = -u*(gamma*E + (1-gamma)*u*u)
	F[3][2] = H + (1-gamma) * u*u
	F[3][3] = gamma * u
	
	local S = sim.eigenvalues[i]
	S[1] = u - Cs
	S[2] = u
	S[3] = u + Cs
	
	local U = sim.eigenvectors[i]
	U[1][1] = 1
	U[1][2] = 1
	U[1][3] = 1
	U[2][1] = u - Cs
	U[2][2] = u
	U[2][3] = u + Cs
	U[3][1] = H - Cs * u
	U[3][2] = .5 * u*u
	U[3][3] = H + Cs * u
	
	-- [[ symbolically
	local V = sim.eigenvectorsInverse[i]
	V[1][1] = (.5 * (gamma - 1) * u^2 + Cs * u) / (2 * Cs^2)
	V[1][2] = -(Cs + (gamma - 1) * u) / (2 * Cs^2)
	V[1][3] = (gamma - 1) / (2 * Cs^2)
	V[2][1] = 1 - (gamma - 1) * u^2 / (2 * Cs^2)
	V[2][2] = (gamma - 1) * u / Cs^2
	V[2][3] = -(gamma - 1) / Cs^2
	V[3][1] = (.5 * (gamma - 1) * u^2 - Cs * u) / (2 * Cs^2)
	V[3][2] = (Cs - (gamma - 1) * u) / (2 * Cs^2)
	V[3][3] = (gamma - 1) / (2 * Cs^2)
	--]]
	--[[ numerically via cramers rule
	local det = U[1][1] * U[2][2] * U[3][3]
			+ U[2][1] * U[3][2] * U[1][3]
			+ U[3][1] * U[1][2] * U[2][3]
			- U[3][1] * U[2][2] * U[1][3]
			- U[2][1] * U[1][2] * U[3][3]
			- U[1][1] * U[3][2] * U[2][3];
	if det == 0 then
		for j=1,3 do
			for k=1,3 do
				console.log('A('+i+','+j+') = '+U[j][k]);
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
			V[k][j] = invDet * (U[j1][k1] * U[j2][k2] - U[j1][k2] * U[j2][k1])
		end
	end
	--]]
end

return Euler1D

