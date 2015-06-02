require 'ext'
local Simulation = require 'simulation'
local Euler1DSimulation = class(Simulation)

Euler1DSimulation.numStates = 3
Euler1DSimulation.gamma = 5/3	

function Euler1DSimulation:init(...)
	Simulation.init(self, ...)

	--index:bind(qs) => f(k) = qs[k] takes 1 arg and returns an array of 3 elements
	--index:bind(qs)[1] => f(k) = qs[k][1] takes 1 arg and returns the 1st element of the 3
	self.graphInfos = {
		{viewport={0/3, 0/2, 1/3, 1/2}, getter=function(i) return self.qs[i][1] end, name='rho', color={1,0,1}},
		{viewport={1/3, 0/2, 1/3, 1/2}, getter=function(i) return self.qs[i][2] / self.qs[i][1] end, name='u', color={0,1,0}},
		{viewport={2/3, 0/2, 1/3, 1/2}, getter=function(i) return self.qs[i][3] / self.qs[i][1] end, name='E', color={.5,.5,1}},
		{viewport={0/3, 1/2, 1/3, 1/2}, getter=log:compose(index:bind(self.eigenbasisErrors)), name='log eigenbasis error', color={1,0,0}, range={-30, 30}},
		{viewport={1/3, 1/2, 1/3, 1/2}, getter=log:compose(index:bind(self.fluxMatrixErrors)), name='log reconstruction error', color={1,0,0}, range={-30, 30}},
	}
end

function Euler1DSimulation:initCell(i)
	local rho = self.xs[i] < 0 and 1 or .1
	local u = 0
	local E = 1	+ .5 * u * u	-- internal + kinetic
	return {rho, rho * u, rho * E}
end

-- used by HLL
function Euler1DSimulation:calcFluxForState(q, flux)
	flux = flux or {}
	local gamma = self.gamma
	flux[1] = q[2]
	flux[2] = (gamma - 1) * q[3] + (3 - gamma) / 2 * q[2] * q[2] / q[1]
	flux[3] = gamma * q[2] * q[3] / q[1] + (1 - gamma) / 2 * q[2] * q[2] * q[2] / (q[1] * q[1])
	return flux
end

-- used by HLL
function Euler1DSimulation:calcInterfaceEigenvalues(qL, qR, S)
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
function Euler1DSimulation:calcInterfaceEigenBasis(i,qL,qR)
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
	
	local F = self.fluxMatrix[i]
	F[1][1] = 0
	F[1][2] = 1
	F[1][3] = 0
	F[2][1] = (gamma-3)/2*u*u
	F[2][2] = (3-gamma)*u
	F[2][3] = gamma-1
	F[3][1] = -u*(gamma*E + (1-gamma)*u*u)
	F[3][2] = H + (1-gamma) * u*u
	F[3][3] = gamma * u
	
	local S = self.eigenvalues[i]
	S[1] = u - Cs
	S[2] = u
	S[3] = u + Cs
	
	local U = self.eigenvectors[i]
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
	local V = self.eigenvectorsInverse[i]
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

return Euler1DSimulation

