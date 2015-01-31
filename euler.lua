require 'ext'
local Simulation = require 'simulation'
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

return EulerSim

