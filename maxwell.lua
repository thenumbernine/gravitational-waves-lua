-- going by Trangenstein
require 'ext'
local Simulation = require 'simulation'
local MaxwellSimulation = class(Simulation)

MaxwellSimulation.numStates = 6	--E,B xyz

local e0 = 1	-- permittivity
local u0 = 1	-- permissivity

function MaxwellSimulation:init(...)
	MaxwellSimulation.super.init(self, ...)
	local Ex = index:bind(self.qs):index(1) / e0
	local Ey = index:bind(self.qs):index(2) / e0
	local Ez = index:bind(self.qs):index(3) / e0
	local ESq = Ex^2 + Ey^2 + Ez^2
	local Bx = index:bind(self.qs):index(4)
	local By = index:bind(self.qs):index(5)
	local Bz = index:bind(self.qs):index(6)
	local BSq = Bx^2 + By^2 + Bz^2
	self.graphInfos = {
		{viewport={0/4, 0/3, 1/4, 1/3}, getter=Ex, name='Ex', color={1,0,0}},
		{viewport={1/4, 0/3, 1/4, 1/3}, getter=Ey, name='Ey', color={0,1,0}},
		{viewport={2/4, 0/3, 1/4, 1/3}, getter=Ez, name='Ez', color={0,0,1}},
		{viewport={3/4, 0/3, 1/4, 1/3}, getter=sqrt:compose(ESq), name='E', color={0,0,1}},
		{viewport={0/4, 1/3, 1/4, 1/3}, getter=Bx, name='Bx', color={0,1,1}},
		{viewport={1/4, 1/3, 1/4, 1/3}, getter=By, name='By', color={1,0,1}},
		{viewport={2/4, 1/3, 1/4, 1/3}, getter=Bz, name='Bz', color={1,1,0}},
		{viewport={3/4, 1/3, 1/4, 1/3}, getter=sqrt:compose(BSq), name='B', color={1,1,0}},
		{viewport={0/4, 2/3, 1/4, 1/3}, getter=log:compose(index:bind(self.eigenbasisErrors)), name='log eigenbasis error', color={1,0,0}, range={-30, 30}},
		{viewport={1/4, 2/3, 1/4, 1/3}, getter=log:compose(index:bind(self.fluxMatrixErrors)), name='log reconstruction error', color={1,0,0}, range={-30, 30}},
		{viewport={3/4, 2/3, 1/4, 1/3}, getter=.5 * (ESq * e0 + BSq / u0), name='energy', color={0,.5,1}},
	}
end

function MaxwellSimulation:initCell(i)
	local x = self.xs[i]
	local Ex = 0
	local Ey = 0
	local Ez = 1
	local Bx = 1
	local By = x < 0 and 1 or -1
	local Bz = 0
	return {Ex * e0, Ey * e0, Ez * e0, Bx, By, Bz}
end

function MaxwellSimulation:calcInterfaceEigenBasis(i)
	local avgQ = {}
	for j=1,self.numStates do
		avgQ[j] = (self.qs[i-1][j] + self.qs[i][j]) / 2
	end
	local lambda = 1/sqrt(e0 * u0)
	self.eigenvalues[i] = {-lambda, -lambda, 0, 0, lambda, lambda}
	local se = sqrt(e0/2)
	local su = sqrt(u0/2)
	local seu = sqrt(e0/u0)/u0
	self.fluxMatrix[i] = {
		{0,0,0,0,0,0},
		{0,0,0,0,0,seu},
		{0,0,0,0,-seu,0},
		{0,0,0,0,0,0},
		{0,0,-1/seu,0,0,0},
		{0,1/seu,0,0,0,0}
	}
	self.eigenvectors[i] = {
		{0,0,-se,se,0,0},
		{0,-se,0,0,se,0},
		{se,0,0,0,0,-se},
		{0,0,su,su,0,0},
		{su,0,0,0,0,su},
		{0,su,0,0,su,0}
	}
	local ise = 1/se
	local isu = 1/su
	self.eigenvectorsInverse[i] = {
		{0,0,ise,0,isu,0},
		{0,-ise,0,0,0,isu},
		{-ise,0,0,isu,0,0},
		{ise,0,0,isu,0,0},
		{0,ise,0,0,0,isu},
		{0,0,-ise,0,isu,0}
	}
end

local sigma = 1	-- conductivity
function MaxwellSimulation:addSourceToDerivCell(dq_dts, i)
	local eEx, eEy, eEz, Bx, By, Bz = unpack(self.qs[i])
	local Ex = eEx / e0
	local Ey = eEy / e0
	local Ez = eEz / e0
	dq_dts[i][1] = dq_dts[i][1] - Ex * sigma
	dq_dts[i][2] = dq_dts[i][2] - Ey * sigma
	dq_dts[i][3] = dq_dts[i][3] - Ez * sigma
end

return MaxwellSimulation

