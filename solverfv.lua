--[[
finite volume solver
has flux
--]]
local class = require 'ext.class'
local Solver = require 'solver' 

local SolverFV = class(Solver)

function SolverFV:init(args)
	SolverFV.super.init(self, args)
	
	-- matches Burgers and HLL
	self.fluxes = {}

	-- matches HLL
	self.eigenvalues = {}
end

function SolverFV:reset()
	SolverFV.super.reset(self)

	-- state interfaces
	for i=1,self.gridsize+1 do
		self.fluxes[i] = {}
		for j=1,self.numStates do
			self.fluxes[i][j] = 0
		end
	end
	
	for i=1,self.gridsize+1 do
		self.eigenvalues[i] = {}
		for j=1,self.numWaves do
			self.eigenvalues[i][j] = 0
		end
	end
end

function SolverFV:calcDT()
	-- [[ using interface
	for i=2,self.gridsize do
		local qL = self:get_qL(i)
		local qR = self:get_qR(i)
		self.equation:calcInterfaceEigenvalues(self, i, qL, qR)
	end
	--]]

	if self.fixed_dt then
		return self.fixed_dt
	else
		local result = math.huge
		for i=1,self.gridsize do
			-- [[ using interface 
			local eigenvaluesL = self.eigenvalues[i]
			local eigenvaluesR = self.eigenvalues[i+1]
			local lambdaMax = eigenvaluesL[#eigenvaluesL] -- math.max(0, unpack(eigenvaluesL))
			local lambdaMin = eigenvaluesR[1] -- math.min(0, unpack(eigenvaluesR))
			--]]
			--[[ using cell
			local lambdaMin, lambdaMax = self.equation:calcCellMinMaxEigenvalues(self, i)
			lambdaMin = math.min(0, lambdaMin)
			lambdaMax = math.max(0, lambdaMax)
			--]]
			local dx = self.ixs[i+1] - self.ixs[i]
			local dum = dx / (math.abs(lambdaMax - lambdaMin) + 1e-9)
			result = math.min(result, dum)
		end
		return result * self.cfl
	end
end

function SolverFV:calcDerivFromFluxes(dt)
	self:calcFluxes(dt)
	
	local dq_dts = self:newState()
	for i=1,self.gridsize do
		local dx = self.ixs[i+1] - self.ixs[i]
		for j=1,self.numStates do
			dq_dts[i][j] = dq_dts[i][j] - (self.fluxes[i+1][j] - self.fluxes[i][j]) / dx
		end
	end
	return dq_dts
end

return SolverFV
