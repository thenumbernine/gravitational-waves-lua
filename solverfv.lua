--[[
finite volume solver
has flux
--]]
local class = require 'ext.class'
local Solver = require 'solver' 
local matrix = require 'matrix'

local SolverFV = class(Solver)

function SolverFV:init(args)
	SolverFV.super.init(self, args)
	
	self.fluxLimiter = assert(args.fluxLimiter)
	
	if self.fluxLimiter.name ~= 'donorCell' then
		self.name = self.name .. ' flux lim.=' .. self.fluxLimiter.name
	end
	
	self.fluxes = matrix()
end

function SolverFV:reset()
	SolverFV.super.reset(self)

	-- state interfaces
	for i=1,self.gridsize+1 do
		self.fluxes[i] = matrix()
		for j=1,self.numStates do
			self.fluxes[i][j] = 0
		end
	end
end

function SolverFV:calcDT()
	if self.fixed_dt then
		return self.fixed_dt
	else
		local result = math.huge
		for i=self.numGhost+1,self.gridsize-self.numGhost do
			--[[ using interface 
			local eigenvaluesL = self.eigenvalues[i]
			local eigenvaluesR = self.eigenvalues[i+1]
			local lambdaMax = eigenvaluesL[#eigenvaluesL] -- math.max(0, unpack(eigenvaluesL))
			local lambdaMin = eigenvaluesR[1] -- math.min(0, unpack(eigenvaluesR))
			--]]
			-- [[ using cell
			local lambdaMin, lambdaMax = self.equation:calcCellMinMaxEigenvalues(self, i)
			local lambdaAbsMax = math.max(
				math.abs(lambdaMin),
				math.abs(lambdaMax),
				1e-9)
			--]]
			local dx = self.ixs[i+1] - self.ixs[i]
			local dum = dx / lambdaAbsMax
			result = math.min(result, dum)
		end
		return result * self.cfl
	end
end

function SolverFV:calcDerivFromFluxes(dt)	
self:applyBoundary()

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
