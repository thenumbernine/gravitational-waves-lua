local class = require 'ext.class'
local SolverFV = require 'solverfv'

local HLL = class(SolverFV)

function HLL:init(args)
	self.equation = assert(self.equation or args.equation)
	self.name = self.equation.name .. ' HLL'
	self.eigenvalues = {}
	
	HLL.super.init(self, args)
end

function HLL:reset()
	HLL.super.reset(self)

	for i=1,self.gridsize+1 do
		self.eigenvalues[i] = {}
		for j=1,self.numWaves do
			self.eigenvalues[i][j] = 0
		end
	end
end

function HLL:calcFluxes(dt)
	-- only the eigenvalues
	self:calcInterfaceEigenvalues()
	
	local gamma = self.gamma
	
	local iqs = self:newState()
	
	for i=2,self.gridsize do
		-- TODO use qL and qR to allow compatability with MUSCL
		local qL = self.qs[i-1]
		local qR = self.qs[i]
		
		local sL = self.eigenvalues[i][1]
		local sR = self.eigenvalues[i][self.numWaves]

		local fluxL = {self.equation:calcFluxForState(qL)}
		local fluxR = {self.equation:calcFluxForState(qR)}

		local flux = self.fluxes[i]
		for i=1,self.numWaves do
			if 0 <= sL then
				flux[i] = fluxL[i]
			elseif sL <= 0 and 0 <= sR then
				flux[i] = (sR * fluxL[i] - sL * fluxR[i] + sL * sR * (qR[i] - qL[i])) / (sR - sL)
			elseif sR <= 0 then
				flux[i] = fluxR[i]
			end
		end
	end
end

return HLL
