local class = require 'ext.class'
local SolverFV = require 'solverfv'

local HLL = class(SolverFV)

function HLL:init(args)
	self.equation = assert(self.equation or args.equation)
	HLL.super.init(self, args)
end

function HLL:calcDT(getLeft, getRight)
	-- matches Roe, except without eigenvectors
	for i=2,self.gridsize do
		local qL = getLeft and getLeft(self,i) or self.qs[i-1]
		local qR = getRight and getRight(self,i) or self.qs[i]
		self.equation:calcInterfaceEigenvalues(self, i, qL, qR, self.eigenvalues[i])
	end

	return HLL.super.calcDT(self)
end
	
function HLL:calcFlux(dt, getLeft, getRight, getLeft2, getRight2)
	local gamma = self.gamma
	
	local iqs = self:newState()
	
	for i=2,self.gridsize do
		-- TODO use qL and qR to allow compatability with MUSCL
		local qL = getLeft and getLeft(self,i) or self.qs[i-1]
		local qR = getRight and getRight(self,i) or self.qs[i]
		
		local sL = self.eigenvalues[i][1]
		local sR = self.eigenvalues[i][self.numStates]

		local fluxL = self.equation:calcFluxForState(self, i, qL)
		local fluxR = self.equation:calcFluxForState(self, i, qR)

		local flux = self.fluxes[i]
		for i=1,self.numStates do
			if 0 <= sL then
				flux[i] = fluxL[i]
			elseif sL <= 0 and 0 <= sR then
				flux[i] = (sR * fluxL[i] - sL * fluxR[i] + sL * sR * (qR[i] - qL[i])) / (sR - sL)
			elseif sR <= 0 then
				flux[i] = fluxR[i]
			end
		end
	end
	
	local dq_dts = self:newState()
	for i=1,self.gridsize do
		local dx = self.ixs[i+1] - self.ixs[i]
		for j=1,self.numStates do
			dq_dts[i][j] = dq_dts[i][j] - (self.fluxes[i+1][j] - self.fluxes[i][j]) / dx
		end
	end

	return dq_dts
end

return HLL

