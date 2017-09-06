local class = require 'ext.class'
local SolverFV = require 'solverfv'

local HLL = class(SolverFV)

function HLL:init(args)
	self.equation = assert(self.equation or args.equation)
	self.name = self.equation.name .. ' HLL'

	self.useDirect = args.useDirect

	HLL.super.init(self, args)
end

function HLL:calcFluxes(dt)
	for i=2,self.gridsize do
		local qL = self:get_qL(i)
		local qR = self:get_qR(i)

		local lambdaInt = {} 
		fill(lambdaInt, self.equation:calcEigenvalues(self.equation:calcRoeValues(qL, qR)))

		local lambdaL = {}
		fill(lambdaL, self.equation:calcEigenvaluesFromCons(table.unpack(qL)))
		
		local lambdaR = {}
		fill(lambdaR, self.equation:calcEigenvaluesFromCons(table.unpack(qR)))

		local sL, sR
		if self.useDirect then
			sL = lambdaL[1]
			sR = lambdaR[self.numWaves]
		else
			sL = math.min(lambdaInt[1], lambdaL[1])
			sR = math.max(lambdaInt[self.numWaves], lambdaR[self.numWaves])
		end

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
