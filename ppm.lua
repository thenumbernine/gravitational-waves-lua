local class = require 'ext.class'
local math = require 'ext.math'
local matrix = require 'matrix'

local function PPMBehavior(parentClass)

	local PPMTemplate = class(parentClass)

	function PPMTemplate:init(args)
		self.equation = assert(args.equation or self.equation)
	
		-- disable flux limiter
		args = table(args)
		local limiter = require 'limiter' 
		args.fluxLimiter = limiter.donorCell
		
		PPMTemplate.super.init(self, args)
		self.name = self.name .. ' PPM'
		
		self.qLs = matrix.zeros{self.gridsize, self.numStates}
		self.qRs = matrix.zeros{self.gridsize, self.numStates}
	end

	function PPMTemplate:reset()
		PPMTemplate.super.reset(self)

		-- centered
		for i=1,self.gridsize do
			for j=1,self.numStates do
				self.qLs[i][j] = 0
				self.qRs[i][j] = 0
			end
		end
	end
	
	function PPMTemplate:calcFluxes(dt)
		local gridsize = self.gridsize
		local numStates = self.numStates
		local eqn = self.equation

		local applyPPMToCons = true

		-- calc prims
		local W
		if not applyPPMToCons then
			W = matrix.zeros{gridsize, numStates}
			for i=1,gridsize do
				W[i] = matrix{eqn:calcPrimFromCons(table.unpack(self.qs[i]))}
			end
		else
			W = self.qs
		end

		-- calc cell-centered eigenbasis wrt prims
		local lambdas = matrix.zeros{gridsize, numStates}
		local evls = matrix.zeros{gridsize, numStates, numStates}
		local evrs = matrix.zeros{gridsize, numStates, numStates}
		for i=1,gridsize do
			if not applyPPMToCons then
				--[[
				eqn:calcEigenBasisWrtPrimFromPrim(lambdas[i], evrs[i], evls[i], W[i])
				--]]
				-- [[
				local gamma = self.equation.gamma
				local rho, v, ETotal = table.unpack(W[i])
				local dUdW = matrix{
					{1, 0, 0},
					{rho, v, 0},
					{.5 * v*v, rho * v, 1/(gamma-1)},
				}
				local dWdU = matrix{
					{1, 0, 0},
					{-v/rho, 1/rho, 0},
					{(gamma-1)*.5*v*v, -(gamma-1)*v, (gamma-1)},
				}

				eqn:calcEigenBasisFromCons(lambdas[i], evrs[i], evls[i], nil, self.qs[i])
				evls[i] = evls[i] * dUdW
				evrs[i] = dWdU * evrs[i]
				--]]	
			else
				eqn:calcEigenBasisFromCons(lambdas[i], evrs[i], evls[i], nil, self.qs[i])
			end
		end

		local dWm = matrix.zeros{gridsize, numStates}
		for i=2,gridsize-1 do
			-- calc prim differences across cell
			local dWc = matrix.zeros{numStates}
			local dWl = matrix.zeros{numStates}
			local dWr = matrix.zeros{numStates}
			local dWg = matrix.zeros{numStates}
			for j=1,numStates do
				dWc[j] = W[i+1][j] - W[i-1][j]
				dWl[j] = W[i][j] - W[i-1][j]
				dWr[j] = W[i+1][j] - W[i][j]
				if dWl[j] * dWr[j] > 0 then
					dWg[j] = 2 * dWl[j] * dWr[j] / (dWl[j] + dWr[j])
				else
					dWg[j] = 0
				end
			end

			-- project prim differences into eigenbasis 
			local dac = matrix()
			local dal = matrix()
			local dar = matrix()
			local dag = matrix()
			for j=1,numStates do
				dac[j] = 0
				dal[j] = 0
				dar[j] = 0
				dag[j] = 0
				for k=1,numStates do
					dac[j] = dac[j] + evls[i][j][k] * dWc[k]
					dal[j] = dal[j] + evls[i][j][k] * dWl[k]
					dar[j] = dar[j] + evls[i][j][k] * dWr[k]
					dag[j] = dag[j] + evls[i][j][k] * dWg[k]
				end
			end
		
			local da = matrix()
			for j=1,numStates do
				da[j]  = 0
				if dal[j] * dar[j] > 0 then
					local lim_slope1 = math.min(math.abs(dal[j]), math.abs(dar[j]))
					local lim_slope2 = math.min(.5 * math.abs(dac[j]), math.abs(dag[j]))
					da[j] = math.sign(dac[j]) * math.min(2 * lim_slope1, lim_slope2)
				end
			end
			
			for j=1,numStates do
				dWm[i][j] = 0
				for k=1,numStates do
					dWm[i][j] = dWm[i][j] + evrs[i][j][k] * da[k]
				end
			end
		end
		
		-- interface primitives based on parabolic reconstruction
		local iW = matrix.zeros{gridsize, numStates}
		for i=2,gridsize do
			for j=1,numStates do
				iW[i][j] = .5 * (W[i][j] + W[i-1][j]) - (dWm[i][j] - dWm[i-1][j]) / 6
			end
		end
		
		for i=2,gridsize-1 do
			-- cell left and right values
			local Wlv = matrix(iW[i])
			local Wrv = matrix(iW[i+1])

			for j=1,numStates do
				local qa = (Wrv[j] - W[i][j]) * (W[i][j] - Wlv[j])
				local qb = Wrv[j] - Wlv[j]
				local qc = 6 * (W[i][j] - .5 * (Wlv[j] + Wrv[j]))
				if qa <= 0 then
					Wlv[j] = W[i][j]
					Wrv[j] = W[i][j]
				elseif qb * qc > qb * qb then
					Wlv[j] = 3 * W[i][j] - 2 * Wrv[j]
				elseif qb * qc < -qb * qb then
					Wrv[j] = 3 * W[i][j] - 2 * Wlv[j]
				end
			end
	
			-- [[
			for j=1,numStates do
				Wlv[j] = math.max(math.min(W[i][j], W[i-1][j]), Wlv[j])
				Wlv[j] = math.min(math.max(W[i][j], W[i-1][j]), Wlv[j])

				Wrv[j] = math.max(math.min(W[i][j], W[i+1][j]), Wrv[j])
				Wrv[j] = math.min(math.max(W[i][j], W[i+1][j]), Wrv[j])
			end
			--]]

			--[[
			local dW = matrix()
			local W6 = matrix()
			for j=1,numStates do
				dW[j] = Wrv[j] - Wlv[j]
				W6[j] = 6 * (W[i][j] - .5 * (Wlv[j] + Wrv[j]))
			end
			--]]
		
			if not applyPPMToCons then
				self.qLs[i+1] = matrix{eqn:calcConsFromPrim(table.unpack(Wrv))}
				self.qRs[i] = matrix{eqn:calcConsFromPrim(table.unpack(Wlv))}
			else
				self.qLs[i+1] = matrix(Wrv)
				self.qRs[i] = matrix(Wlv)
			end
		end
		
		-- now qLs and qRs can be used
		PPMTemplate.super.calcFluxes(self, dt)
	end

	function PPMTemplate:get_qL(i)
		assert(1 <= i and i <= self.gridsize)
		assert(self.qLs[i])
		for j=1,self.equation.numStates do
			assert(self.qLs[i][j])
		end
		return self.qLs[i]
	end

	function PPMTemplate:get_qR(i)
		assert(1 <= i and i <= self.gridsize)
		assert(self.qRs[i])
		for j=1,self.equation.numStates do
			assert(self.qRs[i][j], tostring(
				require 'ext.tolua'{q=self.qRs[i][j], i=i, j=j}
			))
		end
		return self.qRs[i]
	end

	return PPMTemplate
end

return PPMBehavior
