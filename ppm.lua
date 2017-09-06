local class = require 'ext.class'
local math = require 'ext.math'
local matrix = require 'matrix'

local function PPMBehavior(parentClass)

	local PPMTemplate = class(parentClass)

	function PPMTemplate:init(args)
		self.equation = assert(args.equation or self.equation)
		self.name = self.equation.name .. ' PPMTemplate'
	
		-- disable flux limiter
		args = table(args)
		local limiter = require 'limiter' 
		args.fluxLimiter = limiter.donorCell
		
		PPMTemplate.super.init(self, args)
		
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

		--[[ calc prims
		local W = matrix.zeros{gridsize, numStates}
		for i=1,gridsize do
			W[i] = matrix{eqn:calcPrimFromCons(table.unpack(self.qs[i]))}
		end
		--]]
		-- [[
		local W = self.qs
		--]]

		-- calc cell-centered eigenbasis wrt prims
		local lambdas = matrix.zeros{gridsize, numStates}
		local evls = matrix.zeros{gridsize, numStates, numStates}
		local evrs = matrix.zeros{gridsize, numStates, numStates}
		for i=1,gridsize do
			--eqn:calcEigenBasisWrtPrimFromPrim(lambdas[i], evrs[i], evls[i], W[i])
			eqn:calcEigenBasisFromCons(lambdas[i], evrs[i], evls[i], nil, self.qs[i])
		end

		local dWm = matrix.zeros{gridsize, numStates}
		local Wim1h = matrix.zeros{gridsize, numStates}
		for i=3,gridsize-1 do
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
			local dac = matrix.zeros{numStates}
			local dal = matrix.zeros{numStates}
			local dar = matrix.zeros{numStates}
			local dag = matrix.zeros{numStates}
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
		
			local da = matrix.zeros{numStates}
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
		
			for j=1,numStates do
				Wim1h[i][j] = .5 * (W[i][j] + W[i-1][j]) - (dWm[i][j] - dWm[i-1][j]) / 6
			end

			local Wlv = matrix(Wim1h[i-1])	-- notice this depends on loop order 
			local Wrv = matrix(Wim1h[i])

			local gamma_curv = 0
			for j=1,numStates do
				local qa = (Wrv[j] - W[i-1][j]) * (W[i-1][j] - Wlv[j])
				local qb = Wrv[j] - Wlv[j]
				local qc = 6 * (W[i-1][j] - .5 * (Wlv[j] * (1 - gamma_curv) + Wrv[j] * (1 + gamma_curv)))
				if qa < 0 then
					Wlv[j] = W[i-1][j]
					Wrv[j] = W[i-1][j]
				elseif qb * qc > qb * qb then	-- why not divide by qb?
					Wlv[j] = (6 * W[i-1][j] - Wrv[j] * (4 + 3 * gamma_curv)) / (2 - 3 * gamma_curv)
				elseif qb * qc < -qb * qb then
					Wrv[j] = (6 * W[i-1][j] - Wlv[j] * (4 - 3 * gamma_curv)) / (2 + 3 * gamma_curv)
				end
			end
		
			for j=1,numStates do
				Wlv[j] = math.max( math.min( W[i-1][j], W[i-2][j]), Wlv[j])
				Wlv[j] = math.min( math.max( W[i-1][j], W[i-2][j]), Wlv[j])

				Wrv[j] = math.max( math.min( W[i-1][j], W[i][j]), Wrv[j])
				Wrv[j] = math.min( math.max( W[i-1][j], W[i][j]), Wrv[j])
			end
		
			local dW = matrix.zeros{numStates}
			local W6 = matrix.zeros{numStates}
			for j=1,numStates do
				dW[j] = Wrv[j] - Wlv[j]
				W6[j] = 6 * (W[i-1][j] - .5 * (Wlv[j] * (1 - gamma_curv) + Wrv[j] * (1 + gamma_curv)))
			end
		
			-- now we have Wl = Wrv and Wr = Wlv
			--[[
			self.qLs[i+1] = matrix{eqn:calcConsFromPrim(table.unpack(Wrv))}
			self.qRs[i] = matrix{eqn:calcConsFromPrim(table.unpack(Wlv))}
			--]]
			-- [[
			self.qLs[i+1] = matrix(Wrv)
			self.qRs[i] = matrix(Wlv)
			--]]
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
