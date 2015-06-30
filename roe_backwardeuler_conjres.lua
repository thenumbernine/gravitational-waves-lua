local conjres = require 'conjres'

local RoeBackwardEulerConjRes = require 'ext.class'(require 'roe')

function RoeBackwardEulerConjRes:integrateFlux(dt)
	-- the previous call in Solver:iterate is to self:calcDT
	-- which calls self.equation:calcInterfaceEigenBasis, which fills the self.eigenvalues table

	-- TODO shares in common with Roe
	for i=2,self.gridsize do
		local qL = getLeft and getLeft(i) or self.qs[i-1]
		local qR = getRight and getRight(i) or self.qs[i]
		
		local dq = {}
		for j=1,self.numStates do
			dq[j] = qR[j] - qL[j]
		end
		self.deltaQTildes[i] = self:eigenfields(i, dq)
	end

	-- compute Phi vector along interfaces
	-- Phi = diag(sign(v_i) + 1/2 phi * (dt/dx v_i - sign(v_i)))
	-- TODO shares in common with Roe
	local Phi = {}
	for i=2,self.gridsize do
		Phi[i] = {}
		
		for j=1,self.numStates do
			local lambda = self.eigenvalues[i][j]
			local rTilde
			if self.deltaQTildes[i][j] == 0 then
				rTilde = 0
			else
				if lambda >= 0 then
					rTilde = self.deltaQTildes[i-1][j] / self.deltaQTildes[i][j]
				else
					rTilde = self.deltaQTildes[i+1][j] / self.deltaQTildes[i][j]
				end
			end
			local phi = self.fluxLimiter(rTilde)
			local theta = lambda >= 0 and 1 or -1
			local dx = self.xs[i] - self.xs[i-1]
			local epsilon = lambda * dt / dx
			Phi[i][j] = theta + .5 * phi * (epsilon - theta)
		end
	end
	
	local function applyLeft(i,x)
		x = self:eigenfields(i,x)
		for j=1,self.numStates do
			x[j] = .5 * x[j] * self.eigenvalues[i][j] * (1 + Phi[i][j])
		end
		return self:eigenfieldsInverse(i, x)
	end
	local function applyRight(i,x)
		x = self:eigenfields(i,x)
		for j=1,self.numStates do
			x[j] = .5 * x[j] * self.eigenvalues[i][j] * (1 - Phi[i][j])
		end
		return self:eigenfieldsInverse(i, x)
	end

	local q = self.qs
	self.qs = conjres{
		x0 = q:clone(),
		b = q:clone(),
		A = function(x)
			local y = self:newState()
			for i=2,self.gridsize-1 do
				-- ALeft_{i+1/2} = 1/2 Q V (I + Phi) Q^-1
				-- ARight_{i+1/2} = 1/2 Q V (I - Phi) Q^-1
				if i>1 then
					local left = applyLeft(i,x[i-1]) -- apply ALeft_{i-1/2} to q_{i-1}
					for j=1,self.numStates do
						y[i][j] = y[i][j] + left[j]
					end
				end
				local right = applyRight(i,x[i]) -- apply ARight_{i-1/2} to q_{i} 
				local left = applyLeft(i+1,x[i]) -- apply -ALeft_{i+1/2} to q_{i}
				for j=1,self.numStates do
					y[i][j] = y[i][j] + right[j]
					y[i][j] = y[i][j] - left[j]
				end
				if i<self.gridsize then
					local right = applyRight(i+1,x[i+1]) -- apply -ARight_{i+1/2} to q_{i+1}
					for j=1,self.numStates do
						y[i][j] = y[i][j] - right[j]
					end
				end
				local dx = self.ixs[i+1] - self.ixs[i]
				for j=1,self.numStates do
					y[i][j] = x[i][j] - dt/dx * y[i][j]
				end
			end
			return y
		end,
	}
end

return RoeBackwardEulerConjRes

