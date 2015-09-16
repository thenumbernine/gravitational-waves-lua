local class = require 'ext.class'
local Roe = require 'roe'
local RoeBackwardEulerConjRes = class(Roe)

RoeBackwardEulerConjRes.fixed_dt = 1/50

function RoeBackwardEulerConjRes:init(args)
	RoeBackwardEulerConjRes.super.init(self, args)

	self.linearSolver = args.linearSolver or require 'linearsolvers'.conjres

	self.Phis = {}
end

function RoeBackwardEulerConjRes:reset()
	RoeBackwardEulerConjRes.super.reset(self)

	for i=1,self.gridsize do
		self.Phis[i] = {}
		for j=1,self.numStates do
			self.Phis[i][j] = 1
		end
	end
end

function RoeBackwardEulerConjRes:integrateFlux(dt)
	-- the previous call in Solver:iterate is to self:calcDT
	-- which calls self.equation:calcInterfaceEigenBasis, which fills the self.eigenvalues table

	self:calcDeltaQTildes()
	self:calcRTildes()

	-- compute Phi vector along interfaces
	-- Phi = diag(sign(v_i) + 1/2 phi * (dt/dx v_i - sign(v_i)))
	-- TODO shares in common with Roe
	for i=2,self.gridsize do
		for j=1,self.numStates do
			local phi = self.fluxLimiter(self.rTildes[i][j])
			local theta = self.eigenvalues[i][j] >= 0 and 1 or -1
			local dx = self.xs[i] - self.xs[i-1]
			local epsilon = self.eigenvalues[i][j] * dt / dx
			self.Phis[i][j] = theta + .5 * phi * (epsilon - theta)
		end
	end
	
	local function applyLeft(i,x)
		x = self:eigenfields(i,x)
		for j=1,self.numStates do
			x[j] = .5 * x[j] * self.eigenvalues[i][j] * (1 + self.Phis[i][j])
		end
		return self:eigenfieldsInverse(i, x)
	end
	local function applyRight(i,x)
		x = self:eigenfields(i,x)
		for j=1,self.numStates do
			x[j] = .5 * x[j] * self.eigenvalues[i][j] * (1 - self.Phis[i][j])
		end
		return self:eigenfieldsInverse(i, x)
	end

	local q = self.qs
	self.qs = self.linearSolver{
		--maxiter = 1000,
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
		ADiag = (function()
			local n = self:newState()
			for i=1,self.gridsize do
				for j=1,self.numStates do
					n[i][j] = 1
				end
			end
			return n
		end)(),
	}
	self:boundaryMethod()
end

return RoeBackwardEulerConjRes

