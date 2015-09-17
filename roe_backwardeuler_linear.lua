local class = require 'ext.class'
local Roe = require 'roe'
local RoeBackwardEulerConjRes = class(Roe)

--RoeBackwardEulerConjRes.fixed_dt = 1/200

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

function RoeBackwardEulerConjRes:integrateFlux(dt)	--calcDT extra params: getLeft, getRight, getLeft2, getRight2
	-- the previous call in Solver:iterate is to self:calcDT
	-- which calls self.equation:calcInterfaceEigenBasis, which fills the self.eigenvalues table

	self:calcDeltaQTildes()
	self:calcRTildes()

	-- compute Phi vector along interfaces
	-- Phi = diag(sign(v_i) + 1/2 phi * (dt/dx v_i - sign(v_i)))
	-- TODO shares in common with Roe
	for i=2,self.gridsize do
		local dx = self.xs[i] - self.xs[i-1]
		for j=1,self.numStates do
			local phi = self.fluxLimiter(self.rTildes[i][j])
			local theta = self.eigenvalues[i][j] >= 0 and 1 or -1
			local epsilon = self.eigenvalues[i][j] * dt / dx
			self.Phis[i][j] = theta + phi * (epsilon - theta)
		end
	end

	local zero = {}
	for j=1,self.numStates do
		zero[j] = 0
	end
	
	local function dq_dt(qs)
		-- flux = ALeft * qs[i-1] + ARight * qs[i]
		-- for matrixes ALeft & ARight
		local function apply(i, q, s)
			if i == 1 or i == self.gridsize+1 then return zero end
			local qTilde = self:eigenfields(i, q)
			local fluxTilde = {}
			for j=1,self.numStates do
				fluxTilde[j] = .5 * self.eigenvalues[i][j] * qTilde[j] * (1 + s * self.Phis[i][j])
			end
			return self:eigenfieldsInverse(i, fluxTilde)
		end

		local y = self:newState()
		for i=1,self.gridsize do
			local dx = self.ixs[i+1] - self.ixs[i]
			
			local fluxLfromL, fluxLfromR = apply(i, qs[i-1], 1), apply(i, qs[i], -1)
			local fluxRfromL, fluxRfromR = apply(i+1, qs[i], 1), apply(i+1, qs[i+1], -1)
			
			for j=1,self.numStates do
				-- Forward-Euler integration:
				-- (Backward-Euler would be subtracting)
				y[i][j] = (fluxLfromL[j] + fluxLfromR[j] - fluxRfromL[j] - fluxRfromR[j]) / dx
			end
		end
		return y
	end

	-- [[ implicit via some linear solver
	-- boundaries are messing up ... 
	local qs = self.qs
	self.qs = self.linearSolver{
		--maxiter = 1000,
		x0 = qs:clone(),
		--[=[ backward Euler
		epsilon = 1e-50,
		b = qs:clone(),
		A = function(qs)
			qs = qs - dt * dq_dt(qs)
			self.boundaryMethod(qs)
			return qs
		end,
		--]=]
		-- [=[ crank-nicolson - converges faster
		epsilon = 1e-50,
		b = (function(qs)
			qs = qs + .5 * dt * dq_dt(qs)
			self.boundaryMethod(qs)
			return qs
		end)(qs),
		A = function(qs)
			qs = qs - .5 * dt * dq_dt(qs)
			self.boundaryMethod(qs)
			return qs
		end,
		--]=]
		-- mostly true ... mostly ...
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
	--]]
	--[[ explicit - should be equivalent to the explicit, at least
	-- but it's not ... so find out why ...
	self.qs = self.qs + dt * dq_dt(self.qs)
	--]]
end

return RoeBackwardEulerConjRes

