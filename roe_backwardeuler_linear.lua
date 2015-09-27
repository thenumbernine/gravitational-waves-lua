local class = require 'ext.class'
local Roe = require 'roe'
local RoeBackwardEulerLinear = class(Roe)

--RoeBackwardEulerLinear.fixed_dt = 1/200

function RoeBackwardEulerLinear:init(args)
	RoeBackwardEulerLinear.super.init(self, args)
	self.linearSolver = args.linearSolver or require 'linearsolvers'.conjres
	self.errorLogging = args.errorLogging
end

function RoeBackwardEulerLinear:integrateFlux(dt)
	--calcDT extra params: getLeft, getRight, getLeft2, getRight2
	-- the previous call in Solver:iterate is to self:calcDT
	-- which calls self.equation:calcInterfaceEigenBasis, which fills the self.eigenvalues table

	self:calcDeltaQTildes()
	self:calcRTildes()
	self:calcPhis(dt)

	-- function that returns deriv when provided a state vector
	-- TODO make this consider getLeft/getRight ... which themselves are not modular wrt state vector
	local function dq_dt(qs)
		local getLeft = function(i) return qs[i-1] end
		local getRight = function(i) return qs[i] end
		return self:calcDeriv(getLeft, getRight)
	end

	-- [[ implicit via some linear solver
	local qs = self.qs
	self.qs = self.linearSolver{
		--maxiter = 1000,
		x0 = qs:clone(),
		-- [=[ backward Euler
		epsilon = 1e-20,
		maxiter = 200,
		b = qs:clone(),
		A = function(qs)
			qs = qs - dt * dq_dt(qs)
			self.boundaryMethod(qs)
			return qs
		end,
		--]=]
		--[=[ crank-nicolson - converges faster
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
		-- not true for any 2nd derivative terms
		-- this method is only used for Jacobi method, so I don't really care
		ADiag = (function()
			local n = self:newState()
			for i=1,self.gridsize do
				for j=1,self.numStates do
					n[i][j] = 1
				end
			end
			return n
		end)(),
		-- logging:
		errorCallback = self.errorLogging and function(err, convergenceIteration)
			print(self.t, convergenceIteration, err)
		end,
	}
	if self.errorLogging then
		print()
	end
	--]]
	--[[ explicit - forward Euler - for debugging
	self.qs = self.qs + dt * dq_dt(self.qs)
	--]]
end

return RoeBackwardEulerLinear
