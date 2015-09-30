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
		-- this is only the flux deriv...
		-- TODO rename this to "calcFluxDeriv"
		local dq_dt = self:calcDeriv(getLeft, getRight)		
		-- so we gotta add source terms as well ...
		if self.equation.sourceTerm	then
			dq_dt = dq_dt + self.equation:sourceTerm(self, qs)
		end
		if self.equation.postIterate then
			self.equation:postIterate(self, qs)
		end
		return dq_dt
	end

	-- [[ implicit via some linear solver
	local qs = self.qs

	local linearSolverArgs = {
		--maxiter = 1000,
		x0 = qs:clone(),
		epsilon = 1e-20,
		maxiter = 100,
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
	
	-- [=[ backward Euler
	linearSolverArgs.b = qs:clone()
	linearSolverArgs.A = function(qs)
		qs = qs - dt * dq_dt(qs)
		self.boundaryMethod(qs)
		return qs
	end
	--]=]
	--[=[ crank-nicolson - converges faster
	linearSolverArgs.b = (function(qs)
		qs = qs + .5 * dt * dq_dt(qs)
		self.boundaryMethod(qs)
		return qs
	end)(qs)
	linearSolverArgs.A = function(qs)
		qs = qs - .5 * dt * dq_dt(qs)
		self.boundaryMethod(qs)
		return qs
	end
	--]=]

	self.qs = self.linearSolver(linearSolverArgs)
	if self.errorLogging then
		print()
	end
	--]]
	--[[ explicit - forward Euler - for debugging
	self.qs = self.qs + dt * dq_dt(self.qs)
	--]]
end

return RoeBackwardEulerLinear
