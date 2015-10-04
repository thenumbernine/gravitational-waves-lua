local class = require 'ext.class'
local Roe = require 'roe'
local RoeImplicitLinearized = class(Roe)

--RoeImplicitLinearized.fixed_dt = 1/200

function RoeImplicitLinearized:init(args)
	RoeImplicitLinearized.super.init(self, args)
	self.linearSolver = args.linearSolver or require 'linearsolvers'.conjres
	self.errorLogging = args.errorLogging
end

function RoeImplicitLinearized:integrateFlux(dt)
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
		-- this is only the flux deriv... right?  or does it include the source term as well?
		return self:calcDeriv(getLeft, getRight)		
	end

-- [[ implicit via some linear solver
	local qs = self.qs

	local linearSolverArgs = {
		--maxiter = 1000,
		x0 = qs:clone(),
		epsilon = 1e-50,
		maxiter = 10000,
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
		-- if this is a linearized implicit solver
		-- then the matrix should be computed before invoking the iterative solver
		-- which means the matrix coeffiicents shouldn't be changing per-iteration
		-- which means dq_dt() should be based on the initial state and not the iterative state
		--qs = qs - dt * dq_dt(self.qs)
		-- ... but 
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

return RoeImplicitLinearized
