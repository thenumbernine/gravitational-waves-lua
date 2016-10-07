local class = require 'ext.class'
local Roe = require 'roe'

local RoeImplicitLinearized = class(Roe)

function RoeImplicitLinearized:init(args)
	RoeImplicitLinearized.super.init(self, args)
	
	self.linearSolver = args.linearSolver or require 'linearsolvers'.gmres
	self.linearSolverEpsilon = args.linearSolverEpsilon or 1e-10
	self.linearSolverMaxIter = args.linearSolverMaxIter or 10 * self.gridsize * self.numStates 
	self.linearSolverRestart = args.linearSolverRestart or self.gridsize * self.numStates 
	self.errorLogging = args.errorLogging

	self.name = self.equation.name .. ' Roe Implicit Linearized'
end

-- step contains integrating flux and source terms
-- but not post iterate
function RoeImplicitLinearized:step(dt)
	-- the previous call in Solver:iterate is to self:calcDT
	-- which calls self.equation:calcInterfaceEigenBasis, which fills the self.eigenvalues table

	-- function that returns deriv when provided a state vector
	local function calc_dq_dt(qs)
		local oldQs = rawget(self, 'qs')
		self.qs = qs

		-- this is only the flux deriv... right?  or does it include the source term as well?
		local dq_dt = self:calcDerivFromFluxes(dt)
		if self.equation.sourceTerm then
			dq_dt = dq_dt + self.equation:sourceTerm(self, qs)
		end
		
		rawset(self, 'qs', oldQs)
		return dq_dt
	end

-- [[ implicit via some linear solver
	local qs = self.qs

	local linearSolverArgs = {
		--maxiter = 1000,
		x0 = qs:clone(),
		epsilon = self.linearSolverEpsilon, 
		maxiter = self.linearSolverMaxIter,
		restart = self.linearSolverRestart,
		--[=[ mostly true ... mostly ...
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
		--]=]
		-- logging:
		errorCallback = self.errorLogging and function(err, convergenceIteration)
			print(self.t, convergenceIteration, err)
		end,
	}

	--[=[ identity.  do nothing.
	linearSolverArgs.b = qs:clone()
	linearSolverArgs.A = function(qs) return qs end
	--]=]
	-- [=[ backward Euler
	linearSolverArgs.b = qs:clone()
	linearSolverArgs.A = function(qs)
		-- if this is a linearized implicit solver
		-- then the matrix should be computed before invoking the iterative solver
		-- which means the matrix coeffiicents shouldn't be changing per-iteration
		-- which means calc_dq_dt() should be based on the initial state and not the iterative state
		qs = qs - dt * calc_dq_dt(self.qs)
		-- ... but 
		--qs = qs - dt * calc_dq_dt(qs)
		--self.boundaryMethod(qs)
		return qs
	end
	--]=]
	--[=[ crank-nicolson - converges faster
	linearSolverArgs.b = (function(qs)
		qs = qs + .5 * dt * calc_dq_dt(qs)
		self.boundaryMethod(qs)
		return qs
	end)(qs)
	linearSolverArgs.A = function(qs)
		qs = qs - .5 * dt * calc_dq_dt(qs)
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
	self.qs = self.qs + dt * calc_dq_dt(self.qs)
--]]
end

return RoeImplicitLinearized
