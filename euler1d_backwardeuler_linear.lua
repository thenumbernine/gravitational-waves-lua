--[[
q_i(t+dt) = q_i(t) + dq_i/dt(t+dt)
q_i(t+dt) - dt dq_i/dt(t+dt) = q_i(t)
(delta_ij - dt A_ij) q_j(t+dt) = q_i(t)
...for A_ij = linear components of dq_i/dt(t+dt) wrt q_j(t+dt)
q_i(t+dt) = ||(delta_ij - dt A_ij)^-1||_ij q_j(t)
--]]

local class = require 'ext.class'
local table = require 'ext.table'

-- I'm only using SolverFV over Solver for its calcDT
local SolverFV = require 'solverfv'

local EulerBackwardEulerLinear = class(SolverFV)

EulerBackwardEulerLinear.equation = require 'euler1d'()

function EulerBackwardEulerLinear:init(args)
	EulerBackwardEulerLinear.super.init(self, args)
	self.linearSolver = assert(args.linearSolver)
end

function EulerBackwardEulerLinear:iterate()
	self:applyBoundary()
	local dt = self:calcDT()
	local gamma = self.equation.gamma
	local q = self.qs
	self.qs = self.linearSolver{
		b = q:clone(),
		x0 = q:clone(),
		A = function(x)
			local y = x:clone()
			for i=1,self.gridsize do
				local dx = self.ixs[i+1] - self.ixs[i]
				if i>1 then
					y[i][1] = y[i][1] - dt * -( -x[i-1][2] ) / (2*dx)
					y[i][2] = y[i][2] - dt * -( (1-gamma) * x[i-1][3] + (gamma-3)/2*x[i-1][2]*q[i-1][2]/q[i-1][1] ) / (2*dx)
					y[i][3] = y[i][3] - dt * -( .5*(gamma-1) * x[i-1][2]*q[i-1][2]^2/q[i-1][1]^2 - gamma * x[i-1][3]*q[i-1][2]/q[i-1][1] ) / (2*dx)
				end
				if i<self.gridsize then
					y[i][1] = y[i][1] - dt * ( -x[i+1][2] ) / (2*dx)
					y[i][2] = y[i][2] - dt * ( (1 - gamma) * x[i+1][3] + (gamma-3)/2*x[i+1][2]*q[i+1][2]/q[i+1][1] )/(2*dx)
					y[i][3] = y[i][3] - dt * ( .5*(gamma-1) * x[i+1][2]*q[i+1][2]^2/q[i+1][1]^2 - gamma * x[i+1][3]*q[i+1][2]/q[i+1][1] ) / (2*dx)
				end
			end
			return y
		end,
	}
end

return EulerBackwardEulerLinear
