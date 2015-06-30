--[[
q_i(t+dt) = q_i(t) + dq_i/dt(t+dt)
q_i(t+dt) - dt dq_i/dt(t+dt) = q_i(t)
(delta_ij - dt A_ij) q_j(t+dt) = q_i(t)
...for A_ij = linear components of dq_i/dt(t+dt) wrt q_j(t+dt)
q_i(t+dt) = ||(delta_ij - dt A_ij)^-1||_ij q_j(t)
--]]

local class = require 'ext.class'
local table = require 'ext.table'

local conjres = require 'conjres'
local Solver = require 'solver'

local EulerBackwardEulerConjRes = class(Solver)

EulerBackwardEulerConjRes.fixed_dt = 1/10

function EulerBackwardEulerConjRes:init(args)
	args = table(args)
	args.equation = require 'euler1d'(args)
	EulerBackwardEulerConjRes.super.init(self, args)
end

function EulerBackwardEulerConjRes:iterate()
	local dt = self.fixed_dt
	local gamma = self.equation.gamma
	local q = self.qs
	self.qs = conjres{
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

return EulerBackwardEulerConjRes

