local table = require 'ext.table'
local class = require 'ext.class'
local complex = require 'symmath.complex'
local matrix = require 'matrix'

local Equation = require 'equation'
local NLSEqn = class(Equation)
NLSEqn.numStates = 2

do
	local q = function(self, i) return self.qs[i] end
	NLSEqn:buildGraphInfos{
		{re = q:_(1)},
		{im = q:_(2)},
		{norm = math.sqrt:o(q:_(1)^2 + q:_(2)^2)},
	}
end

function NLSEqn:initCell(sim, i)
	local r = sim.xs[i]
	local A = 10
	q = A * complex.exp(-r^2)
	return {q:unpack()}
end


local Solver = require 'solver'
local NonLinearSchrodinger = class(Solver)
NonLinearSchrodinger.name = 'NonLinearSchrodinger'

NonLinearSchrodinger.fixed_dt = .000001

function NonLinearSchrodinger:init(...)
	NonLinearSchrodinger.equation = NLSEqn(self)
	NonLinearSchrodinger.super.init(self, ...)
end

local i = complex(0,1)

function NonLinearSchrodinger:calcDeriv(dt)
	local qs = self.qs
	local dq_dts = self:newState()
	for j=3,self.gridsize-3 do
		local r = self.xs[j]
		local dr = self.ixs[j+1] - self.ixs[j]
		local q_2L = complex(unpack(qs[j-2]))
		local q_1L = complex(unpack(qs[j-1]))
		local q = complex(unpack(qs[j]))
		local q_1R = complex(unpack(qs[j+1]))
		local q_2R = complex(unpack(qs[j+2]))
		local d2q_dx2 = (-q_2L + 16 * q_1L - 30 * q + 16 * q_1R - q_2R) / (12 * dr * dr)
		local dq_dx = (q_2L - 8 * q_1L + 8 * q_1R - q_2R) / (12 * dr)
		local q2 = q:norm()
		local q4 = q2 * q2
		dq_dts[j] = {(-i * (-d2q_dx2 - 4 / r * dq_dx + q4 * q)):unpack()} 
	end
	return dq_dts 
end

function NonLinearSchrodinger:step(dt)
	self:integrate(dt, function()
		return self:calcDeriv(dt)
	end)
end

function NonLinearSchrodinger:applyBoundary()
	-- freeflow on the left
	self.qs[1] = {unpack(self.qs[5])}
	self.qs[2] = {unpack(self.qs[4])}
	self.qs[#self.qs] = {0,0}
	self.qs[#self.qs-1] = {0,0}
end

return NonLinearSchrodinger 
