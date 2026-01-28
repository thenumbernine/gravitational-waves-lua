local table = require 'ext.table'
local class = require 'ext.class'
local complex = require 'complex'
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

--[[
the eqn is 
i u,t + u,jj - |u|^(p-1) u = 0
<=> u,t = i (u,jj - |u|^(p-1) u)
where p is always chosen to be odd

u,jj is evaluated in radial coordinates, which is
1/r^2 (r^2 u,r),r
= u,rr + 2/r u,r

substitute ...
u,t = i (u,rr + 2/r u,r - |u|^(p-1) u)
... so why is the discretized expression in the paper 4/r ?

... substitute finite difference coefficients:

u,rr = (-1/12 u[j+2] + 4/3 u[j+1] - 5/2 u[j] + 4/3 u[j-1] - 1/12 u[j-2]) / h^2
= (-u[j+2] + 16 u[j+1] - 30 u[j] + 16 u[j-1] - u[j-2]) / (12 h^2)

u,r = (-1/12 u[j+2] + 2/3 u[j+1] - 2/3 u[j-1] + 1/12 u[j-2]) / h
= (-u[j+2] + 8 u[j+1] - 8 u[j-1] + u[j-2]) / (12 h)

... still looks like the 4 in the paper should be a 2 

--]]

-- derivCoeffs[derivative][accuracy] = {coeffs...}
local derivCoeffs = {
	-- antisymmetric coefficients 
	{
		[2] = {.5},
		[4] = {2/3, -1/12},
		[6] = {3/4, -3/20, 1/60},
		[8] = {4/5, -1/5, 4/105, -1/280},
	},
	-- symmetric
	{
		[2] = {[0] = -2, 1},
		[4] = {[0] = -5/2, 4/3, -1/12},
		[6] = {[0] = -49/18, 3/2, -3/20, 1/90},
		[8] = {[0] = -205/72, 8/5, -1/5, 8/315, -1/560},
	},
	-- antisymmetric 
	{
		[2] = {-1, 1/2},
		[4] = {-13/8, 1, -1/8},
		[6] = {-61/30, 169/120, -3/10, 7/240},
	},
	-- symmetric
	{
		[2] = {[0] = 6, -4, 1},
		[4] = {[0] = 28/3, -13/2, 2, -1/6},
		[6] = {[0] = 91/8, -122/15, 169/60, -2/5, 7/240},
	},
	-- antisymmetric 
	{
		[2] = {5/2, -2, 1/2},
	},
	-- symmetric
	{
		[2] = {[0] = -20, 15, -6, 1},
	},
}

-- finite difference operator
local function D(qs, j, degree, order)
	local c = assert(derivCoeffs[degree][order], "don't have coefficients for degree "..degree.." order "..order)
	if degree % 2 == 0 then	-- symmetric
		local sum = complex(unpack(qs[j])) * c[0]
		for i=1,#c do
			sum = sum + (complex(unpack(qs[j-i])) + complex(unpack(qs[j+i]))) * c[i]
		end
		return sum
	else
		local sum = complex()
		for i=1,#c do
			sum = sum + (complex(unpack(qs[j+i])) - complex(unpack(qs[j-i]))) * c[i]
		end
		return sum
	end
end

function NonLinearSchrodinger:calcDeriv(dt)
	local qs = self.qs
	local dq_dts = self:newState()
	for j=3,self.gridsize-3 do
		local r = self.xs[j]
		local dr = self.ixs[j+1] - self.ixs[j]
		-- [[
		local q_2L = complex(unpack(qs[j-2]))
		local q_1L = complex(unpack(qs[j-1]))
		local q = complex(unpack(qs[j]))
		local q_1R = complex(unpack(qs[j+1]))
		local q_2R = complex(unpack(qs[j+2]))
		local d2q_dx2 = (-q_2L + 16 * q_1L - 30 * q + 16 * q_1R - q_2R) / (12 * dr * dr)
		local dq_dx = (q_2L - 8 * q_1L + 8 * q_1R - q_2R) / (12 * dr)
		--]]
		--[[
		local q = complex(unpack(qs[j]))
		local dq_dx = D(qs, j, 1, 4)
		local d2q_dx2 = D(qs, j, 2, 4)
		--]]
		
		local q2 = q:norm()
		local q4 = q2 * q2
		dq_dts[j] = {(i * (d2q_dx2 + 4 / r * dq_dx - q4 * q)):unpack()} 
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
