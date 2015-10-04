local class = require 'ext.class'
local integrators = require 'integrators'

-- base functions
local Solver = class()

function Solver:init(args)
	self.equation = assert(args.equation or self.equation)
	
	self.numStates = self.equation.numStates
	
	self.gridsize = assert(args.gridsize)
	self.domain = assert(args.domain)
	self.boundaryMethod = assert(args.boundaryMethod)
	self.fluxLimiter = assert(args.fluxLimiter)

	self.integrator = args.integrator or integrators.ForwardEuler()
	self.t = args.t or 0
	self.iteration = args.iteration or 0
	self.cfl = args.cfl or .5
	self.fixed_dt = args.fixed_dt
	self.stopAtTime = args.stopAtTime
	
	self.xs = {}
	self.ixs = {}
	self.qs = self:newState()
end

function Solver:newState()
	return self.equation.State(self.gridsize, self.numStates)
end

function Solver:reset()
	for i=1,self.gridsize do
		self.xs[i] = (i-.5)/self.gridsize*(self.domain.xmax - self.domain.xmin) + self.domain.xmin
	end
	-- interface center of coordinate system
	for i=1,self.gridsize+1 do
		self.ixs[i] = (i-1)/self.gridsize*(self.domain.xmax - self.domain.xmin) + self.domain.xmin
	end

	-- state at cell centers
	for i=1,self.gridsize do
		self.qs[i] = self.equation:initCell(self,i)
	end
end

function Solver:integrate(dt, dq_dts)
	self.qs = self.integrator:integrate(self.qs, dt, dq_dts)
end

function Solver:integrateFlux(dt, getLeft, getRight)
	self:integrate(dt, function()
		return self:calcFlux(dt, getLeft, getRight)
	end)
end

function Solver:applyBoundary()
	self.boundaryMethod(self.qs)
end

function Solver:iterate()
	self:applyBoundary()
--[[
print('new iter:')
for i=1,self.gridsize do
	io.write(self.xs[i])
	for j=1,self.numStates do
		io.write('\t',self.qs[i][j])
	end
	print()
end
--]]
	local getLeft = function(sim,i) return sim.qs[i-1] end
	local getRight = function(sim,i) return sim.qs[i] end
	
	local dt = self:calcDT(getLeft, getRight)

	self:integrateFlux(dt, getLeft, getRight)

	if self.equation.sourceTerm then
		self:integrate(dt, function()
			return self.equation:sourceTerm(self, self.qs)
		end)
	end

	if self.postIterate then
		self:postIterate(dt)
	end
	if self.equation.postIterate then
		self.equation:postIterate(self, self.qs)
	end

	self.t = self.t + dt
	self.iteration = self.iteration + 1
end

return Solver

