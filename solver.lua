local class = require 'ext.class'
local integrators = require 'integrators'

-- base functions
local Solver = class()

function Solver:init(args)
	self.equation = assert(args.equation or self.equation)
	
	self.numStates = self.equation.numStates
	self.numWaves = self.equation.numWaves or self.numStates

	self.gridsize = assert(args.gridsize)
	self.domain = assert(args.domain)
	self.boundaryMethod = assert(args.boundaryMethod)

	self.integrator = args.integrator or integrators.ForwardEuler()
	
	self.name = self.name .. ' int.=' .. self.integrator.name 

	self.t0 = args.t or 0
	self.iteration = args.iteration or 0
	self.cfl = args.cfl or .5
	self.fixed_dt = args.fixed_dt
	self.stopAtTimes = args.stopAtTimes

	self.xs = {}
	self.ixs = {}
	self.qs = self:newState()
end

local matrix = require 'matrix'
function Solver:newState()
	return matrix.zeros(self.gridsize, self.numStates)
end

function Solver:reset()
	self.t = self.t0
	
	for i=1,self.gridsize do
		self.xs[i] = ((i-2)-.5)/(self.gridsize-4)*(self.domain.xmax - self.domain.xmin) + self.domain.xmin
	end
	-- interface center of coordinate system
	for i=1,self.gridsize+1 do
		self.ixs[i] = ((i-2)-1)/(self.gridsize-4)*(self.domain.xmax - self.domain.xmin) + self.domain.xmin
	end

	-- state at cell centers
	for i=1,self.gridsize do
		self.qs[i] = self.equation:initCell(self,i)
	end
end

function Solver:integrate(dt, dq_dts)
	self.qs = self.integrator:integrate(self.qs, dt, dq_dts)
end

function Solver:applyBoundary()
	self.boundaryMethod(self.qs)
end

function Solver:step(dt)
	self:integrate(dt, function()
		local dq_dt = self:calcDerivFromFluxes(dt)
		if self.equation.sourceTerm then
			dq_dt = dq_dt + self.equation:sourceTerm(self, self.qs, dt)
		end
		return dq_dt
	end)
end

function Solver:calcDT()
	return assert(self.fixed_dt)
end

function Solver:iterate()
	self:applyBoundary()

	local dt = self:calcDT()
	self:step(dt)

	if self.postIterate then
		self:postIterate(dt)
	end
	if self.equation.postIterate then
		self.equation:postIterate(self, self.qs)
	end
	
	self.t = self.t + dt
	self.iteration = self.iteration + 1
end

-- get the q at the left side of the interface
function Solver:get_qL(i)
	return self.qs[i-1]
end

-- get the q at the right side of the interface
function Solver:get_qR(i)
	return self.qs[i]
end


return Solver
