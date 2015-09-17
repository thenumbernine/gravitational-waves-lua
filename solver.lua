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
	self.t = 0
	self.cfl = .5
	self.xs = {}
	self.ixs = {}
	self.qs = self:newState()
end

Solver.State = require 'state' 

function Solver:newState()
	return self.State(self.gridsize, self.numStates)
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

function Solver:addSourceToDeriv()
	local dq_dts = self:newState()
	for i=1,self.gridsize do
		self:addSourceToDerivCell(dq_dts, i)
	end
	return dq_dts
end

function Solver:integrate(dt, dq_dts)
	self.qs = self.integrator:integrate(self.qs, dt, dq_dts)
end

function Solver:integrateFlux(dt)
	self:integrate(dt, function()
		return self:calcFlux(dt)
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
	local dt = self:calcDT()

	self:integrateFlux(dt)

	self:integrate(dt, function()
		return self:addSourceToDeriv()
	end)
	
	if self.postIterate then
		self:postIterate(dt)
	end
	
	self.t = self.t + dt
end

function Solver:addSourceToDerivCell() end

return Solver

