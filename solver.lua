local class = require 'ext.class'
local integrators = require 'integrators'
local matrix = require 'matrix'

-- base functions
local Solver = class()

-- TODO use this everywhere
Solver.numGhost = 2

function Solver:init(args)
	self.equation = assert(args.equation or self.equation)
	
	self.numStates = self.equation.numStates
	self.numWaves = self.equation.numWaves or self.numStates

	self.gridsize = assert(args.gridsize) + 2 * self.numGhost
	self.domain = assert(args.domain)
	self.boundaryMethod = assert(args.boundaryMethod)

	self.integrator = args.integrator or integrators.ForwardEuler()
	
	self.name = self.name .. ' int.=' .. self.integrator.name 

	self.t0 = args.t or 0
	self.iteration = args.iteration or 0
	self.cfl = args.cfl or .6
	self.fixed_dt = args.fixed_dt
	self.stopAtTimes = args.stopAtTimes

	self.xs = matrix()
	self.ixs = matrix()
	self.qs = self:newState()
end

function Solver:newState()
	return matrix.zeros(self.gridsize, self.numStates)
end

function Solver:reset()
	self.t = self.t0

	local xmin = self.domain.xmin
	local xmax = self.domain.xmax
	local width = xmax - xmin

	if not self.useXEdgeCentered then
		for i=1,self.gridsize do
			self.xs[i] = ((i-self.numGhost)-.5)/(self.gridsize-2*self.numGhost)*width + xmin
		end
		-- interface center of coordinate system
		for i=1,self.gridsize+1 do
			self.ixs[i] = ((i-self.numGhost)-1)/(self.gridsize-2*self.numGhost)*width + xmin
		end
	else
		local avg = .5 * (xmin + xmax)
		for i=1,self.gridsize do
			local f = math.floor(i-self.numGhost-1)/math.floor(self.gridsize-2*self.numGhost-1)
			self.xs[i] = f*width + xmin
			--self.xs[i] = f * xmax + (1-f) * xmin
			--self.xs[i] = avg + (f - .5) * width
		end
		-- interface center of coordinate system
		for i=1,self.gridsize+1 do
			local f = (i-self.numGhost-1.5)/(self.gridsize-2*self.numGhost-1)
			self.ixs[i] = f*width + xmin
			--self.ixs[i] = f * xmax + (1-f) * xmin
			--self.ixs[i] = avg + (f - .5) * width
		end
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
	self.boundaryMethod(self.qs, self.numGhost)
end

function Solver:step(dt)
	self:integrate(dt, function()
		self:applyBoundary()
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
--[[
print'\nself.UBufObj:'
local max = #tostring(self.gridsize)
for i=1,self.gridsize do
	io.write((' '):rep(max-#tostring(i)), i, ':')
	print(('\t[%.50f,'):format(self.qs[i][1]))
	print(('\t%.50f,'):format(self.qs[i][2]))
	print(('\t%.50f,'):format(0))
	print(('\t%.50f,'):format(0))
	print(('\t%.50f]'):format(self.qs[i][3]))
end
--]]	
	local dt = self:calcDT()

	self:step(dt)

	if self.postIterate then
		self:applyBoundary()
		self:postIterate(dt)
	end
	if self.equation.postIterate then
		self:applyBoundary()
		self.equation:postIterate(self, self.qs)
	end
	
	self.t = self.t + dt
	self.iteration = self.iteration + 1

	-- boundary before exiting the loop so rendering looks correct
	self:applyBoundary()
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
