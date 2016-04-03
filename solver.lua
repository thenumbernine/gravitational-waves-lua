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
	self.t0 = args.t or 0
	self.iteration = args.iteration or 0
	self.cfl = args.cfl or .5
	self.fixed_dt = args.fixed_dt
	self.stopAtTime = args.stopAtTime

	self.usePPM = args.usePPM

	self.xs = {}
	self.ixs = {}
	self.qs = self:newState()
end

function Solver:newState()
	return self.equation.State(self.gridsize, self.numStates)
end

function Solver:reset()
	self.t = self.t0
	
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

function Solver:applyBoundary()
	self.boundaryMethod(self.qs)
end

-- [[ PPM hack
function Solver:getPPM(xi,j)
	if not self.ppm_qRs then return end
	for i=1,#self.ixs-1 do
		--self.ixs[i] <= xi < self.ixs[i+1]
		if self.ixs[i+1] > xi then
			local x = (xi - self.ixs[i]) / (self.ixs[i+1] - self.ixs[i])
			-- DeltaA[j] = aR[j] - aL[j]
			local DeltaA = self.ppm_qRs[i][j] - self.ppm_qLs[i][j]
			-- a6[j] = 6 * (a[j] - .5 * (aL[j] + aR[j]))
			local a6 = 6 * (self.qs[i][j] - .5 * (self.ppm_qLs[i][j] + self.ppm_qRs[i][j]))
			-- a(x) = aL[j] + (x - ix[j]) / dx[j] * (DeltaA[j] + a6[j] * (1 - x))
			return self.ppm_qLs[i][j] + x * (DeltaA + a6 * (1 - x))
		end
	end	
end
--]]

function Solver:step(dt, getLeft, getRight)
	self:integrate(dt, function()
		local dq_dt = self:calcFlux(dt, getLeft, getRight)
		if self.equation.sourceTerm then
			dq_dt = dq_dt + self.equation:sourceTerm(self, self.qs)
		end
		return dq_dt
	end)
end

function Solver:iterate()
	self:applyBoundary()

	local getLeft, getRight
-- [[ my attempt at PPM
if self.usePPM then
	self.ppm_iqs = self.ppm_iqs or self:newState()
	-- left and right parabolic extrapolated values at each cell
	self.ppm_qLs = self.ppm_qLs or self:newState()
	self.ppm_qRs = self.ppm_qRs or self:newState()
	for i=3,self.gridsize-1 do
		for j=1,self.numStates do
			self.ppm_iqs[i][j] = (7*(self.qs[i-1][j] + self.qs[i][j]) - (self.qs[i+1][j] + self.qs[i-2][j]))/12
		end
	end
	for i=1,self.gridsize-2 do
		for j=1,self.numStates do
			self.ppm_qLs[i][j] = self.ppm_iqs[i][j]
			self.ppm_qRs[i][j] = self.ppm_iqs[i+1][j]
			-- [=[
			if (self.ppm_qRs[i][j] - self.qs[i][j]) * (self.qs[i][j] - self.ppm_qLs[i][j]) <= 0 then
				self.ppm_qLs[i][j] = self.qs[i][j]
				self.ppm_qRs[i][j] = self.qs[i][j]
			else
				local deltaRL = self.ppm_qRs[i][j] - self.ppm_qLs[i][j]
				local deltaRL2 = deltaRL * deltaRL  
				local upperBound = deltaRL2 / 6
				local lowerBound = -upperBound
				local samplePoint = deltaRL * (self.qs[i][j] - 1/2 * (self.ppm_qLs[i][j] + self.ppm_qRs[i][j]))
				if samplePoint > upperBound then
					self.ppm_qLs[i][j] = 3 * self.qs[i][j] - 2 * self.ppm_qRs[i][j]
				elseif lowerBound > samplePoint then
					self.ppm_qRs[i][j] = 3 * self.qs[i][j] - 2 * self.ppm_qLs[i][j]
				end
			end
			--]=]
		end
	end
	-- get left and right at each interface
	getLeft = function(sim,i) return sim.ppm_qRs[i-1] end
	getRight = function(sim,i) return sim.ppm_qLs[i] end
getLeft = function(sim,i) return sim.qs[i-1] end
getRight = function(sim,i) return sim.qs[i] end
else
--]]
-- [[ ordinary	
	getLeft = function(sim,i) return sim.qs[i-1] end
	getRight = function(sim,i) return sim.qs[i] end
--]]
end

	local dt = self:calcDT(getLeft, getRight)

	self:step(dt, getLeft, getRight)

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
