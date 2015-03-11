require 'ext'
local calcFluxSchemes = require 'scheme'

-- base functions
local Simulation = class()

function Simulation:init(args)
	args = args or {}
	self.gridsize = assert(args.gridsize)
	self.domain = assert(args.domain)
	self.boundaryMethod = assert(args.boundaryMethod)
	self.slopeLimiter = assert(args.slopeLimiter)

	self.calcFlux = calcFluxSchemes.Roe 
	self.t = 0
	self.cfl = .5
	self.xs = {}
	self.ixs = {}
	self.qs = {}
	self.deltaQTildes = {}
	self.fluxes = {}
	self.dq_dts = {}
	-- used by Roe
	self.fluxMatrix = {}
	self.eigenvalues = {}
	self.eigenvectors = {}
	self.eigenvectorsInverse = {}
	self.eigenbasisErrors = {}
	self.fluxMatrixErrors = {}
end

local function buildField(matrixField)
	return function(self, i, v)
		local m = self[matrixField][i]
		local result = {}
		for j=1,self.numStates do
			local sum = 0
			for k=1,self.numStates do
				sum = sum + m[j][k] * v[k]
			end
			result[j] = sum
		end
		return result 
	end
end

--[[
default implementation will dot with j'th row of eigenvectorsInverse[i]
subclasses with sparse matrices (like ADM) will be able to override this and optimize away (those 37x37 matrices)

another note: eigenfields never have input vectors.  they are made of state vaules, and their input is state values, so there's no need to define an inner product.
...except the fact that some of the state variables are on the i'th entry, and some are of the i+1/2'th entry...
--]]
Simulation.fluxTransform = buildField'fluxMatrix'
Simulation.eigenfields = buildField'eigenvectorsInverse'
Simulation.eigenfieldsInverse = buildField'eigenvectors'

function Simulation:reset()
	for i=1,self.gridsize do
		self.xs[i] = (i-.5)/self.gridsize*(self.domain.xmax - self.domain.xmin) + self.domain.xmin
	end
	-- interface center of coordinate system
	for i=1,self.gridsize+1 do
		self.ixs[i] = (i-1)/self.gridsize*(self.domain.xmax - self.domain.xmin) + self.domain.xmin
	end

	-- state at cell centers
	for i=1,self.gridsize do
		self.qs[i] = self:initCell(i)
	end

	-- state interfaces
	for i=1,self.gridsize+1 do
		self.deltaQTildes[i] = {}
		self.fluxes[i] = {}
		for j=1,self.numStates do
			self.deltaQTildes[i][j] = 0
			self.fluxes[i][j] = 0
		end
	end

	-- integration vars
	for i=1,self.gridsize do
		self.dq_dts[i] = {}
		for j=1,self.numStates do
			self.dq_dts[i][j] = 0
		end
	end
	
	for i=1,self.gridsize+1 do
		self.fluxMatrix[i] = {}
		self.eigenvalues[i] = {}
		self.eigenvectors[i] = {}
		self.eigenvectorsInverse[i] = {}
		for j=1,self.numStates do
			self.fluxMatrix[i][j] = {}
			self.eigenvectors[i][j] = {}
			self.eigenvectorsInverse[i][j] = {}
		end
		self.eigenbasisErrors[i] = 0
		self.fluxMatrixErrors[i] = 0
	end
end

function Simulation:zeroDeriv(dq_dts)
	-- zero deriv
	for i=1,self.gridsize do
		for j=1,self.numStates do
			dq_dts[i][j] = 0
		end
	end
end

function Simulation:addSourceToDeriv(dq_dts)
	for i=1,self.gridsize do
		self:addSourceToDerivCell(dq_dts, i)
	end
end

function Simulation:addFluxToDeriv(dq_dts)
	-- 5) integrate
	for i=1,self.gridsize do
		local dx = self.ixs[i+1] - self.ixs[i]
		for j=1,self.numStates do
			dq_dts[i][j] = dq_dts[i][j] - (self.fluxes[i+1][j] - self.fluxes[i][j]) / dx
		end
	end
end

function Simulation:integrateDeriv(dq_dts, dt)
	for i=1,self.gridsize do
		for j=1,self.numStates do
			self.qs[i][j] = self.qs[i][j] + dt * dq_dts[i][j]	
		end
	end
	self.t = self.t + dt
end

function Simulation:iterate()
	self:boundaryMethod()
	
	local dt = self:calcFlux()

	-- TODO create this up front, and create as many as needed for a particular integrator
	local dq_dts = {}
	for i=1,self.gridsize do
		dq_dts[i] = {}
		for j=1,self.numStates do
			dq_dts[i][j] = 0
		end
	end

	self:zeroDeriv(dq_dts)
	self:addSourceToDeriv(dq_dts)
	self:addFluxToDeriv(dq_dts)
	self:integrateDeriv(dq_dts, dt)
end

function Simulation:addSourceToDerivCell() end

return Simulation

