require 'ext'
local schemes = require 'scheme'


local State = class()

function State:init(h, w)
	for i=1,h do
		self[i] = {}
		for j=1,w do
			self[i][j] = 0
		end
	end
end

function State.__add(a,b)
	local c = State(#a, #a[1])
	for i=1,#a do
		for j=1,#a[1] do
			c[i][j] = a[i][j] + b[i][j]
		end
	end
	return c
end

function State.__mul(a,b)
	local function is(x) return type(x) == 'table' and x.isa and x:isa(State) end
	local src = is(a) and a or b
	local c = State(#src, #src[1])
	for i=1,#src do
		for j=1,#src[1] do
			local aij = type(a) == 'number' and a or a[i][j]
			local bij = type(b) == 'number' and b or b[i][j]
			c[i][j] = aij * bij
		end
	end
	return c
end

-- base functions
local Simulation = class()

function Simulation:init(args)
	args = args or {}
	self.gridsize = assert(args.gridsize)
	self.domain = assert(args.domain)
	self.boundaryMethod = assert(args.boundaryMethod)
	self.slopeLimiter = assert(args.slopeLimiter)

	self.scheme = args.scheme or schemes.Roe 
	self.t = 0
	self.cfl = .5
	self.xs = {}
	self.ixs = {}
	self.qs = self:newState()
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

Simulation.State = State

function Simulation:newState()
	return self.State(self.gridsize, self.numStates)
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
		self.qs[i] = self:initCell(i) -- adm1d3var requires direct assignment here
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
			dq_dts[i][j] = dq_dts[i][j] or 0
		end
	end
end

function Simulation:addSourceToDeriv()
	local dq_dts = self:newState()
	if self.scheme.addSourceToDeriv then
		self.scheme.addSourceToDeriv(self, dq_dts)
	end
	for i=1,self.gridsize do
		self:addSourceToDerivCell(dq_dts, i)
	end
	return dq_dts
end

function Simulation:integrate(dt, dq_dts)
	-- [[ Euler
	self.qs = self.qs + dt * dq_dts()
	--]]

	--[[ RK4
	local k1 = dq_dts(self.qs)
	local k2 = dq_dts(self.qs + .5 * k1)
	local k3 = dq_dts(self.qs + .5 * k2)
	local k4 = dq_dts(self.qs + k3)
	self.qs = self.qs + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
	--]]
end

function Simulation:iterate()
	self:boundaryMethod()
	
	local dt = self.scheme.calcDT(self)

	self:integrate(dt, function()
		return self.scheme.calcFlux(self, dt)
	end)

	self:integrate(dt, function()
		return self:addSourceToDeriv()
	end)
	
	if self.scheme.postIterate then
		self.scheme.postIterate(self, dt)
	end
	
	self.t = self.t + dt
end

function Simulation:addSourceToDerivCell() end

return Simulation

