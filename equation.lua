local class = require 'ext.class'

local Equation = class()

function Equation:init()
	-- create transform functions based on class members

	--[[
	default implementation will dot with j'th row of eigenvectorsInverse[i]
	subclasses with sparse matrices (like ADM) will be able to override this and optimize away (those 37x37 matrices)

	another note: eigenfields never have input vectors.  they are made of state vaules, and their input is state values, so there's no need to define an inner product.
	...except the fact that some of the state variables are on the i'th entry, and some are of the i+1/2'th entry...
	--]]
	self.fluxTransform = self.createTransformFunc('fluxMatrix', true, true)
	self.applyLeftEigenvectors = self.createTransformFunc('eigenvectorsInverse', true, false)
	self.applyRightEigenvectors = self.createTransformFunc('eigenvectors', false, true)
end

Equation.State = require 'state' 

-- note this is to be used on the child class object
-- so 'self' is the subclass.  or the object.  either works.
function Equation:buildGraphInfos(getters)
	-- TODO give this functionality to Equation
	self.graphInfos = table()
	for _,getter in ipairs(getters) do
		local name = next(getter)
		local func = getter[name]
		self.graphInfos:insert{
			getter = getter[name],
			name = name,
		}
	end
	self.graphInfoForNames = self.graphInfos:map(function(info,i)
		return info, info.name
	end)
end

--[[
from = whether we are transforming from state vector (true) or wave vector (false)
to = same
--]]
function Equation.createTransformFunc(matrixField, from, to)
	return function(self, solver, i, v)
		local matrixWidth = from and solver.numStates or solver.numWaves
		local matrixHeight = to and solver.numStates or solver.numWaves
		local m = solver[matrixField][i]
		local result = {}
		for j=1,matrixHeight do
			local sum = 0
			for k=1,matrixWidth do
				sum = sum + m[j][k] * v[k]
			end
			result[j] = sum
		end
		return result 
	end
end

-- this 'state' vs 'cons' distinction only exists because I was messing with the 'euler1dquasilinear' which was stupid
function Equation:calcEigenvaluesFromState(...)
	return self:calcEigenvaluesFromCons(self:calcConsFromState(...))
end
function Equation:calcMinMaxEigenvaluesFromState(...)
	if self.calcMinMaxEigenvaluesFromCons then
		return self:calcMinMaxEigenvaluesFromCons(self:calcConsFromState(...))
	else
		return firstAndLast(self:calcEigenvaluesFromState(...))
	end
end

-- functions that use sim:

-- used by SolverFV
function Equation:calcCellMinMaxEigenvalues(sim, i)
	return self:calcMinMaxEigenvaluesFromState(table.unpack(sim.qs[i]))
end

-- by default, assume the state is conservative variables
function Equation:calcConsFromState(...) return ... end

return Equation
