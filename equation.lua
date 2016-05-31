local class = require 'ext.class'

local Equation = class()

Equation.State = require 'state' 

-- note this is to be used on the child class object
-- so 'self' is the subclass.  or the object.  either works.
function Equation:buildGraphInfos(getters)
	-- TODO give this functionality to Equation
	self.graphInfos = table()
	local w = math.ceil(math.sqrt(#getters))
	local h = math.ceil(#getters/w)
	local i,j = 0,0
	for _,getter in ipairs(getters) do
		local name = next(getter)
		local func = getter[name]
		self.graphInfos:insert{
			viewport = {i/w, j/h, 1/w, 1/h},
			getter = getter[name],
			name = name,
		}
		i = i + 1
		if i == w then
			i = 0
			j = j + 1
		end
	end
	self.graphInfoForNames = self.graphInfos:map(function(info,i)
		return info, info.name
	end)
end
-- future TODO: build this once all graphs are collected, 
-- so the graphs don't only have to match the first sim of the running set
-- further future TODO: everything with ImGUI, and open and close windows and stuff

local solverLinearFunc = require 'solverlinearfunc'
function Equation.buildField(matrixField)
	local linearFunc = solverLinearFunc(matrixField)
	return function(self, ...)	-- self isn't needed for linear systems.  but it is needed for some subclasses. 
		return linearFunc(...)
	end
end

--[[
default implementation will dot with j'th row of eigenvectorsInverse[i]
subclasses with sparse matrices (like ADM) will be able to override this and optimize away (those 37x37 matrices)

another note: eigenfields never have input vectors.  they are made of state vaules, and their input is state values, so there's no need to define an inner product.
...except the fact that some of the state variables are on the i'th entry, and some are of the i+1/2'th entry...
--]]
Equation.fluxTransform = Equation.buildField'fluxMatrix'
Equation.applyLeftEigenvectors = Equation.buildField'eigenvectorsInverse'
Equation.applyRightEigenvectors = Equation.buildField'eigenvectors'

function Equation:calcMinMaxEigenvaluesFromCons(...)
	return firstAndLast(self:calcEigenvaluesFromCons(...))
end

-- functions that use sim:

-- used by SolverFV
function Equation:calcCellMinMaxEigenvalues(sim, i)
	return self:calcMinMaxEigenvaluesFromCons(table.unpack(sim.qs[i]))
end


return Equation
