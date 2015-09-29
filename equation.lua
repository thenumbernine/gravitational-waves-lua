local class = require 'ext.class'

local Equation = class()

local function buildField(matrixField)
	return function(self, sim, i, v)
		local m = sim[matrixField][i]
		local result = {}
		for j=1,sim.numStates do
			local sum = 0
			for k=1,sim.numStates do
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
Equation.fluxTransform = buildField'fluxMatrix'
Equation.eigenfields = buildField'eigenvectorsInverse'
Equation.eigenfieldsInverse = buildField'eigenvectors'

return Equation
