local class = require 'ext.class'

local Equation = class()

-- note this is to be used on the child class object
-- so 'self' is the subclass.  or the object.  either works.
function Equation:buildGraphInfos(getters)
	-- TODO give this functionality to Equation
	self.graphInfos = table()
	for _,getter in ipairs(getters) do
		local name = next(getter)
		local func = assert(getter[name], "failed to get graphInfo for "..name)
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
apply any arbitrary left/right eigenvector matrix, or flux to any vector 
from/true = true if transforming from states, false if transforming from waves
from, to:
	true, false = left transform
	false, true = right transform
	true, true = flux transform
--]]
function Equation:eigenTransform(solver, m, v, from, to)
	local matrixWidth = from and solver.numStates or solver.numWaves
	local matrixHeight = to and solver.numStates or solver.numWaves
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

-- helpers for the above from/to 
-- any arbitrary left eigenvector matrix to any vector 
function Equation:eigenLeftTransform(solver, m, v)
	return self:eigenTransform(solver, m, v, true, false)
end
-- any arbitrary right eigenvector matrix to any vector 
function Equation:eigenRightTransform(solver, m, v)
	return self:eigenTransform(solver, m, v, false, true)
end
-- any arbitrary flux matrix to any vector
function Equation:fluxMatrixTransform(solver, m, v)
	return self:eigenTransform(solver, m, v, true, true)
end

-- default implementation is arithmetic
function Equation:calcRoeValues(qL, qR)
	local q = {}
	for i=1,self.numStates do
		q[i] = .5 * (qL[i] + qR[i])
	end
	return table.unpack(q)
end

-- default implementation just treats L and R the same
function Equation:calcCellCenterRoeValues(solver, i)
	local q = solver.qs[i]
	return self:calcRoeValues(q, q)
end

-- functions that use sim:

-- used by Roe, passed to the input arguments calcEigenBasis
-- calculates the values used by calcEigenBasis to compute the eigenvalues, eigenvectors, and dF/dU matrices
function Equation:calcInterfaceRoeValues(solver, i)
	return self:calcRoeValues(solver:get_qL(i), solver:get_qR(i))
end

-- used by SolverFV
function Equation:calcCellMinMaxEigenvalues(sim, i)
	if self.calcMinMaxEigenvaluesFromCons then
		return self:calcMinMaxEigenvaluesFromCons(table.unpack(sim.qs[i]))
	else
		return firstAndLast(self:calcEigenvaluesFromCons(table.unpack(sim.qs[i])))
	end
end

return Equation
