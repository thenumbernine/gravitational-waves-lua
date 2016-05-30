local function solverLinearFunc(matrixField)
	return function(solver, i, v)
		local m = solver[matrixField][i]
		local result = {}
		for j=1,solver.numStates do
			local sum = 0
			for k=1,solver.numStates do
				sum = sum + m[j][k] * v[k]
			end
			result[j] = sum
		end
		return result 
	end
end

return solverLinearFunc
