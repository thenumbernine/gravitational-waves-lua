local mat33det = function(A)
	return A[1][1] * A[2][2] * A[3][3]
		+ A[2][1] * A[3][2] * A[1][3]
		+ A[3][1] * A[1][2] * A[2][3]
		- A[3][1] * A[2][2] * A[1][3]
		- A[2][1] * A[1][2] * A[3][3]
		- A[1][1] * A[3][2] * A[2][3]
end

local function mat33inv(A,det)
	det = det or mat33det(A)
	if det == 0 then
		print('A=')
		for j=1,3 do
			print(table.concat(A[j], ', '))
		end
		error 'singular!'
	end
	local invDet = 1 / det
	local V = {{},{},{}}
	for j=1,3 do
		local j1 = j % 3 + 1 
		local j2 = j1 % 3 + 1
		for k=1,3 do
			local k1 = k % 3 + 1
			local k2 = k1 % 3 + 1
			V[k][j] = invDet * (A[j1][k1] * A[j2][k2] - A[j1][k2] * A[j2][k1])
		end
	end
	return V
end

return {
	det = mat33det,
	inv = mat33inv,
}
