local mat33det = function(U)
	return U[1][1] * U[2][2] * U[3][3]
		+ U[2][1] * U[3][2] * U[1][3]
		+ U[3][1] * U[1][2] * U[2][3]
		- U[3][1] * U[2][2] * U[1][3]
		- U[2][1] * U[1][2] * U[3][3]
		- U[1][1] * U[3][2] * U[2][3]
end

local function mat33inv(U,det)
	det = det or mat33det(U)
	if det == 0 then
		for j=1,3 do
			for k=1,3 do
				print('A('+i+','+j+') = '+U[j][k])
			end
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
			V[k][j] = invDet * (U[j1][k1] * U[j2][k2] - U[j1][k2] * U[j2][k1])
		end
	end
	return V
end

return {
	det = mat33det,
	inv = mat33inv,
}
