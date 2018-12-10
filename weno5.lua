local matrix = require 'matrix'

local coeffs = {
	fd = {
		-- 1998 Shu, table 2.1, sc[k][r][j+1]
		c = {
			[3] = {
				{11/6, -7/6,  1/3},
				{ 1/3,  5/6, -1/6},
				{-1/6,  5/6,  1/3},
				{ 1/3, -7/6, 11/6},
			},
		},
		-- 1998 Shu, between eqns 2.54 and 2.55
		-- d[r+1]
		d = {
			3/10, 3/5, 1/10, 
		},
	},
	fv = {
		c = {
			{15/8, -5/4,  3/8},
			{ 3/8,  3/4, -1/8},
			{-1/8,  3/4,  3/8},
			{ 3/8, -5/4, 15/8},
		},
		d = {
			l = {   1./ 16.,   5./  8.,   5./ 16. },
			r = {   5./ 16.,   5./  8.,   1./ 16. },
		},
	},
}

local function sqr(x) return x*x end

local function weno5(v, i, fv_or_fd, l_or_r, numWaves)
	local cofs = l_or_r == 'l' and 1 or 2
	local c = coeffs[fv_or_fd].c[3]
	local d = coeffs[fv_or_fd].d
	
	local epsilon = 1e-14
	local result = matrix()
	for j=1,numWaves do
		-- 1998 Shu eqn 2.63
		-- cites 1996 Jiang, Shu for the coefficient source
		local beta = {
			(13/12) * sqr( v[i+0][j] - 2 * v[i+1][j] + v[i+2][j]) + (1/4) * sqr(3 * v[i+0][j] - 4 * v[i+1][j] +     v[i+2][j]),
			(13/12) * sqr( v[i-1][j] - 2 * v[i+0][j] + v[i+1][j]) + (1/4) * sqr(    v[i-1][j]                 -     v[i+1][j]),
			(13/12) * sqr( v[i-2][j] - 2 * v[i-1][j] + v[i+0][j]) + (1/4) * sqr(    v[i-2][j] - 4 * v[i-1][j] + 3 * v[i+0][j]),
		}

		-- 1998 Shu, eqn 2.59
		local alpha = l_or_r == 'l' 
			and {
				d[3] / sqr(epsilon + beta[1]),
				d[2] / sqr(epsilon + beta[2]),
				d[1] / sqr(epsilon + beta[3]),
			}
			or {
				d[1] / sqr(epsilon + beta[1]),
				d[2] / sqr(epsilon + beta[2]),
				d[3] / sqr(epsilon + beta[3]),
			}
		
		local alphasum = alpha[1] + alpha[2] + alpha[3]
		
		local vs = {
			c[cofs+0][1] * v[i+0][j] + c[cofs+0][2] * v[i+1][j] + c[cofs+0][3] * v[i+2][j],
			c[cofs+1][1] * v[i-1][j] + c[cofs+1][2] * v[i+0][j] + c[cofs+1][3] * v[i+1][j],
			c[cofs+2][1] * v[i-2][j] + c[cofs+2][2] * v[i-1][j] + c[cofs+2][3] * v[i+0][j],
		}
		
		result[j] = (alpha[1] * vs[1] + alpha[2] * vs[2] + alpha[3] * vs[3]) / alphasum
	end
	return result
end

return weno5
