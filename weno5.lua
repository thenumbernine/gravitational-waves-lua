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
		d = {3/10, 3/5, 1/10},
	},
	fv = {
		c = {
			[3] = {
				{15/8, -5/4,  3/8},
				{ 3/8,  3/4, -1/8},
				{-1/8,  3/4,  3/8},
				{ 3/8, -5/4, 15/8},
			},
		},
		d = {5./ 16.,   5./  8.,   1./ 16.},
	},
}

local function sqr(x) return x*x end

local function weno5(v, i, fv_or_fd, l_or_r, weno5method, numWaves)
	local cofs = l_or_r == 'l' and 1 or 2
	local c = coeffs[fv_or_fd].c[3]
	local d = coeffs[fv_or_fd].d
	local d0 = l_or_r == 'l' and 3 or 1
	local dd = l_or_r == 'l' and -1 or 1
	
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
		local w
		if weno5method == '1996 Jiang Shu' then	
			w = {
				d[d0 + 0 * dd] / sqr(epsilon + beta[1]),
				d[d0 + 1 * dd] / sqr(epsilon + beta[2]),
				d[d0 + 2 * dd] / sqr(epsilon + beta[3]),
			}
		elseif weno5method == '2008 Borges' then
			local tau5 = math.abs(beta[1] - beta[3])
			w = {
				d[d0 + 0 * dd] * (1 + tau5 / (epsilon + beta[1])),
				d[d0 + 1 * dd] * (1 + tau5 / (epsilon + beta[2])),
				d[d0 + 2 * dd] * (1 + tau5 / (epsilon + beta[3])),
			}
		elseif weno5method == '2010 Shen Zha' then
			local epsilon = 1e-10
			local shen_zha_A = 50
			local minB = math.min(beta[1], beta[2], beta[3])
			local maxB = math.max(beta[1], beta[2], beta[3])
			local R0 = minB / (maxB + epsilon)
			beta[1] = beta[1] + R0 * shen_zha_A * minB;
			beta[2] = beta[2] + R0 * shen_zha_A * minB;
			beta[3] = beta[3] + R0 * shen_zha_A * minB;
			w = {
				(d[d0 + 0 * dd]) / sqr(epsilon + beta[1]),
				(d[d0 + 1 * dd]) / sqr(epsilon + beta[2]),
				(d[d0 + 2 * dd]) / sqr(epsilon + beta[3]),
			}
		end
		
		local vs = {
			c[cofs+0][1] * v[i+0][j] + c[cofs+0][2] * v[i+1][j] + c[cofs+0][3] * v[i+2][j],
			c[cofs+1][1] * v[i-1][j] + c[cofs+1][2] * v[i+0][j] + c[cofs+1][3] * v[i+1][j],
			c[cofs+2][1] * v[i-2][j] + c[cofs+2][2] * v[i-1][j] + c[cofs+2][3] * v[i+0][j],
		}
		
		local wsum = w[1] + w[2] + w[3]
		result[j] = (w[1] * vs[1] + w[2] * vs[2] + w[3] * vs[3]) / wsum
	end
	return result
end

return weno5
