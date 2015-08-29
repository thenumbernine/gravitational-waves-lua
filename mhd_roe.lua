local class = require 'ext.class'
local Roe = require 'roe'

local MHDRoe = class(Roe)

MHDRoe.equation = require 'mhd'()

-- a,b,c,d 3x1 columns, A = [a,b,c], solves Ax=d
local function cramerSolve(a,b,c,d)
end

-- A x = b
-- x = A^-1 b
-- returns x
local function linearSolve(A, b)
	local result = {unpack(b)}

	local n = #A

	local originalA = A
	do
		local cloneA = {}
		for i=1,n do
			assert(#A[i] == n, "exepcted A to be a square matrix")
			cloneA[i] = {unpack(A[i])}
		end
		A = cloneA
	end

	for i=1,n do
		if A[i][i] == 0 then
			-- pivot with a row beneath this one
			local found = false
			for j=i+1,n do
				if A[j][i] ~= 0 then
					for k=1,n do
						A[j][k], A[i][k] = A[i][k], A[j][k]
					end
					result[j], result[i] = result[i], result[j]
					found = true
					break
				end
			end
			if not found then
				error("couldn't find a row to pivot for matrix:\n"..table.map(originalA,
					function(row) return '[' .. table.concat(row, ',') .. ']' end
				):concat'\n')
			end
		end
		-- rescale diagonal
		if A[i][i] ~= 1 then
			-- rescale column
			local s = A[i][i]
			for j=1,n do
				A[i][j] = A[i][j] / s
			end
			result[i] = result[i] / s
		end
		-- eliminate columns apart from diagonal
		for j=1,n do
			if j ~= i then
				if A[j][i] ~= 0 then
					local s = A[j][i]
					for k=1,n do
						A[j][k] = A[j][k] - s * A[i][k]
					end
					result[j] = result[j] - s * result[i]
				end
			end
		end
	end

	return result
end

--[[
-- testing
print(unpack(linearSolve(
	-- store row-major so Lua indexing matches math indexing
	{
		{1,0,0},
		{0,1,2},
		{0,1,0},
	},
		{1,2,3}
)))
os.exit()
--]]

function MHDRoe:eigenfields(i, v)

	-- eigenvector transform:
	-- y = R*v; return y
	-- inverse transform:
	-- x = R^-1 * v; return x;
	-- v = R*x
	-- solve linear system above, use right-eigenvectors as matrix and v as solution
	return linearSolve(self.eigenvectors[i], v)

--[=[ me trying to decypher Brio & Wu's paper

	local qL = self.qs[i-1]
	local qR = self.qs[i]

-- begin block similar to calcInterfaceEigenBasis
	local gamma = self.gamma	
	local gammaMinusOne = gamma - 1
	
	local rhoL, uxL, uyL, uzL, BxL, ByL, BzL, pL = self.equation:stateToPrims(unpack(qL))
	local rhoR, uxR, uyR, uzR, BxR, ByR, BzR, pR = self.equation:stateToPrims(unpack(qR))

	local rho = .5*(rhoL + rhoR)
	local ux = .5*(uxL + uxR)
	local uy = .5*(uyL + uyR)
	local uz = .5*(uzL + uzR)
	local Bx = .5*(BxL + BxR)
	local By = .5*(ByL + ByR)
	local Bz = .5*(BzL + BzR)
	local p = .5*(pL + pR)

	local mu = self.equation.mu
	local BSq = Bx*Bx + By*By + Bz*Bz
	local uSq = ux*ux + uy*uy + uz*uz
	local BdotU = Bx*ux + By*uy + Bz*uz
	local vaxSq = Bx*Bx / (mu*rho)
	local vax = sqrt(vaxSq)
	local CsSq = gamma*p / rho
	local Cs = sqrt(CsSq)
	local vaSq = BSq / (mu*rho)
	local va = sqrt(vaSq)
	local cStarSq = vaSq + CsSq
	local discr = sqrt(cStarSq*cStarSq - 4 * vaxSq*CsSq)
	local vfSq = .5 * (cStarSq + discr)
	local vsSq = .5 * (cStarSq - discr)
	local vf = sqrt(vfSq)
	local vs = sqrt(vsSq)
	local sgnBx = Bx >= 0 and 1 or -1

	local ETotal = sim.qs[i][8]
	local PStar = p + .5 * BSq
	local H = (ETotal + PStar) / rho
-- end block

	-- transform v by eigenbasis inverse
	--[[
	Brio & Wu 51 & 52
	--]]
	local d547 = cramerSolve(
		-- column 1
		{	alpha_f * vf,
			-beta_y * alpha_s * b_x * vf,
			-beta_z * alpha_s * b_x * vf	},
		-- column 2
		{	alpha_s * vs,
			beta_y * alpha_s * sgnBx / vf,
			beta_z * alpha_s * sgnBx / vf	},
		-- column 3
		{	0,
			-beta_z * sgnBx,
			beta_y * sgnBx	},
		-- soln - should be formed from the input vector (i.e. the L and R states, not the interface state)
		{	v[2]/v[1] * rhoStar,
			v[3]/v[1] * rhoStar,
			v[4]/v[1] * rhoStar	}
	)
	local d126 = cramerSolve(
		-- column 1
		{	alpha_s * beta_y * vfSq / sqrtRho,
			alpha_s * beta_z * vfSq / sqrtRho,
			alpha_f * vfSq / (gamma - 1) + alpha_s * (gamma - 2) / (gamma - 1) * vfSq	},
		-- column 2
		{	-alpha_f * a^2 * beta_y / (vfSq * sqrtRho),
			-alpha_f * a^2 * beta_z / (vfSq * sqrtRho),
			alpha_s * bx^2 * a^2 / (vfSq * (gamma - 1)) + alpha_f * (gamma - 2) / (gamma - 1) * a^2 / vfSq	},
		-- column 3
		{	beta_z / sqrtRho,
			-beta_y / sqrtRho,
			0	},
		-- soln
		{	v[5],
			v[6],
			((v[8] - .5/mu * (v[5]^2 + v[6]^2 + v[7]^2)) / v[1] - .5 * (v[2]^2 + v[3]^2 + v[4]^2) / v[1]^2) / (gamma - 1)^2 - .5 * (v[5]^2 + v[6]^2)	}
	)

	local c4 = Lambda * rho - alpha_f * d126[1] - alpha_s * d126[2]

	c7 = (d1 + d5) / 2
	c1 = (d1 - d5) / 2
	c5 = (d2 + d4) / 2
	c3 = (d2 - d4) / 2
	c2 = (d6 + d7) / 2
	c6 = (d6 - d7) / 2

	return {c1, c2, c3, c4, c5, c6, c7}
--]=]
end

return MHDRoe
