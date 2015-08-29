local class = require 'ext.class'

local function isnan(x) return x ~= x end
local function isinf(x) return x == math.huge or x == -math.huge end
local function finite(x) return not isnan(x) and not isinf(x) end

local MHD = class()

MHD.numStates = 8
MHD.gamma = 5/3	
MHD.mu = 1

do
	local rho = function(self,i) return self.qs[i][1] end
	local ux = function(self,i) return self.qs[i][2]/rho(self,i) end
	local uy = function(self,i) return self.qs[i][3]/rho(self,i) end
	local uz = function(self,i) return self.qs[i][4]/rho(self,i) end
	local Bx = function(self,i) return self.qs[i][5]/self.equation.mu end
	local By = function(self,i) return self.qs[i][6]/self.equation.mu end
	local Bz = function(self,i) return self.qs[i][7]/self.equation.mu end
	local ETotal = function(self,i) return self.qs[i][8] end
	local eTotal = ETotal / rho
	local EMag = .5*(Bx*Bx + By*By + Bz*Bz) / function(self,i) return self.equation.mu end
	local EHydro = ETotal - EMag
	local eHydro = EHydro / rho
	local eKin = function(self,i) return .5*(ux(self,i)^2 + uy(self,i)^2 + uz(self,i)^2) end
	local eInt = eHydro - eKin
	local p = function(self,i) return eInt(self,i) / (self.equation.gamma - 1) end
	MHD.graphInfos = table{
		{viewport={0/5, 0/4, 1/5, 1/4}, getter=function(self,i) return self.eigenbasisErrors[i] end, name='eigenbasis error', color={1,0,0}, range={-30, 30}},
		{viewport={0/5, 1/4, 1/5, 1/4}, getter=function(self,i) return self.fluxMatrixErrors[i] end, name='reconstruction error', color={1,0,0}, range={-30, 30}},
		{viewport={1/5, 0/4, 1/5, 1/4}, getter=function(self,i) return rho(self,i) end, name='rho', color={1,0,1}},
		{viewport={1/5, 1/4, 1/5, 1/4}, getter=function(self,i) return p(self,i) end, name='p', color={1,0,1}},
		{viewport={2/5, 0/4, 1/5, 1/4}, getter=function(self,i) return ux(self,i) end, name='ux', color={0,1,0}},
		{viewport={2/5, 1/4, 1/5, 1/4}, getter=function(self,i) return uy(self,i) end, name='uy', color={0,1,0}},
		{viewport={2/5, 2/4, 1/5, 1/4}, getter=function(self,i) return uz(self,i) end, name='uz', color={0,1,0}},
		{viewport={3/5, 0/4, 1/5, 1/4}, getter=function(self,i) return Bx(self,i) end, name='Bx', color={.5,.5,1}},
		{viewport={3/5, 1/4, 1/5, 1/4}, getter=function(self,i) return By(self,i) end, name='By', color={.5,.5,1}},
		{viewport={3/5, 2/4, 1/5, 1/4}, getter=function(self,i) return Bz(self,i) end, name='Bz', color={.5,.5,1}},
		{viewport={4/5, 0/4, 1/5, 1/4}, getter=function(self,i) return eTotal(self,i) end, name='eTotal', color={1,1,0}},
		{viewport={4/5, 1/4, 1/5, 1/4}, getter=function(self,i) return eKin(self,i) end, name='eKin', color={1,1,0}},
		{viewport={4/5, 2/4, 1/5, 1/4}, getter=function(self,i) return eInt(self,i) end, name='eInt', color={1,1,0}},
		{viewport={4/5, 3/4, 1/5, 1/4}, getter=function(self,i) return eHydro(self,i) end, name='eHydro', color={1,1,0}},
		{viewport={3/5, 3/4, 1/5, 1/4}, getter=function(self,i) return EMag(self,i) end, name='EMag', color={1,1,0}},
	}
end
MHD.graphInfoForNames = MHD.graphInfos:map(function(info,i)
	return info, info.name
end)

function MHD:initCell(sim,i)
	local x = sim.xs[i]
	local rho = x < 0 and .1 or 1
	local ux, uy, uz = 0, 0, 0
	local Bx = 0	-- .75
	local By = 0	-- x < 0 and 1 or -1
	local Bz = 0
	local p = 1		-- x < 0 and 1 or .1
	local eInt = p / (self.gamma - 1)
	local eKin = .5*(ux*ux + uy*uy + uz*uz)
	local EMag = .5*(Bx*Bx + By*By + Bz*Bz) / self.mu
	local ETotal = rho*(eInt + eKin) + EMag 
	return {rho, rho*ux, rho*uy, rho*uz, Bx, By, Bz, ETotal}
end

function MHD:stateToPrims(rho, mx, my, mz, Bx, By, Bz, ETotal)
	local ux, uy, uz = mx / rho, my / rho, mz / rho
	local EMag = .5*(Bx*Bx + By*By + Bz*Bz) / self.mu
	local EHydro = ETotal - EMag
	local eHydro = EHydro / rho
	local eKin = .5*(ux*ux + uy*uy + uz*uz)
	local eInt = eHydro - eKin
	local p = eInt / (self.gamma - 1)
	return rho, ux, uy, uz, Bx, By, Bz, p
end

-- using "An Upwind Differing Scheme for the Equations of Ideal Magnetohydrodynamics" by Brio & Wu
function MHD:calcInterfaceEigenBasis(sim,i,qL,qR)
	local gamma = self.gamma	
	local gammaMinusOne = gamma - 1

	local rhoL, uxL, uyL, uzL, BxL, ByL, BzL, pL = self:stateToPrims(unpack(qL))
	local rhoR, uxR, uyR, uzR, BxR, ByR, BzR, pR = self:stateToPrims(unpack(qR))

	local rho = .5*(rhoL + rhoR)
	local ux = .5*(uxL + uxR)
	local uy = .5*(uyL + uyR)
	local uz = .5*(uzL + uzR)
	local Bx = .5*(BxL + BxR)
	local By = .5*(ByL + ByR)
	local Bz = .5*(BzL + BzR)
	local p = .5*(pL + pR)

	local mu = self.mu
	local BSq = Bx*Bx + By*By + Bz*Bz
	local uSq = ux*ux + uy*uy + uz*uz
	local BdotU = Bx*ux + By*uy + Bz*uz
	local caxSq = Bx*Bx / (mu*rho)
	local cax = sqrt(caxSq)
	local aSq = gamma*p / rho
	local a = sqrt(aSq)
	local caSq = BSq / (mu*rho)
	local ca = sqrt(caSq)
	local aStarSq = aSq + caSq
	local discr = sqrt(aStarSq * aStarSq - 4 * caxSq * aSq)
	local cfSq = .5 * (aStarSq + discr)
	local csSq = .5 * (aStarSq - discr)
	local cf = sqrt(cfSq)
	local cs = sqrt(csSq)
	local sgnBx = Bx >= 0 and 1 or -1
	local sqrtRho = sqrt(rho)

	local ETotal = sim.qs[i][8]
	local PStar = p + .5 * BSq
	local H = (ETotal + PStar) / rho

	-- Brio & Wu
	sim.fluxMatrix[i][1][1] = 0
	sim.fluxMatrix[i][1][2] = 1
	sim.fluxMatrix[i][1][3] = 0
	sim.fluxMatrix[i][1][4] = 0
	sim.fluxMatrix[i][1][5] = 0
	sim.fluxMatrix[i][1][6] = 0
	sim.fluxMatrix[i][1][7] = 0
	sim.fluxMatrix[i][1][8] = 0
	sim.fluxMatrix[i][2][1] = (gamma-3)/2*ux^2 + (gamma-1)/2*(uy^2+uz^2)
	sim.fluxMatrix[i][2][2] = (3-gamma)*ux
	sim.fluxMatrix[i][2][3] = (1-gamma)*uy
	sim.fluxMatrix[i][2][4] = (1-gamma)*uz
	sim.fluxMatrix[i][2][5] = 0
	sim.fluxMatrix[i][2][6] = (2-gamma)*By
	sim.fluxMatrix[i][2][7] = (2-gamma)*Bz
	sim.fluxMatrix[i][2][8] = gamma-1
	sim.fluxMatrix[i][3][1] = -ux * uy
	sim.fluxMatrix[i][3][2] = uy
	sim.fluxMatrix[i][3][3] = ux
	sim.fluxMatrix[i][3][4] = 0
	sim.fluxMatrix[i][3][5] = 0
	sim.fluxMatrix[i][3][6] = -By
	sim.fluxMatrix[i][3][7] = 0
	sim.fluxMatrix[i][3][8] = 0
	sim.fluxMatrix[i][4][1] = -ux * uz
	sim.fluxMatrix[i][4][2] = uz
	sim.fluxMatrix[i][4][3] = 0
	sim.fluxMatrix[i][4][4] = ux
	sim.fluxMatrix[i][4][5] = 0
	sim.fluxMatrix[i][4][6] = 0
	sim.fluxMatrix[i][4][7] = -Bx
	sim.fluxMatrix[i][4][8] = 0
	sim.fluxMatrix[i][5][1] = 0
	sim.fluxMatrix[i][5][2] = 0
	sim.fluxMatrix[i][5][3] = 0
	sim.fluxMatrix[i][5][4] = 0
	sim.fluxMatrix[i][5][5] = 1
	sim.fluxMatrix[i][5][6] = 0
	sim.fluxMatrix[i][5][7] = 0
	sim.fluxMatrix[i][5][8] = 0
	sim.fluxMatrix[i][6][1] = (-By * ux + Bx * uy) / rho
	sim.fluxMatrix[i][6][2] = By / rho
	sim.fluxMatrix[i][6][3] = -Bx / rho
	sim.fluxMatrix[i][6][4] = 0
	sim.fluxMatrix[i][6][5] = 0
	sim.fluxMatrix[i][6][6] = ux
	sim.fluxMatrix[i][6][7] = 0
	sim.fluxMatrix[i][6][8] = 0
	sim.fluxMatrix[i][7][1] = (-Bz * ux + Bx * uz) / rho
	sim.fluxMatrix[i][7][2] = Bz / rho
	sim.fluxMatrix[i][7][3] = 0
	sim.fluxMatrix[i][7][4] = -Bx / rho
	sim.fluxMatrix[i][7][5] = 0
	sim.fluxMatrix[i][7][6] = 0
	sim.fluxMatrix[i][7][7] = ux
	sim.fluxMatrix[i][7][8] = 0
	sim.fluxMatrix[i][8][1] = -ux * (H * (gamma - 1) / 2 * uSq + Bx / rho * BdotU)
	sim.fluxMatrix[i][8][2] = H - Bx^2 / rho - (gamma - 1) * ux^2
	sim.fluxMatrix[i][8][3] = (1 - gamma) * ux * uy - Bx * By / rho
	sim.fluxMatrix[i][8][4] = (1 - gamma) * ux * uz - Bx * Bz / rho
	sim.fluxMatrix[i][8][5] = 0
	sim.fluxMatrix[i][8][6] = (2 - gamma) * By * ux - Bx * uy
	sim.fluxMatrix[i][8][7] = (2 - gamma) * Bz * ux - Bx * uz
	sim.fluxMatrix[i][8][8] = gamma * ux

	sim.eigenvalues[i][1] = ux - cf
	sim.eigenvalues[i][2] = ux - cax
	sim.eigenvalues[i][3] = ux - cs
	sim.eigenvalues[i][4] = ux
	sim.eigenvalues[i][5] = ux
	sim.eigenvalues[i][6] = ux + cs
	sim.eigenvalues[i][7] = ux + cax
	sim.eigenvalues[i][8] = ux + cf

	local bx = Bx / sqrtRho
	local by = By / sqrtRho
	local bz = Bz / sqrtRho
	local bxSq = bx * bx
	local bySq = by * by
	local bzSq = bz * bz

	local epsilon = 1e-20
	local BT = By^2 + Bz^2

	-- when by^2 + bz^2 ~= 0 ... (i.e. magnetic field tangent to normal exists)
	local alpha_f, alpha_s, beta_y, beta_z
	if BT > epsilon then
		alpha_f = sqrt((cfSq - bxSq) / (cfSq - csSq))
		--alpha_s = sqrt((bxSq - csSq) / (cfSq - csSq))/abs(bx) 
		alpha_s = sqrt((cfSq - aSq) / (cfSq - csSq))/cf 
		beta_y = By / sqrt(bySq + bzSq)
		beta_z = Bz / sqrt(bySq + bzSq)
		-- now TODO - multiply fast, slow, and Alfven eigenvectors by alpha_f, alpha_s, 1/sqrt(by^2 + bz^2)
	else
		alpha_f = 1
		alpha_s = 1
		beta_y = sqrt(.5)
		beta_z = sqrt(.5)
		-- and Cs^2 - bxSq = 0
	end

	-- h fast/slow, plus/minus
	local hfp = cfSq / (gamma - 1) + cf * ux - (By * uy + Bz * uz) / rho * Bx * cf / (cfSq - bxSq) + (gamma - 2) / (gamma - 1) * (cfSq - aSq)
	local hfm = cfSq / (gamma - 1) - cf * ux + (By * uy + Bz * uz) / rho * Bx * cf / (cfSq - bxSq) + (gamma - 2) / (gamma - 1) * (cfSq - aSq)
	local alpha_s_hsp = alpha_s * (csSq / (gamma - 1) + cs * ux + (gamma - 2) / (gamma - 1) * (csSq - aSq)) + (beta_y * uy + beta_z * uz) * a / cf * alpha_f * sgnBx		-- - alpha_s * (By * uy + Bz * uz) / rho * Bx * cs / (csSq - bxSq)
	local alpha_s_hsm = alpha_s * (csSq / (gamma - 1) - cs * ux + (gamma - 2) / (gamma - 1) * (csSq - aSq)) - (beta_y * uy + beta_z * uz) * a / cf * alpha_f * sgnBx		-- alpha_s * (By * uy + Bz * uz) / rho * Bx * cs / (csSq - bxSq)
	-- g plus/minus
	local gp_over_btSq = -(beta_z * uy - beta_y * uz) * sgnBx
	local gm_over_btSq = (beta_z * uy + beta_y * uz) * sgnBx

	--right eigenvectors
	--fast -
	sim.eigenvectors[i][1][1] = 1
	sim.eigenvectors[i][2][1] = ux - cf
	sim.eigenvectors[i][3][1] = uy + Bx * By * cf / (rho * (cfSq - bxSq))
	sim.eigenvectors[i][4][1] = uz + Bx * Bz * cf / (rho * (cfSq - bxSq))
	sim.eigenvectors[i][5][1] = 0
	sim.eigenvectors[i][6][1] = By * cfSq / (rho * (cfSq - bxSq))
	sim.eigenvectors[i][7][1] = Bz * cfSq / (rho * (cfSq - bxSq))
	sim.eigenvectors[i][8][1] = .5 * uSq + hfm
	for j=1,8 do
		sim.eigenvectors[i][j][1] = sim.eigenvectors[i][j][1] * alpha_f
	end
	--alfven -
	sim.eigenvectors[i][1][2] = 0
	sim.eigenvectors[i][2][2] = 0
	sim.eigenvectors[i][3][2] = beta_z * sgnBx
	sim.eigenvectors[i][4][2] = -beta_y * sgnBx
	sim.eigenvectors[i][5][2] = 0
	sim.eigenvectors[i][6][2] = beta_z / sqrtRho
	sim.eigenvectors[i][7][2] = -beta_y / sqrtRho
	sim.eigenvectors[i][8][2] = gm_over_btSq
	--slow -
	sim.eigenvectors[i][1][3] = alpha_s
	sim.eigenvectors[i][2][3] = alpha_s * (ux - cs)
	sim.eigenvectors[i][3][3] = alpha_s * uy + rho * (-a / cf * beta_y * alpha_f * sgnBx)	-- (alpha_s * by * bx * cs / rho / (csSq - bxSq))
	sim.eigenvectors[i][4][3] = alpha_s * uz + rho * (-a / cf * beta_z * alpha_f * sgnBx)	-- (alpha_s * bz * bx * cs / rho / (csSq - bxSq))
	sim.eigenvectors[i][5][3] = 0
	sim.eigenvectors[i][6][3] = -aSq / cfSq * beta_y / sqrtRho * alpha_f	-- alpha_s * by * csSq / sqrtRho / (csSq - bxSq)
	sim.eigenvectors[i][7][3] = -aSq / cfSq * beta_z / sqrtRho * alpha_f	-- alpha_s * bz * csSq / sqrtRho / (csSq - bxSq)
	sim.eigenvectors[i][8][3] = alpha_s * .5 * uSq + alpha_s_hsm
	--entropy
	sim.eigenvectors[i][1][4] = 1
	sim.eigenvectors[i][2][4] = ux
	sim.eigenvectors[i][3][4] = uy
	sim.eigenvectors[i][4][4] = uz
	sim.eigenvectors[i][5][4] = 0
	sim.eigenvectors[i][6][4] = 0
	sim.eigenvectors[i][7][4] = 0
	sim.eigenvectors[i][8][4] = .5 * uSq
	--zero
	sim.eigenvectors[i][1][5] = 0
	sim.eigenvectors[i][2][5] = 0
	sim.eigenvectors[i][3][5] = 0
	sim.eigenvectors[i][4][5] = 0
	sim.eigenvectors[i][5][5] = 1
	sim.eigenvectors[i][6][5] = 0
	sim.eigenvectors[i][7][5] = 0
	sim.eigenvectors[i][8][5] = 0
	--slow +
	sim.eigenvectors[i][1][6] = alpha_s
	sim.eigenvectors[i][2][6] = alpha_s * (ux + cs)
	sim.eigenvectors[i][3][6] = alpha_s * uy - rho * (-a / cf * beta_y * alpha_f * sgnBx)	-- (alpha_s * by * bx * cs / rho / (csSq - bxSq))
	sim.eigenvectors[i][4][6] = alpha_s * uz - rho * (-a / cf * beta_z * alpha_f * sgnBx)	-- (alpha_s * bz * bx * cs / rho / (csSq - bxSq))
	sim.eigenvectors[i][5][6] = 0
	sim.eigenvectors[i][6][6] = -aSq / cfSq * beta_y / sqrtRho * alpha_f	-- alpha_s * by * csSq / sqrtRho / (csSq - bxSq)
	sim.eigenvectors[i][7][6] = -aSq / cfSq * beta_z / sqrtRho * alpha_f	-- alpha_s * bz * csSq / sqrtRho / (csSq - bxSq)
	sim.eigenvectors[i][8][6] = alpha_s * .5 * uSq + alpha_s_hsp
	--alfven +
	sim.eigenvectors[i][1][7] = 0
	sim.eigenvectors[i][2][7] = 0
	sim.eigenvectors[i][3][7] = -beta_z * sgnBx
	sim.eigenvectors[i][4][7] = beta_y * sgnBx
	sim.eigenvectors[i][5][7] = 0
	sim.eigenvectors[i][6][7] = beta_z / sqrtRho
	sim.eigenvectors[i][7][7] = -beta_y / sqrtRho
	sim.eigenvectors[i][8][7] = gp_over_btSq
	--fast +
	sim.eigenvectors[i][1][8] = 1
	sim.eigenvectors[i][2][8] = ux + cf
	sim.eigenvectors[i][3][8] = uy - Bx * By * cf / (rho * (cfSq - bxSq))
	sim.eigenvectors[i][4][8] = uz - Bx * Bz * cf / (rho * (cfSq - bxSq))
	sim.eigenvectors[i][5][8] = 0
	sim.eigenvectors[i][6][8] = By * cfSq / (rho * (cfSq - bxSq))
	sim.eigenvectors[i][7][8] = Bz * cfSq / (rho * (cfSq - bxSq))
	sim.eigenvectors[i][8][8] = .5 * uSq + hfp
	for j=1,8 do
		sim.eigenvectors[i][j][8] = sim.eigenvectors[i][j][8] * alpha_f
	end

	assert(finite(rho), "rho is not finite")
	assert(finite(ux), "ux is not finite")
	assert(finite(uy), "uy is not finite")
	assert(finite(uz), "uz is not finite")
	assert(finite(Bx), "Bx is not finite")
	assert(finite(By), "By is not finite")
	assert(finite(Bz), "Bz is not finite")
	assert(finite(bx), "bx is not finite")
	assert(finite(by), "by is not finite")
	assert(finite(bz), "bz is not finite")
	assert(finite(cs), "cs is not finite")
	assert(finite(alpha_s_hsm), "alpha_s_hsm is not finite")
	assert(finite(alpha_s_hsp), "alpha_s_hsp is not finite")

	for j=1,8 do
		for k=1,8 do
			assert(finite(sim.eigenvectors[i][j][k]), "eigenvectors["..i.."]["..j.."]["..k.."] is not finite")
		end
	end

	-- use linear solver for eigenvector inverse
	sim.eigenvectorsInverse[i] = nil
end

return MHD

