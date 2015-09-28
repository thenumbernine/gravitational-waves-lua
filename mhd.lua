local class = require 'ext.class'
local Equation = require 'equation'

local function isnan(x) return x ~= x end
local function isinf(x) return x == math.huge or x == -math.huge end
local function finite(x) return not isnan(x) and not isinf(x) end

local MHD = class(Equation)

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
	local EMag = .5*(Bx*Bx + By*By + Bz*Bz) / function(self,i) return self.equation.mu end
	local EHydro = ETotal - EMag
	local EKin = .5 * rho * (ux^2 + uy^2 + uz^2)
	local EInt = EHydro - EKin
	local eInt = EInt / rho
	local gamma = function(self,i) return (self.equation.gamma - 1) end
	local p = eInt * gamma
	local pStar = p + EMag
	-- [[ full mhd
	MHD.graphInfos = table{
		{viewport={0/4, 2/4, 1/4, 1/4}, getter=function(self,i) return self.eigenbasisErrors and math.log(self.eigenbasisErrors[i]) end, name='eigenbasis error', color={1,0,0}, range={-30, 30}},
		{viewport={0/4, 3/4, 1/4, 1/4}, getter=function(self,i) return self.fluxMatrixErrors and math.log(self.fluxMatrixErrors[i]) end, name='reconstruction error', color={1,0,0}, range={-30, 30}},
		{viewport={0/4, 0/4, 1/4, 1/4}, getter=rho, name='rho', color={1,0,1}},
		{viewport={0/4, 1/4, 1/4, 1/4}, getter=p, name='p', color={1,0,1}},
		{viewport={1/4, 0/4, 1/4, 1/4}, getter=ux, name='ux', color={0,1,0}},
		{viewport={1/4, 1/4, 1/4, 1/4}, getter=uy, name='uy', color={0,1,0}},
		{viewport={1/4, 2/4, 1/4, 1/4}, getter=uz, name='uz', color={0,1,0}},
		{viewport={1/4, 3/4, 1/4, 1/4}, getter=pStar, name='pStar', color={0,1,0}},
		{viewport={2/4, 0/4, 1/4, 1/4}, getter=Bx, name='Bx', color={.5,.5,1}},
		{viewport={2/4, 1/4, 1/4, 1/4}, getter=By, name='By', color={.5,.5,1}},
		{viewport={2/4, 2/4, 1/4, 1/4}, getter=Bz, name='Bz', color={.5,.5,1}},
		{viewport={3/4, 0/4, 1/4, 1/4}, getter=ETotal, name='ETotal', color={1,1,0}},
		{viewport={3/4, 1/4, 1/4, 1/4}, getter=EKin, name='EKin', color={1,1,0}},
		{viewport={3/4, 2/4, 1/4, 1/4}, getter=EInt, name='EInt', color={1,1,0}},
		{viewport={3/4, 3/4, 1/4, 1/4}, getter=EHydro, name='EHydro', color={1,1,0}},
		{viewport={2/4, 3/4, 1/4, 1/4}, getter=EMag, name='EMag', color={1,1,0}},
	}
	--]]
	--[[ matching Euler
	MHD.graphInfos = table{
		{viewport={0/3, 0/2, 1/3, 1/2}, getter=function(self,i) return self.qs[i][1] end, name='rho', color={1,0,1}},
		{viewport={1/3, 0/2, 1/3, 1/2}, getter=function(self,i) return self.qs[i][2] / self.qs[i][1] end, name='u', color={0,1,0}},
		{viewport={2/3, 0/2, 1/3, 1/2}, getter=function(self,i) return self.qs[i][8] / self.qs[i][1] end, name='E', color={.5,.5,1}},
		{viewport={0/3, 1/2, 1/3, 1/2}, getter=function(self,i) return self.eigenbasisErrors and math.log(self.eigenbasisErrors[i]) end, name='log eigenbasis error', color={1,0,0}, range={-30, 30}},
		{viewport={1/3, 1/2, 1/3, 1/2}, getter=function(self,i) return self.fluxMatrixErrors and math.log(self.fluxMatrixErrors[i]) end, name='log reconstruction error', color={1,0,0}, range={-30, 30}},
	}
	--]]
end
MHD.graphInfoForNames = MHD.graphInfos:map(function(info,i)
	return info, info.name
end)

function MHD:initCell(sim,i)
	local gamma = self.gamma
	local x = sim.xs[i]
	local rho = x < 0 and 1 or .125
	local ux, uy, uz = 0, 0, 0
	-- [[ Brio & Wu
	local Bx = .75
	local By = x < 0 and 1 or -1
	local Bz = 0
	--]]
	--[[ some other tests
	local Bx, By, Bz = 0, sin(pi/2*x), 0
	--local Bx, By, Bz = 0, 1, 0	-- constant field works
	--]]
	--[[ Sod
	local Bx, By, Bz = 0, 0, 0	-- zero field works
	--]]
	local p = x < 0 and 1 or .1
	local eInt = p / (gamma-1)
	local EInt = rho * eInt
	local eKin = .5*(ux*ux + uy*uy + uz*uz)
	local EKin = rho * eKin
	local EMag = .5*(Bx*Bx + By*By + Bz*Bz) / self.mu
	local ETotal = EInt + EKin + EMag 
	local mx, my, mz = rho * ux, rho * uy, rho * uz
	return {rho, mx, my, mz, Bx, By, Bz, ETotal}
end

function MHD:stateToPrims(rho, mx, my, mz, Bx, By, Bz, ETotal)
	local ux, uy, uz = mx / rho, my / rho, mz / rho
	local EMag = .5*(Bx*Bx + By*By + Bz*Bz) / self.mu
	assert(EMag < ETotal, "magnetic energy ("..EMag..") >= total energy ("..ETotal..")")
	local EHydro = ETotal - EMag
	assert(EHydro > 0, "densitized hydro energy is negative: "..EHydro)
	local eKin = .5*(ux*ux + uy*uy + uz*uz)
	local EKin = rho * eKin
	local EInt = EHydro - EKin
	assert(EInt > 0, "densitized internal energy is negative: "..EInt)
	local p = EInt*(self.gamma-1)
	assert(p > 0, "static pressure is negative: "..p)
	local pStar = p + EMag
	assert(pStar > 0, "full pressure is negative: "..pStar)
	return rho, ux, uy, uz, Bx, By, Bz, pStar
end

-- using "An Upwind Differing Scheme for the Equations of Ideal Magnetohydrodynamics" by Brio & Wu
function MHD:calcInterfaceEigenBasis(sim,i,qL,qR)
	local gamma = self.gamma	
	local gammaMinusOne = gamma - 1

	local rhoL, uxL, uyL, uzL, BxL, ByL, BzL, pStarL = self:stateToPrims(unpack(qL))
	local rhoR, uxR, uyR, uzR, BxR, ByR, BzR, pStarR = self:stateToPrims(unpack(qR))

	-- [[ average
	local rho = .5*(rhoL + rhoR)
	local ux = .5*(uxL + uxR)
	local uy = .5*(uyL + uyR)
	local uz = .5*(uzL + uzR)
	local Bx = .5*(BxL + BxR)
	local By = .5*(ByL + ByR)
	local Bz = .5*(BzL + BzR)
	local pStar = .5*(pStarL + pStarR)
	--]]
	--[[ Roe averaging ... seems to be more stable when dealing with negative pressure (or magnetic energy getting too high)
	local sqrtRhoL = sqrt(rhoL)
	local sqrtRhoR = sqrt(rhoR)
	local denom = 1 / (sqrtRhoL + sqrtRhoR)
	local rho = sqrtRhoL * sqrtRhoR
	local ux = (sqrtRhoL*uxL + sqrtRhoR*uxR)*denom
	local uy = (sqrtRhoL*uyL + sqrtRhoR*uyR)*denom
	local uz = (sqrtRhoL*uzL + sqrtRhoR*uzR)*denom
	local Bx = (sqrtRhoL*BxL + sqrtRhoR*BxR)*denom
	local By = (sqrtRhoL*ByL + sqrtRhoR*ByR)*denom
	local Bz = (sqrtRhoL*BzL + sqrtRhoR*BzR)*denom
	local pStar = (sqrtRhoL*pStarL + sqrtRhoR*pStarR)*denom
	--]]
	
	local ETotal = sim.qs[i][8]
	local H = (ETotal + pStar) / rho

	local mu = self.mu
	local BSq = Bx*Bx + By*By + Bz*Bz
	local uSq = ux*ux + uy*uy + uz*uz
	local BdotU = Bx*ux + By*uy + Bz*uz
	local p = pStar - .5 * BSq
	local caxSq = Bx*Bx / (mu*rho)
	local cax = sqrt(caxSq)	-- positive
	--local aSq = gamma * p / rho
	local aSq = (gamma - 1) * (H - uSq - BSq / rho) - (gamma - 2) * ((ByL-ByR)^2 + (BzL-BzR)^2)/(2*(math.sqrt(rhoL) + math.sqrt(rhoR)))
	local a = sqrt(aSq)
	local caSq = BSq / (mu*rho)
	local ca = sqrt(caSq)
	local aStarSq = (gamma * p + BSq) / rho		-- aSq + caSq
	local discr = sqrt(aStarSq * aStarSq - 4 * caxSq * aSq)
	local cfSq = .5 * (aStarSq + discr)
	local csSq = .5 * (aStarSq - discr)
	local cf = sqrt(cfSq)
	local cs = sqrt(csSq)
	local sgnBx = Bx >= 0 and 1 or -1
	local sqrtRho = sqrt(rho)

	-- Brio & Wu
	sim.fluxMatrix[i][1][1] = 0
	sim.fluxMatrix[i][1][2] = 1
	sim.fluxMatrix[i][1][3] = 0
	sim.fluxMatrix[i][1][4] = 0
	sim.fluxMatrix[i][1][5] = 0
	sim.fluxMatrix[i][1][6] = 0
	sim.fluxMatrix[i][1][7] = 0
	sim.fluxMatrix[i][1][8] = 0
	sim.fluxMatrix[i][2][1] = (gamma-3)/2*ux^2 + (gamma-1)/2*(uy^2+uz^2) - (gamma-2)*((ByL-ByR)^2 + (BzL-BzR)^2)/(2*(math.sqrt(rhoL) + math.sqrt(rhoR)))
	sim.fluxMatrix[i][2][2] = (3-gamma)*ux
	sim.fluxMatrix[i][2][3] = (1-gamma)*uy
	sim.fluxMatrix[i][2][4] = (1-gamma)*uz
	sim.fluxMatrix[i][2][5] = 0
	sim.fluxMatrix[i][2][6] = (2-gamma)*By*(rhoL+rhoR)/(2*rho)
	sim.fluxMatrix[i][2][7] = (2-gamma)*Bz*(rhoL+rhoR)/(2*rho)
	sim.fluxMatrix[i][2][8] = gamma-1
	sim.fluxMatrix[i][3][1] = -ux * uy
	sim.fluxMatrix[i][3][2] = uy
	sim.fluxMatrix[i][3][3] = ux
	sim.fluxMatrix[i][3][4] = 0
	sim.fluxMatrix[i][3][5] = 0
	sim.fluxMatrix[i][3][6] = -Bx
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
	sim.fluxMatrix[i][8][1] = -ux * H + (gamma - 1) * ux * uSq / 2 + Bx * BdotU / rho - ux * (gamma - 2) * ((ByL-ByR)^2 + (BzL-BzR)^2)/(2*(math.sqrt(rhoL) + math.sqrt(rhoR)))
	sim.fluxMatrix[i][8][2] = H - Bx^2 / rho - (gamma - 1) * ux^2
	sim.fluxMatrix[i][8][3] = (1 - gamma) * ux * uy - Bx * By / rho
	sim.fluxMatrix[i][8][4] = (1 - gamma) * ux * uz - Bx * Bz / rho
	sim.fluxMatrix[i][8][5] = 0
	sim.fluxMatrix[i][8][6] = -Bx * uy - (gamma - 2) * By * ux * (rhoL + rhoR) / (2 * rho)
	sim.fluxMatrix[i][8][7] = -Bx * uz - (gamma - 2) * Bz * ux * (rhoL + rhoR) / (2 * rho)
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

	local epsilon = 1e-12
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
--	for j=1,8 do
--		sim.eigenvectors[i][j][1] = sim.eigenvectors[i][j][1] * alpha_f
--	end
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
--	for j=1,8 do
--		sim.eigenvectors[i][j][8] = sim.eigenvectors[i][j][8] * alpha_f
--	end

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
	assert(p >= 0, "pressure is negative: "..tostring(p))
	assert(finite(aStarSq), "aStarSq is not finite")
	assert(finite(caxSq), "caxSq is not finite")
	assert(finite(aSq), "aSq is not finite")
	assert(aStarSq * aStarSq - 4 * caxSq * aSq >= 0, "discr is imaginary")
	assert(finite(discr), "discr is not finite: "..tostring(discr))
	assert(aStarSq >= discr, "aStarSq ("..tostring(aStarSq)..") is not >= discr ("..tostring(discr)..")")
	assert(finite(csSq) and csSq >= 0, "csSq is not finite and positive: "..tostring(csSq))
	assert(finite(cs), "cs is not finite: "..tostring(cs))
	assert(finite(alpha_s_hsm), "alpha_s_hsm is not finite")
	assert(finite(alpha_s_hsp), "alpha_s_hsp is not finite")

	for j=1,8 do
		for k=1,8 do
			assert(finite(sim.eigenvectors[i][j][k]), "eigenvectors["..i.."]["..j.."]["..k.."] is not finite")
		end
	end

	-- use linear solver for eigenvector inverse
--	sim.eigenvectorsInverse[i] = nil
end




-- A x = b
-- x = A^-1 b
-- returns x
-- uses Gauss-Jordan method
local function linearSolveGaussJordan(A, b)
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
print(unpack(linearSolveGaussJordan(
	-- store row-major so Lua indexing matches math indexing
	{
		{3,0,0},
		{2,1,0},
		{1,0,1},
	},
		{1,2,3}
)))
os.exit()
--]]

-- solve for x in Ax=b using Gauss-Seidel iterative method
local function linearSolveGaussSeidel(A, b)
	local result = {unpack(b)}
	local maxiter = 20
	local n = #A
	for iter=1,maxiter do
		local lastResult = {unpack(result)}
		for i=1,n do
			local sum = 0
			for j=1,n do
				if j ~= i then
					sum = sum + A[i][j] * result[j]
				end
			end
			-- what if A[i][i] == 0 ?
			result[i] = (b[i] - sum) / A[i][i]
		end
		local err = 0
		for i=1,n do
			local delta = result[i] - lastResult[i]
			err = err + .5 * delta * delta
		end
		print('iter',i,'error',err)
	end
	return result
end

local function linearSolveConjGrad(A, b)
	local n = #A
	local function dot(a,b) return range(n):map(function(i) return a[i] * b[i] end):sum() end
	local function norm(a) return dot(a,a) end
	local maxiter = 100
	local epsilon = 1e-10
	local x = {unpack(b)}
	local r = range(n):map(function(i) return b[i] - dot(A[i], x) end)
	local r2 = norm(r)
	if r2 < epsilon then return x end
	local p = {unpack(r)}
	for iter=1,maxiter do
		local Ap = range(n):map(function(i) return dot(A[i], p) end)
		local alpha = r2 / dot(p, Ap)
		x = range(n):map(function(i) return x[i] + alpha * p[i] end)
		local nr = range(n):map(function(i) return r[i] - alpha * Ap[i] end)
		local nr2 = norm(nr)
		local beta = nr2 / r2
		if nr2 < epsilon then break end
		r = nr
		r2 = nr2
		p = range(n):map(function(i) return r[i] + beta * p[i] end)
	end
	return x
end

-- TODO use linearsolvers.lua
local linearSolve = linearSolveGaussJordan
--local linearSolve = linearSolveGaussSeidel
--local linearSolve = linearSolveConjGrad

--[==[
function MHD:eigenfields(sim, i, v)

	-- eigenvector transform:
	-- y = R*v; return y
	-- inverse transform:
	-- x = R^-1 * v; return x;
	-- v = R*x
	-- solve linear system above, use right-eigenvectors as matrix and v as solution
	return linearSolve(sim.eigenvectors[i], v)

--[=[ me trying to decypher Brio & Wu's paper

	local qL = sim.qs[i-1]
	local qR = sim.qs[i]

-- begin block similar to calcInterfaceEigenBasis
	local gamma = sim.gamma	
	local gammaMinusOne = gamma - 1
	
	local rhoL, uxL, uyL, uzL, BxL, ByL, BzL, pL = sim.equation:stateToPrims(unpack(qL))
	local rhoR, uxR, uyR, uzR, BxR, ByR, BzR, pR = sim.equation:stateToPrims(unpack(qR))

	local rho = .5*(rhoL + rhoR)
	local ux = .5*(uxL + uxR)
	local uy = .5*(uyL + uyR)
	local uz = .5*(uzL + uzR)
	local Bx = .5*(BxL + BxR)
	local By = .5*(ByL + ByR)
	local Bz = .5*(BzL + BzR)
	local p = .5*(pL + pR)

	local mu = sim.equation.mu
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
--]==]


return MHD

