local class = require 'ext.class'
local Equation = require 'equation'

local function isnan(x) return x ~= x end
local function isinf(x) return x == math.huge or x == -math.huge end
local function isfinite(x) return not isnan(x) and not isinf(x) end

local MHD = class(Equation)

MHD.numStates = 8
MHD.gamma = 5/3	
MHD.mu = 1

do
	local rho = function(self,i) return self.qs[i][1] end
	local vx = function(self,i) return self.qs[i][2]/rho(self,i) end
	local vy = function(self,i) return self.qs[i][3]/rho(self,i) end
	local vz = function(self,i) return self.qs[i][4]/rho(self,i) end
	local bx = function(self,i) return self.qs[i][5] end
	local by = function(self,i) return self.qs[i][6] end
	local bz = function(self,i) return self.qs[i][7] end
	local ETotal = function(self,i) return self.qs[i][8] end
	local EMag = .5*(bx*bx + by*by + bz*bz) / function(self,i) return self.equation.mu end
	local EHydro = ETotal - EMag
	local EKin = .5 * rho * (vx^2 + vy^2 + vz^2)
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
		{viewport={1/4, 0/4, 1/4, 1/4}, getter=vx, name='vx', color={0,1,0}},
		{viewport={1/4, 1/4, 1/4, 1/4}, getter=vy, name='vy', color={0,1,0}},
		{viewport={1/4, 2/4, 1/4, 1/4}, getter=vz, name='vz', color={0,1,0}},
		{viewport={1/4, 3/4, 1/4, 1/4}, getter=pStar, name='pStar', color={0,1,0}},
		{viewport={2/4, 0/4, 1/4, 1/4}, getter=bx, name='bx', color={.5,.5,1}},
		{viewport={2/4, 1/4, 1/4, 1/4}, getter=by, name='by', color={.5,.5,1}},
		{viewport={2/4, 2/4, 1/4, 1/4}, getter=bz, name='bz', color={.5,.5,1}},
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
	local vx, vy, vz = 0, 0, 0
	--[[ Brio & Wu
	local bx = 1
	local by = x < 0 and 1 or -1
	local bz = 0
	--]]
	--[[ some other tests
	local bx, by, bz = 0, sin(pi/2*x), 0
	--local bx, by, bz = 0, 1, 0	-- constant field works
	--]]
	-- [[ Sod
	local bx, by, bz = 0, 0, 0	-- zero field works
	--]]
	local p = x < 0 and 1 or .1
	local eInt = p / (gamma-1)
	local EInt = rho * eInt
	local eKin = .5*(vx*vx + vy*vy + vz*vz)
	local EKin = rho * eKin
	local bSq = bx*bx + by*by + bz*bz
	local EMag = .5*bSq
	local ETotal = EInt + EKin + EMag 
	local mx, my, mz = rho * vx, rho * vy, rho * vz
	return {rho, mx, my, mz, bx, by, bz, ETotal}
end

function MHD:stateToPrims(rho, mx, my, mz, bx, by, bz, ETotal)
	local gamma = self.gamma
	local vx, vy, vz = mx / rho, my / rho, mz / rho
	local vSq = vx*vx + vy*vy + vz*vz
	local bSq = bx*bx + by*by + bz*bz
	local EKin = 1/2 * rho * vSq
	local EMag = bSq/2
	local P = (ETotal - EKin - EMag) * (gamma - 1)		
	local PStar = P + EMag								
	local H = (ETotal + PStar) / rho					
	return rho, vx, vy, vz, bx, by, bz, H
end

function MHD:calcInterfaceEigenBasis(sim,i,qL,qR)
	local gamma = self.gamma	
	local gammaMinusOne = gamma - 1

	local rhoL, vxL, vyL, vzL, bxL, byL, bzL, HL = self:stateToPrims(unpack(qL))
	local rhoR, vxR, vyR, vzR, bxR, byR, bzR, HR = self:stateToPrims(unpack(qR))

	-- eqn 56 - averaging
	local sqrtRhoL = sqrt(rhoL)
	local sqrtRhoR = sqrt(rhoR)
	local invDenom = 1 / (sqrtRhoL + sqrtRhoR)
	local rho = sqrtRhoL * sqrtRhoR
	local vx = (sqrtRhoL*vxL + sqrtRhoR*vxR)*invDenom
	local vy = (sqrtRhoL*vyL + sqrtRhoR*vyR)*invDenom
	local vz = (sqrtRhoL*vzL + sqrtRhoR*vzR)*invDenom
	local bx = (sqrtRhoL*bxL + sqrtRhoR*bxR)*invDenom
	local by = (sqrtRhoL*byL + sqrtRhoR*byR)*invDenom
	local bz = (sqrtRhoL*bzL + sqrtRhoR*bzR)*invDenom
	local H = (sqrtRhoL*HL + sqrtRhoR*HR)*invDenom

	-- eqn 56 says to average H
	-- but B8 says H = (E + P + b^2/2)/rho

	local vSq = vx*vx + vy*vy + vz*vz

	local sqrtRho = sqrt(rho)

	-- TODO divide by mu instead?
	local oneOverSqrtMu = 1 / sqrt(4 * pi)
	local bSq = bx*bx + by*by + bz*bz
	local bDotV = bx*vx + by*vy + bz*vz

	local gammaPrime = gamma - 1
	local X = ((byL - byR)^2 + (bzL - bzR)^2) / (2 * (sqrt(rhoL) + sqrt(rhoR)))
	local XPrime = (gamma - 2) * X
	local Y = (rhoL + rhoR) / (2 * rho)
	local YPrime = (gamma - 2) * Y

	local fluxMatrix = sim.fluxMatrix[i]
	fluxMatrix[1][1] = 0
	fluxMatrix[1][2] = 1
	fluxMatrix[1][3] = 0
	fluxMatrix[1][4] = 0
	fluxMatrix[1][5] = 0
	fluxMatrix[1][6] = 0
	fluxMatrix[1][7] = 0
	fluxMatrix[1][8] = 0
	fluxMatrix[2][1] = -vx^2 + gammaPrime * vSq / 2 - XPrime
	fluxMatrix[2][2] = -(gamma-3)*vx
	fluxMatrix[2][3] = -gammaPrime*vy
	fluxMatrix[2][4] = -gammaPrime*vz
	fluxMatrix[2][5] = 0
	fluxMatrix[2][6] = -by*YPrime
	fluxMatrix[2][7] = -bz*YPrime
	fluxMatrix[2][8] = gammaPrime
	fluxMatrix[3][1] = -vx * vy
	fluxMatrix[3][2] = vy
	fluxMatrix[3][3] = vx
	fluxMatrix[3][4] = 0
	fluxMatrix[3][5] = 0
	fluxMatrix[3][6] = -bx
	fluxMatrix[3][7] = 0
	fluxMatrix[3][8] = 0
	fluxMatrix[4][1] = -vx * vz
	fluxMatrix[4][2] = vz
	fluxMatrix[4][3] = 0
	fluxMatrix[4][4] = vx
	fluxMatrix[4][5] = 0
	fluxMatrix[4][6] = 0
	fluxMatrix[4][7] = -bx
	fluxMatrix[4][8] = 0
	fluxMatrix[5][1] = 0
	fluxMatrix[5][2] = 0
	fluxMatrix[5][3] = 0
	fluxMatrix[5][4] = 0
	fluxMatrix[5][5] = 1
	fluxMatrix[5][6] = 0
	fluxMatrix[5][7] = 0
	fluxMatrix[5][8] = 0
	fluxMatrix[6][1] = (bx * vy - by * vx) / rho
	fluxMatrix[6][2] = by / rho
	fluxMatrix[6][3] = -bx / rho
	fluxMatrix[6][4] = 0
	fluxMatrix[6][5] = 0
	fluxMatrix[6][6] = vx
	fluxMatrix[6][7] = 0
	fluxMatrix[6][8] = 0
	fluxMatrix[7][1] = (bx * vz - bz * vx) / rho
	fluxMatrix[7][2] = bz / rho
	fluxMatrix[7][3] = 0
	fluxMatrix[7][4] = -bx / rho
	fluxMatrix[7][5] = 0
	fluxMatrix[7][6] = 0
	fluxMatrix[7][7] = vx
	fluxMatrix[7][8] = 0
	fluxMatrix[8][1] = -vx * H + gammaPrime * vx * vSq / 2 + bx * bDotV / rho - vx * XPrime
	fluxMatrix[8][2] = -gammaPrime * vx * vx + H - bx * bx / rho
	fluxMatrix[8][3] = -gammaPrime * vx * vy - bx * by / rho
	fluxMatrix[8][4] = -gammaPrime * vx * vz - bx * bz / rho
	fluxMatrix[8][5] = 0
	fluxMatrix[8][6] = -(bx * vy - by * vx * YPrime)
	fluxMatrix[8][7] = -(bx * vz - bz * vx * YPrime)
	fluxMatrix[8][8] = gamma * vx

	local aTildeSq = gammaPrime * (H - vSq/2 - bSq/rho) - XPrime
	local CAxSq = bx*bx / rho
	local bPerpSq = by*by + bz*bz
	local bStarPerpSq = (gammaPrime - YPrime) * bPerpSq
	local CATildeSq = CAxSq + bStarPerpSq / rho
	local CAx = sqrt(CAxSq)
	local CfSq = .5 * ((aTildeSq + CATildeSq) + sqrt((aTildeSq + CATildeSq)^2 - 4 * aTildeSq * CAxSq))
	local CsSq = .5 * ((aTildeSq + CATildeSq) - sqrt((aTildeSq + CATildeSq)^2 - 4 * aTildeSq * CAxSq))
	local Cf = sqrt(CfSq)
	local Cs = sqrt(CsSq)

	local eigenvalues = sim.eigenvalues[i]
	eigenvalues[1] = vx - Cf
	eigenvalues[2] = vx - CAx
	eigenvalues[3] = vx - Cs
	eigenvalues[4] = vx
	eigenvalues[5] = 0
	eigenvalues[6] = vx + Cs
	eigenvalues[7] = vx + CAx
	eigenvalues[8] = vx + Cf

	-- to prevent divide-by-zero errors
	local epsilon = 1e-20

	-- eqn A13-A17, replace 'a' with 'aTilde'
	local aTilde = sqrt(aTildeSq)
	local S = bx >= 0 and 1 or -1
	local alpha_f, alpha_s
	if CaSq == aTildeSq and CAxSq == aTildeSq then
		alpha_f = 1
		alpha_s = 0
	else
		alpha_f = (aTildeSq - CsSq) / (CfSq - CsSq)
		alpha_s = (CfSq - aTildeSq) / (CfSq - CsSq)
	end
	local Cff = Cf * alpha_f
	local Css = Cs * alpha_s
	local Qf = Cf * alpha_f * S
	local Qs = Cs * alpha_s * S
	local Af = aTilde * alpha_f * sqrtRho
	local As = aTilde * alpha_s * sqrtRho
	local bPerp = sqrt(bPerpSq)
	local beta_y, beta_z
	if bPerp < epsilon then
		beta_y = 1
		beta_z = 0
	else
		beta_y = by / bPerp
		beta_z = bz / bPerp
	end

	-- V[xyz][fs] = v[xyz] * alpha_[fs]
	local Vxf = vx * alpha_f
	local Vyf = vy * alpha_f
	local Vzf = vz * alpha_f
	local Vxs = vx * alpha_s
	local Vys = vy * alpha_s
	local Vzs = vz * alpha_s

	local bStarPerp = sqrt(bStarPerpSq)
	local betaStar_y, betaStar_z
	if bStarPerp < epsilon then
		betaStar_y = 1
		betaStar_z = 0
	else
		betaStar_y = by / bStarPerp
		betaStar_z = bz / bStarPerp
	end
	local betaStarPerpSq = betaStar_y * betaStar_y + betaStar_z * betaStar_z

	local HPrime = H - bSq / rho		-- (densitized) total gas enthalpy (internal energy + kinetic energy + pressure)

	local eigenvectors = sim.eigenvectors[i]
	--right eigenvectors
	--fast -
	eigenvectors[1][1] = alpha_f
	eigenvectors[2][1] = Vxf - Cff
	eigenvectors[3][1] = Vyf + Qs * betaStar_y
	eigenvectors[4][1] = Vzf + Qs * betaStar_z
	eigenvectors[5][1] = 0
	eigenvectors[6][1] = As * betaStar_y / rho
	eigenvectors[7][1] = As * betaStar_z / rho
	eigenvectors[8][1] = alpha_f * (HPrime - vx * Cf) + Qs * (vy * betaStar_y + vz * betaStar_z) + As * bStarPerp * betaStarPerpSq / rho
	--alfven -
	eigenvectors[1][2] = 0
	eigenvectors[2][2] = 0
	eigenvectors[3][2] = -beta_z
	eigenvectors[4][2] = beta_y
	eigenvectors[5][2] = 0
	eigenvectors[6][2] = -beta_z * S / sqrtRho
	eigenvectors[7][2] = beta_y * S / sqrtRho
	eigenvectors[8][2] = -(vy * beta_z - vz * beta_y)
	--slow -
	eigenvectors[1][3] = alpha_s
	eigenvectors[2][3] = Vxs - Css
	eigenvectors[3][3] = Vys - Qf * betaStar_y
	eigenvectors[4][3] = Vzs - Qf * betaStar_z
	eigenvectors[5][3] = 0
	eigenvectors[6][3] = -Af * betaStar_y / rho
	eigenvectors[7][3] = -Af * betaStar_z / rho
	eigenvectors[8][3] = alpha_s * (HPrime - vx * Cs) - Qf * (vy * betaStar_y + vz * betaStar_z) - Af * bStarPerp * betaStarPerpSq / rho
	--entropy
	eigenvectors[1][4] = 1
	eigenvectors[2][4] = vx
	eigenvectors[3][4] = vy
	eigenvectors[4][4] = vz
	eigenvectors[5][4] = 0
	eigenvectors[6][4] = 0
	eigenvectors[7][4] = 0
	eigenvectors[8][4] = vSq/2 + XPrime / gammaPrime
	--zero
	eigenvectors[1][5] = 0
	eigenvectors[2][5] = 0
	eigenvectors[3][5] = 0
	eigenvectors[4][5] = 0
	eigenvectors[5][5] = 1
	eigenvectors[6][5] = 0
	eigenvectors[7][5] = 0
	eigenvectors[8][5] = 0
	--slow +
	eigenvectors[1][6] = alpha_s
	eigenvectors[2][6] = Vxs + Css
	eigenvectors[3][6] = Vys + Qf * betaStar_y
	eigenvectors[4][6] = Vzs + Qf * betaStar_z
	eigenvectors[5][6] = 0
	eigenvectors[6][6] = -Af * betaStar_y / rho
	eigenvectors[7][6] = -Af * betaStar_z / rho
	eigenvectors[8][6] = alpha_s * (HPrime + vx * Cs) + Qf * (vy * betaStar_y + vz * betaStar_z) - Af * bStarPerp * betaStarPerpSq / rho
	--alfven +
	eigenvectors[1][7] = 0
	eigenvectors[2][7] = 0
	eigenvectors[3][7] = beta_z
	eigenvectors[4][7] = -beta_y
	eigenvectors[5][7] = 0
	eigenvectors[6][7] = -beta_z * S / sqrtRho
	eigenvectors[7][7] = beta_y * S / sqrtRho
	eigenvectors[8][7] = vy * beta_z - vz * beta_y
	--fast +
	eigenvectors[1][8] = alpha_f
	eigenvectors[2][8] = Vxf + Cff
	eigenvectors[3][8] = Vyf - Qs * betaStar_y
	eigenvectors[4][8] = Vzf - Qs * betaStar_z
	eigenvectors[5][8] = 0
	eigenvectors[6][8] = As * betaStar_y / rho
	eigenvectors[7][8] = As * betaStar_z / rho
	eigenvectors[8][8] = alpha_f * (HPrime + vx * Cf) - Qs * (vy * betaStar_y + vz * betaStar_z) + As * bStarPerp * betaStarPerpSq / rho

	local QStar_y, QStar_z
	if betaStarPerpSq < epsilon then
		QStar_y = 1
		QStar_z = 0
	else
		QStar_y = betaStar_y / betaStarPerpSq
		QStar_z = betaStar_z / betaStarPerpSq
	end

	--local aSq = gamma * P / rho	-- adiabatic sound speed
	local aSq = gammaPrime * (HPrime - .5 * vSq)
	local hatScalar = 1 / (2 * aSq)
	local QHat_f = Qf * hatScalar
	local QHat_s = Qs * hatScalar
	local AHat_f = Af * hatScalar	
	local AHat_s = As * hatScalar	
	local CHat_ff = Cff * hatScalar
	local CHat_ss = Css * hatScalar
	local XHatPrime = XPrime * hatScalar
	
	local barScalar = gammaPrime * hatScalar
	local alphaBar_f = alpha_f * barScalar
	local alphaBar_s = alpha_s * barScalar
	local VBar_xf = Vxf * barScalar
	local VBar_yf = Vyf * barScalar
	local VBar_zf = Vzf * barScalar
	local VBar_xs = Vxs * barScalar
	local VBar_ys = Vys * barScalar
	local VBar_zs = Vzs * barScalar
	local vBar_x = vx * barScalar
	local vBar_y = vy * barScalar
	local vBar_z = vz * barScalar
	local bBar_y = by * barScalar
	local bBar_z = bz * barScalar
	
	local vBarSq = vSq * barScalar * barScalar

	-- use linear solver for eigenvector inverse
	local eigenvectorsInverse = sim.eigenvectorsInverse[i]
	--fast -
	eigenvectorsInverse[1][1] = alphaBar_f * (vSq - HPrime) + CHat_ff * (Cf + vx) - QHat_s * (vy * QStar_y + vz * QStar_z) - AHat_s * bPerp / rho
	eigenvectorsInverse[1][2] = -VBar_xf - CHat_ff
	eigenvectorsInverse[1][3] = -VBar_yf + QHat_s * QStar_y
	eigenvectorsInverse[1][4] = -VBar_zf + QHat_s * QStar_z
	eigenvectorsInverse[1][5] = 0
	eigenvectorsInverse[1][6] = AHat_s * QStar_y - alphaBar_f * by
	eigenvectorsInverse[1][7] = AHat_s * QStar_z - alphaBar_f * bz
	eigenvectorsInverse[1][8] = alphaBar_f
	--alfven -
	eigenvectorsInverse[2][1] = (vy * beta_z - vz * beta_y) / 2
	eigenvectorsInverse[2][2] = 0
	eigenvectorsInverse[2][3] = -beta_z / 2
	eigenvectorsInverse[2][4] = beta_y / 2
	eigenvectorsInverse[2][5] = 0
	eigenvectorsInverse[2][6] = -beta_z * S * sqrtRho / 2
	eigenvectorsInverse[2][7] = beta_y * S * sqrtRho / 2
	eigenvectorsInverse[2][8] = 0
	--slow - 
	eigenvectorsInverse[3][1] = alphaBar_s * (vSq - HPrime) + CHat_ss * (Cs + vx) + QHat_f * (vy * QStar_y + vz * QStar_z) + AHat_f * bPerp / rho
	eigenvectorsInverse[3][2] = -VBar_xs - CHat_ss
	eigenvectorsInverse[3][3] = -VBar_ys - QHat_f * QStar_y
	eigenvectorsInverse[3][4] = -VBar_zs - QHat_f * QStar_z
	eigenvectorsInverse[3][5] = 0
	eigenvectorsInverse[3][6] = -AHat_f * QStar_y - alphaBar_s * by
	eigenvectorsInverse[3][7] = -AHat_f * QStar_z - alphaBar_s * bz
	eigenvectorsInverse[3][8] = alphaBar_s
	--entropy
	eigenvectorsInverse[4][1] = 1 - vBarSq + 2 * XHatPrime
	eigenvectorsInverse[4][2] = 2 * vBar_x
	eigenvectorsInverse[4][3] = 2 * vBar_y
	eigenvectorsInverse[4][4] = 2 * vBar_z
	eigenvectorsInverse[4][5] = 0
	eigenvectorsInverse[4][6] = 2 * bBar_y
	eigenvectorsInverse[4][7] = 2 * bBar_z
	eigenvectorsInverse[4][8] = -gammaPrime / aSq
	--zero
	eigenvectorsInverse[5][1] = 0
	eigenvectorsInverse[5][2] = 0
	eigenvectorsInverse[5][3] = 0
	eigenvectorsInverse[5][4] = 0
	eigenvectorsInverse[5][5] = 1
	eigenvectorsInverse[5][6] = 0
	eigenvectorsInverse[5][7] = 0
	eigenvectorsInverse[5][8] = 0
	--slow +
	eigenvectorsInverse[6][1] = alphaBar_s * (vSq - HPrime) + CHat_ss * (Cs - vx) - QHat_f * (vy * QStar_y + vz * QStar_z) + AHat_f * bPerp / rho
	eigenvectorsInverse[6][2] = -VBar_xs + CHat_ss
	eigenvectorsInverse[6][3] = -VBar_ys + QHat_f * QStar_y
	eigenvectorsInverse[6][4] = -VBar_zs + QHat_f * QStar_z
	eigenvectorsInverse[6][5] = 0
	eigenvectorsInverse[6][6] = -AHat_f * QStar_y - alphaBar_s * by
	eigenvectorsInverse[6][7] = -AHat_f * QStar_z - alphaBar_s * bz
	eigenvectorsInverse[6][8] = alphaBar_s
	--alfven +
	eigenvectorsInverse[7][1] = -(vy * beta_z - vz * beta_y) / 2
	eigenvectorsInverse[7][2] = 0
	eigenvectorsInverse[7][3] = beta_z / 2
	eigenvectorsInverse[7][4] = -beta_y / 2
	eigenvectorsInverse[7][5] = 0
	eigenvectorsInverse[7][6] = -beta_z * S * sqrtRho / 2
	eigenvectorsInverse[7][7] = beta_y * S * sqrtRho / 2
	eigenvectorsInverse[7][8] = 0
	--fast +
	eigenvectorsInverse[8][1] = alphaBar_f * (vSq - HPrime) + CHat_ff * (Cf - vx) + QHat_s * (vy * QStar_y + vz * QStar_z) - AHat_s * bPerp / rho
	eigenvectorsInverse[8][2] = -VBar_xf + CHat_ff
	eigenvectorsInverse[8][3] = -VBar_yf - QHat_s * QStar_y
	eigenvectorsInverse[8][4] = -VBar_zf - QHat_s * QStar_z
	eigenvectorsInverse[8][5] = 0
	eigenvectorsInverse[8][6] = AHat_s * QStar_y - alphaBar_f * by
	eigenvectorsInverse[8][7] = AHat_s * QStar_z - alphaBar_f * bz
	eigenvectorsInverse[8][8] = alphaBar_f

	print()
	print('error of eigenbasis '..i)
	local evErr = {}
	for i=1,8 do
		evErr[i] = table()
		for j=1,8 do
			local sum = 0
			for k=1,8 do
				sum = sum + eigenvectorsInverse[i][k] * eigenvectors[k][j]
			end
			evErr[i][j] = sum
		end
		print(evErr[i]:concat', ')
	end
end

return MHD

