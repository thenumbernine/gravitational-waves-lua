-- https://arxiv.org/pdf/0804.0402v1.pdf
-- based on Athena's eigenvectors of derivative of adiabatic MHD flux wrt conservative variables
local class = require 'ext.class'
local Equation = require 'equation'

local MHD = class(Equation)

-- whether to use the 7x7 system
local use7x7

MHD.numStates = 8
MHD.gamma = 5/3	
MHD.mu = 1

do
	local q = function(self,i) return self.qs[i] end
	local mu = function(self,i) return self.equation.mu end
	local gamma = function(self,i) return (self.equation.gamma - 1) end
	local rho = q:_(1)
	local mx, my, mz = q:_(2), q:_(3), q:_(4)
	local bx, by, bz = q:_(5), q:_(6), q:_(7)
	local ETotal = q:_(8)
	local vx, vy, vz = mx/rho, my/rho, mz/rho
	local EMag = .5*(bx^2 + by^2 + bz^2) / mu
	local EHydro = ETotal - EMag
	local EKin = .5 * rho * (vx^2 + vy^2 + vz^2)
	local EInt = EHydro - EKin
	local P = (gamma - 1) * EInt
	local PStar = P + EMag
	local S = P / rho^gamma	-- mhd entropy the same as non-mhd?
	-- [[ full mhd
	MHD:buildGraphInfos{
		{rho=rho},
		{vx=vx}, {vy=vy}, {vz=vz},
		{bx=bx}, {by=by}, {bz=bz},
		{P=P}, {PStar=PStar},
		{ETotal=ETotal}, {EKin=EKin}, {EInt=EInt}, {EHydro=EHydro}, {EMag=EMag},
		{['eigenbasis error'] = function(self,i) return self.eigenbasisErrors and math.log(self.eigenbasisErrors[i]) end},
		{['reconstruction error'] = function(self,i) return self.fluxMatrixErrors and math.log(self.fluxMatrixErrors[i]) end},
	}
	--]]
	--[[ matching Euler
	MHD:buildGraphInfos{
		{rho=rho},
		{vx=vx},
		{P=P},
		{S=S},
		{mom=mx},
		{ETotal=ETotal},
		{['eigenbasis error'] = function(self,i) return self.eigenbasisErrors and math.log(self.eigenbasisErrors[i]) end},
		{['reconstruction error'] = function(self,i) return self.fluxMatrixErrors and math.log(self.fluxMatrixErrors[i]) end},
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
	-- [[ some other tests
	--local bx, by, bz = 0, math.sin(math.pi/2*x), 0
	--local bx, by, bz = 1, 0, 0	-- constant x field 
	--local bx, by, bz = 0, 1, 0	-- constant y field
	local bx, by, bz = 0, 1, 0	-- constant z field
	--]]
	--[[ Sod
	local bx, by, bz = 0, 0, 0	-- zero field works ... sort of.
	--]]
	local P = x < 0 and 1 or .1
	local EInt = P / (gamma-1)
	local eKin = .5*(vx*vx + vy*vy + vz*vz)
	local EKin = rho * eKin
	local BSq = bx*bx + by*by + bz*bz
	local EMag = .5*BSq
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

-- 
function MHD:calcMinMaxEigenvaluesFromCons(rho, mx, my, mz, bx, by, bz, ETotal)
	local rho, vx, vy, vz, bx, by, bz, H = self:stateToPrims(rho, mx, my, mz, bx, by, bz, ETotal)
	local vSq = vx*vx + vy*vy + vz*vz	
	local bSq = bx*bx + by*by + bz*bz
	local gamma = self.gamma	
	local gammaPrime = gamma - 1
	-- this usually goes on at interfaces, hence the X and Y calculations using L and R states
	-- but the same paper that says how to calc the eigenvectors, also says to use cell-centered eigenvalues for the cfl timestep
	local X = 0--((byL - byR)^2 + (bzL - bzR)^2) / (2 * (math.sqrt(rhoL) + math.sqrt(rhoR)))
	local XPrime = (gamma - 2) * X
	local Y = 1--(rhoL + rhoR) / (2 * rho)
	local YPrime = (gamma - 2) * Y
	local aTildeSq = gammaPrime * (H - vSq/2 - bSq/rho) - XPrime
	local CAxSq = bx*bx / rho
	local bPerpSq = by*by + bz*bz
	local bStarPerpSq = (gammaPrime - YPrime) * bPerpSq
	local CATildeSq = CAxSq + bStarPerpSq / rho
	local CAx = math.sqrt(CAxSq)
	local CfSq = .5 * ((aTildeSq + CATildeSq) + math.sqrt((aTildeSq + CATildeSq)^2 - 4 * aTildeSq * CAxSq))
	local Cf = math.sqrt(CfSq)
	return vx - Cf, vx + Cf
end

function MHD:calcInterfaceEigenBasis(sim,i,qL,qR)
	local gamma = self.gamma	
	local gammaMinusOne = gamma - 1

	local rhoL, vxL, vyL, vzL, bxL, byL, bzL, HL = self:stateToPrims(unpack(qL))
	local rhoR, vxR, vyR, vzR, bxR, byR, bzR, HR = self:stateToPrims(unpack(qR))

	-- eqn 56 - averaging
	local sqrtRhoL = math.sqrt(rhoL)
	local sqrtRhoR = math.sqrt(rhoR)
	local invDenom = 1 / (sqrtRhoL + sqrtRhoR)
	local rho = sqrtRhoL * sqrtRhoR
	local vx = (sqrtRhoL*vxL + sqrtRhoR*vxR)*invDenom
	local vy = (sqrtRhoL*vyL + sqrtRhoR*vyR)*invDenom
	local vz = (sqrtRhoL*vzL + sqrtRhoR*vzR)*invDenom
	local bx = (sqrtRhoL*bxL + sqrtRhoR*bxR)*invDenom
	local by = (sqrtRhoL*byL + sqrtRhoR*byR)*invDenom
	local bz = (sqrtRhoL*bzL + sqrtRhoR*bzR)*invDenom
	local H = (sqrtRhoL*HL + sqrtRhoR*HR)*invDenom
assertfinite(sqrtRhoL)
assertfinite(sqrtRhoR)
assertfinite(invDenom)
assertfinite(rho)
assertfinite(vx)
assertfinite(vy)
assertfinite(vz)
assertfinite(bx)
assertfinite(by)
assertfinite(bz)
assertfinite(H)

	-- eqn 56 says to average H
	-- but B8 says H = (E + P + b^2/2)/rho

	local vSq = vx*vx + vy*vy + vz*vz
assertfinite(vSq)

	local sqrtRho = math.sqrt(rho)
assertfinite(sqrtRho)

	-- TODO divide by mu instead?
	local oneOverSqrtMu = 1 / math.sqrt(4 * math.pi)
assertfinite(oneOverSqrtMu)
	local bSq = bx*bx + by*by + bz*bz
assertfinite(bSq)
	local bDotV = bx*vx + by*vy + bz*vz
assertfinite(bDotV)

	local gammaPrime = gamma - 1
	local X = ((byL - byR)^2 + (bzL - bzR)^2) / (2 * (math.sqrt(rhoL) + math.sqrt(rhoR)))
assertfinite(X)
	local XPrime = (gamma - 2) * X
assertfinite(XPrime)
	local Y = (rhoL + rhoR) / (2 * rho)
assertfinite(Y)
	local YPrime = (gamma - 2) * Y
assertfinite(YPrime)

if use7x7 then
	-- , paper-order
	-- storing these matrices as 7x7 just as they are in the paper (i.e. rho, v, E, by, bz)
	-- then in the mat mul function I'm going to rearrange the conservative vector before and after applying it 
	local A51 = -vx*H + gammaPrime * vx * vSq/2 + bx * (bx * vx + by * vy + bz * vz) / rho - vx * XPrime
	local A52 = -gammaPrime * vx^2 + H - bx^2 - bx^2/rho
	local A53 = -gammaPrime * vx * vy - bx * by / rho
	local A54 = -gammaPrime * vx * vz - bx * bz / rho
	local A56 = -(bx * vy + by * vx * YPrime)
	local A57 = -(bx * vz + bz * vx * YPrime)
	sim.fluxMatrix[i] = {
		{0, 1, 0, 0, 0, 0, 0},
		{-vx^2+gammaPrime*vSq/2-XPrime, -(gamma-3)*vx, -gammaPrime*vy, -gammaPrime*vz, gammaPrime, -by*YPrime, -bz*YPrime},
		{-vx*vy, vy, vx, 0, 0, -bx, 0},
		{-vx*vz, vz, 0, vx, 0, 0, -bx},
		{A51, A52, A53, A54, gamma*vx, A56, A57},
		{(bx*vy-by*vx)/rho, by/rho, -bx/rho, 0, 0, vx, 0},
		{(bx*vz-bz*vx)/rho, bz/rho, 0, -bx/rho, 0, 0, vx},
	}
else
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
end

	local aTildeSq = gammaPrime * (H - vSq/2 - bSq/rho) - XPrime
aTildeSq = math.max(1e-7, aTildeSq)
--if aTildeSq < 0 then error(tolua({aTildeSq=aTildeSq, H=H, vSq=vSq, bSq=bSq, rho=rho, XPrime=XPrime}, {indent=true})) end
	local CAxSq = bx*bx / rho
assertfinite(CAxSq)
	local bPerpSq = by*by + bz*bz
assertfinite(bPerpSq)
	local bStarPerpSq = (gammaPrime - YPrime) * bPerpSq
assertfinite(bStarPerpSq)
	local CATildeSq = CAxSq + bStarPerpSq / rho
assertfinite(CATildeSq)
	local CAx = math.sqrt(CAxSq)
assertfinite(CAx)
	local CfSq = .5 * ((aTildeSq + CATildeSq) + math.sqrt((aTildeSq + CATildeSq)^2 - 4 * aTildeSq * CAxSq))
assertfinite(CfSq)
	local CsSq = .5 * ((aTildeSq + CATildeSq) - math.sqrt((aTildeSq + CATildeSq)^2 - 4 * aTildeSq * CAxSq))
assertfinite(CsSq)
	local Cf = math.sqrt(CfSq)
assertfinite(Cf)
	local Cs = math.sqrt(CsSq)
assertfinite(Cs)

	local eigenvalues = sim.eigenvalues[i]
if use7x7 then -- paper-order
	eigenvalues[1] = vx - Cf
	eigenvalues[2] = vx - CAx
	eigenvalues[3] = vx - Cs
	eigenvalues[4] = vx
	eigenvalues[5] = 0
	eigenvalues[5] = vx + Cs
	eigenvalues[6] = vx + CAx
	eigenvalues[7] = vx + Cf
else -- my-order
	eigenvalues[1] = vx - Cf
	eigenvalues[2] = vx - CAx
	eigenvalues[3] = vx - Cs
	eigenvalues[4] = vx
	eigenvalues[5] = vx
	eigenvalues[6] = vx + Cs
	eigenvalues[7] = vx + CAx
	eigenvalues[8] = vx + Cf
end

for j=1,#eigenvalues do
assertfinite(eigenvalues[j])
end

	-- to prevent divide-by-zero errors
	local epsilon = 1e-20

	-- eqn A13-A17, replace 'a' with 'aTilde'
	local aTilde = math.sqrt(aTildeSq)
if not math.isfinite(aTilde) then error(tolua({aTilde=aTilde, aTildeSq=aTildeSq}, {indent=true})) end
	local S = bx >= 0 and 1 or -1
assertfinite(S)
	local alpha_f, alpha_s
	if CfSq == CsSq then-- CaSq == aTildeSq and CAxSq == aTildeSq then
		alpha_f = 1
		alpha_s = 0
	else
		alpha_f = (aTildeSq - CsSq) / (CfSq - CsSq)
		alpha_s = (CfSq - aTildeSq) / (CfSq - CsSq)
	end
assertfinite(alpha_s)
assertfinite(alpha_f)
	local Cff = Cf * alpha_f
assertfinite(Cff)
	local Css = Cs * alpha_s
assertfinite(Css)
	local Qf = Cf * alpha_f * S
assertfinite(Qf)
	local Qs = Cs * alpha_s * S
assertfinite(Qs)
	local Af = aTilde * alpha_f * sqrtRho
assertfinite(Af)
	local As = aTilde * alpha_s * sqrtRho
assertfinite(As)
	local bPerp = math.sqrt(bPerpSq)
assertfinite(bPerp)
	local beta_y, beta_z
	if bPerp < epsilon then
		beta_y = 1
		beta_z = 0
	else
		beta_y = by / bPerp
		beta_z = bz / bPerp
	end
assertfinite(beta_y)
assertfinite(beta_z)

	-- V[xyz][fs] = v[xyz] * alpha_[fs]
	local Vxf = vx * alpha_f
assertfinite(Vxf)
	local Vyf = vy * alpha_f
assertfinite(Vyf)
	local Vzf = vz * alpha_f
assertfinite(Vzf)
	local Vxs = vx * alpha_s
assertfinite(Vxs)
	local Vys = vy * alpha_s
assertfinite(Vys)
	local Vzs = vz * alpha_s
assertfinite(Vzs)

	local bStarPerp = math.sqrt(bStarPerpSq)
assertfinite(bStarPerp)
	local betaStar_y, betaStar_z
	if bStarPerp < epsilon then
		betaStar_y = 1
		betaStar_z = 0
	else
		betaStar_y = by / bStarPerp
		betaStar_z = bz / bStarPerp
	end
assertfinite(betaStar_y)
assertfinite(betaStar_z)
	local betaStarPerpSq = betaStar_y * betaStar_y + betaStar_z * betaStar_z
assertfinite(betaStarPerpSq)
	local HPrime = H - bSq / rho		-- (densitized) total gas enthalpy (internal energy + kinetic energy + pressure)
assertfinite(HPrime)

if use7x7 then
	local R51 = alpha_f * (HPrime - vx * Cf) + Qs * (vy * betaStar_y + vz * betaStar_z) + As * bStarPerp * betaStarPerpSq / rho
	local R52 = -(vy * beta_z - vz * beta_y)
	local R53 = alpha_s * (HPrime - vx * Cs) - Qf * (vy * betaStar_y + vz * betaStar_z) - Af * bStarPerp * betaStarPerpSq / rho
	local R54 = vSq/2 + XPrime / gammaPrime
	local R55 = alpha_s * (HPrime + vx * Cs) + Qf * (vy * betaStar_y + vz * betaStar_z) - Af * bStarPerp * betaStarPerpSq / rho
	local R56 = -R52
	local R57 = alpha_f * (HPrime + vx * Cf) - Qs * (vy * betaStar_y + vz * betaStar_z) + As * bStarPerp * betaStarPerpSq / rho
	sim.eigenvectors[i] = {
		{alpha_f, 0, alpha_s, 1, alpha_s, 0, alpha_f},
		{Vxf-Cff, 0, Vxs-Css, vx, Vxs+Css, 0, Vxf+Cff},
		{Vyf+Qs*betaStar_y, -beta_z, Vys-Qf*betaStar_y, vy, Vys+Qf*betaStar_y, beta_z, Vyf-Qs*betaStar_y},
		{Vzf+Qs*betaStar_z, beta_y, Vzs-Qf*betaStar_z, vz, Vzs+Qf*betaStar_z, -beta_y, Vzf-Qs*betaStar_z},
		{R51, R52, R53, R54, R55, R56, R57},
		{As*betaStar_y/rho, -beta_z*S/sqrtRho, -Af*betaStar_y/rho, 0, -Af*betaStar_y/rho, -beta_z*S/sqrtRho, As*betaStar_y/rho},
		{As*betaStar_z/rho, beta_y*S/sqrtRho, -Af*betaStar_z/rho, 0, -Af*betaStar_z/rho, beta_y*S/sqrtRho, As*betaStar_z/rho},
	}
else
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
end
	local evR = sim.eigenvectors[i]
	
	for j=1,#evR do
		for k=1,#evR do
			assertfinite(evR[j][k])
		end
	end

	local QStar_y, QStar_z
	if betaStarPerpSq < epsilon then
		QStar_y = 1
		QStar_z = 0
	else
		QStar_y = betaStar_y / betaStarPerpSq
		QStar_z = betaStar_z / betaStarPerpSq
	end
assertfinite(QStar_y)
assertfinite(QStar_z)

	--local aSq = gamma * P / rho	-- adiabatic sound speed
	local aSq = gammaPrime * (HPrime - .5 * vSq)
assertfinite(aSq)
	local hatScalar = 1 / (2 * aSq)
assertfinite(hatScalar)
	local QHat_f = Qf * hatScalar
assertfinite(QHat_f)
	local QHat_s = Qs * hatScalar
assertfinite(QHat_s)
	local AHat_f = Af * hatScalar	
assertfinite(AHat_f)
	local AHat_s = As * hatScalar	
assertfinite(AHat_s)
	local CHat_ff = Cff * hatScalar
assertfinite(CHat_ff)
	local CHat_ss = Css * hatScalar
assertfinite(CHat_ss)
	local XHatPrime = XPrime * hatScalar
assertfinite(XHatPrime)
	
	local barScalar = gammaPrime * hatScalar
assertfinite(barScalar)
	local alphaBar_f = alpha_f * barScalar
assertfinite(alphaBar_f)
	local alphaBar_s = alpha_s * barScalar
assertfinite(alphaBar_s)
	local VBar_xf = Vxf * barScalar
assertfinite(VBar_xf)
	local VBar_yf = Vyf * barScalar
assertfinite(VBar_yf)
	local VBar_zf = Vzf * barScalar
assertfinite(VBar_zf)
	local VBar_xs = Vxs * barScalar
assertfinite(VBar_xs)
	local VBar_ys = Vys * barScalar
assertfinite(VBar_ys)
	local VBar_zs = Vzs * barScalar
assertfinite(VBar_zs)
	local vBar_x = vx * barScalar
assertfinite(vBar_x)
	local vBar_y = vy * barScalar
assertfinite(vBar_y)
	local vBar_z = vz * barScalar
assertfinite(vBar_z)
	local bBar_y = by * barScalar
assertfinite(bBar_y)
	local bBar_z = bz * barScalar
assertfinite(bBar_z)
	
	local vBarSq = vSq * barScalar * barScalar
assertfinite(vBarSq)

if use7x7 then
	local L11 = alphaBar_f*(vSq-HPrime)+CHat_ff*(Cf+vx)-QHat_s*(vy*QStar_y+vz*QStar_z)-AHat_s*bPerp/rho
	local L21 = (vy*beta_z-vz*beta_y)/2
	local L31 = alphaBar_s*(vSq-HPrime)+CHat_ss*(Cs+vx)+QHat_f*(vy*QStar_y+vz*QStar_z)+AHat_f*bPerp/rho
	local L41 = 1 - vBarSq + 2 * XHatPrime^2
	local L51 = alphaBar_s*(vSq-HPrime)+CHat_ss*(Cs-vx)-QHat_f*(vy*QStar_y+vz*QStar_z)+AHat_f*bPerp/rho
	local L61 = -L21
	local L71 = alphaBar_f*(vSq-HPrime)+CHat_ff*(Cf-vx)+QHat_s*(vy*QStar_y+vz*QStar_z)-AHat_s*bPerp/rho
	sim.eigenvectorsInverse[i] = {
		{L11, -VBar_xf-CHat_ff, -VBar_yf+QHat_s*QStar_y, -VBar_zf+QHat_s*QStar_z, alphaBar_f, AHat_s*QStar_y-alphaBar_f*by, AHat_s*QStar_z-alphaBar_f*bz},
		{L21, 0, -beta_z/2, beta_y/2, 0, -beta_z*S*sqrtRho/2, beta_y*S*sqrtRho/2},
		{L31, -VBar_xs-CHat_ss, -VBar_ys-QHat_f*QStar_y, -VBar_zs-QHat_f*QStar_z, alphaBar_s, -AHat_f*QStar_y-alphaBar_s*by, -AHat_f*QStar_z-alphaBar_s*bz},
		{L41, 2*vBar_x, 2*vBar_y, 2*vBar_z, -gammaPrime/aSq, 2*bBar_y, 2*bBar_z},
		{L51, -VBar_xs+CHat_ss, -VBar_ys+QHat_f*QStar_y, -VBar_zs+QHat_f*QStar_z, alphaBar_s, -AHat_f*QStar_y-alphaBar_s*by, -AHat_f*QStar_z-alphaBar_s*bz},
		{L61, 0, beta_z/2, -beta_y/2, 0, -beta_z*S*sqrtRho/2, beta_y*S*sqrtRho/2},
		{L71, -VBar_xf+CHat_ff, -VBar_yf-QHat_s*QStar_y, -VBar_zf-QHat_s*QStar_z, alphaBar_f, AHat_s*QStar_y-alphaBar_f*by, AHat_s*QStar_z-alphaBar_f*bz},
	}
else
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
end
	local evL = sim.eigenvectorsInverse[i]

	for j=1,#evR do
		for k=1,#evR do
			assertfinite(evL[j][k])
		end
	end

-- [[
	print()
	print('error of eigenbasis '..i)
	local errTotal = 0
	local evErr = {}
	for j=1,#evR do
		evErr[j] = table()
		for k=1,#evR do
			local sum = 0
			for l=1,#evR do
				sum = sum + evL[j][l] * evR[l][k]
			end
			evErr[j][k] = sum
			errTotal = errTotal + math.abs(sum - (j == k and 1 or 0))
		end
		print(evErr[j]:concat', ')
	end
	print('error total',errTotal)
--]]
end

if use7x7 then
local function permute(cons8)
	local rho, mx, my, mz, bx, by, bz, E = table.unpack(cons8)
	return {rho, mx, my, mz, E, by, bz}, bx
end

local function unpermute(cons7, bx)
	local rho, mx, my, mz, E, by, bz = table.unpack(cons7)
	return {rho, mx, my, mz, bx, by, bz, E}
end

local function applyPermuteMatrixField(field)
	return function(self, solver, i, vOrig)
		local v, bx = permute(vOrig)
		local m = solver[field][i]
		local result = {}
		for j=1,7 do
			local sum = 0
			for k=1,7 do
				sum = sum + m[j][k] * v[k]
			end
			result[j] = sum
		end
		return unpermute(result, bx)
	end
end

MHD.fluxTransform = applyPermuteMatrixField'fluxMatrix'
MHD.applyLeftEigenvectors = applyPermuteMatrixField'eigenvectorsInverse'
MHD.applyRightEigenvectors = applyPermuteMatrixField'eigenvectors'
end

function MHD:postIterate(sim)
	if MHD.super.postIterate then MHD.super.postIterate(self, sim) end
	for i=1,sim.gridsize do
		local q = sim.qs[i]
		q[1] = math.max(q[1], 1e-7)	-- min bound for rho
	end
end

return MHD
