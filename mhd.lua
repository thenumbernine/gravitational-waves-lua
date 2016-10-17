-- https://arxiv.org/pdf/0804.0402v1.pdf
-- based on Athena's version of eigenvectors of derivative of adiabatic MHD flux wrt primitives
local class = require 'ext.class'
local Equation = require 'equation'

local MHD = class(Equation)
MHD.name = 'MHD'

MHD.numStates = 8	-- including bx 
MHD.numWaves = 7	-- excluding bx
MHD.gamma = 5/3	
MHD.mu = 1

do
	local q = function(self,i) return self.qs[i] end
	local mu = function(self,i) return self.equation.mu end
	local gamma = function(self,i) return self.equation.gamma end
	local rho = q:_(1)
	local mx, my, mz = q:_(2), q:_(3), q:_(4)
	local bx, by, bz = q:_(5), q:_(6), q:_(7)
	local ETotal = q:_(8)
	local vx, vy, vz = mx/rho, my/rho, mz/rho
	local EMag = .5*(bx^2 + by^2 + bz^2)
	local EHydro = ETotal - EMag
	local EKin = .5 * rho * (vx^2 + vy^2 + vz^2)
	local EInt = EHydro - EKin
	local P = (gamma - 1) * EInt	-- hydro pressure
	local H = EInt + P
	local PMag = EMag	-- magnetic pressure
	local PTotal = P + PMag	-- total pressure
	local S = P / rho^gamma	-- mhd entropy the same as non-mhd?
	MHD:buildGraphInfos{
		{rho=rho},
		{vx=vx}, {vy=vy}, {vz=vz},
		--{mx=mx}, {my=my}, {mz=mz},
		{bx=bx}, {by=by}, {bz=bz},
		
		{P=P}, 
		--{PTotal=PTotal},
		{S=S},
		{H=H},
		
		{ETotal=ETotal}, {EKin=EKin}, {EInt=EInt}, {EHydro=EHydro}, {EMag=EMag},
	}
end

function MHD:initCell(sim,i)
	local x = sim.xs[i]
	local rho = x < 0 and 1 or .125
	local vx, vy, vz = 0, 0, 0
	-- [[ Brio & Wu
	self.gamma = 2
	local bx = .75
	local by = x < 0 and 1 or -1
	local bz = 0
	--]]
	--[[ some other tests
	local bx = 1 - math.sin(math.pi/2*x)	-- sine waves break
	local by = math.cos(math.pi/2*x)
	local bz = math.sin(math.pi/2*x)
	--local bx, by, bz = 0, math.sin(math.pi/2*x), 0
	--local bx, by, bz = 0, 1, 0	-- constant field works
	--]]
	--[[ Sod
	local bx, by, bz = 0, 0, 0	-- zero field works ... sort of.
	--]]
	local P = x < 0 and 1 or .1
	return {self:calcConsFromPrim(rho, vx, vy, vz, bx, by, bz, P)}
end

function MHD:calcConsFromPrim(rho, vx, vy, vz, bx, by, bz, P)
	local gamma = self.gamma
	local EInt = P / (gamma-1)
	local eKin = .5*(vx*vx + vy*vy + vz*vz)
	local EKin = rho * eKin
	local bSq = bx*bx + by*by + bz*bz
	local EMag = .5*bSq
	local ETotal = EInt + EKin + EMag 
	local mx, my, mz = rho * vx, rho * vy, rho * vz
	return rho, mx, my, mz, bx, by, bz, ETotal
end

function MHD:calcPrimFromCons(rho, mx, my, mz, bx, by, bz, ETotal)
	local gamma = self.gamma
	local vx, vy, vz = mx / rho, my / rho, mz / rho
	local vSq = vx*vx + vy*vy + vz*vz
	local bSq = bx*bx + by*by + bz*bz
	local EKin = .5 * rho * vSq
	local EMag = .5 * bSq
	local P = (ETotal - EKin - EMag) * (gamma - 1)		
	rho = math.max(rho, 1e-7)
	P = math.max(P, 1e-7)
	return rho, vx, vy, vz, bx, by, bz, P
end

--[[
function MHD:calcMinMaxEigenvaluesFromCons(...)
	local gamma = self.gamma	
	local rho, vx, vy, vz, bx, by, bz, P = self:calcPrimFromCons(...)
	
	local bSq = bx*bx + by*by + bz*bz
	local _1_rho = 1/rho

	local aSq = gamma * P * _1_rho
	local CaxSq = bx * bx * _1_rho
	local CaSq = bSq * _1_rho
	
	local CStarSq = .5 * (CaSq + aSq)
	local sqrtCfsDiscr = math.sqrt(math.max(0, CStarSq * CStarSq - aSq * CaxSq))
	
	local CfSq = CStarSq + sqrtCfsDiscr
	local CsSq = CStarSq - sqrtCfsDiscr

	local Cf = math.sqrt(CfSq)
	local Cs = math.sqrt(math.max(CsSq,0))
	return vx - Cf, vx + Cf	
end
--]]
-- [[
function MHD:calcEigenvaluesFromCons(...)
	local gamma = self.gamma
	local rho, vx, vy, vz, bx, by, bz, P = self:calcPrimFromCons(...)
	local hTotal = .5 * (vx * vx + vy * vy + vz * vz) + (gamma / (gamma - 1) * P + bx * bx + by * by + bz * bz) / rho
	local lambdas = {}
	self:calcEigenBasis(lambdas, nil, nil, nil, rho, vx, vy, vz, bx, by, bz, hTotal, 0, 1)
	return table.unpack(lambdas)
end
--]]

function MHD:calcFluxForState(q)
	local rho, mx, my, mz, bx, by, bz, ETotal = table.unpack(q)
	local rho, vx, vy, vz, bx, by, bz, P = self:calcPrimFromCons(table.unpack(q))
	local bSq = bx*bx + by*by + bz*bz
	local bDotV = bx*vx + by*vy + bz*vz
	local PMag = .5 * bSq
	local PTotal = P + PMag
	local HTotal = ETotal + PTotal
	return 
		mx,
		mx * vx + PTotal - bx * bx,
		mx * vy - bx * by,
		mx * vz - bx * bz,
		0,
		by * vx - bx * vy,
		bz * vx - bx * vz,
		HTotal * vx - bx * bDotV
end

function MHD:calcRoeValues(qL, qR)
	-- should I use bx, or bxL/R, for calculating the PMag at the L and R states?
	local rhoL, vxL, vyL, vzL, bxL, byL, bzL, PL = self:calcPrimFromCons(unpack(qL))
	local ETotalL = qL[8]
	local sqrtRhoL = math.sqrt(rhoL)
	local PMagL = .5 * (bxL * bxL + byL * byL + bzL * bzL)
	local hTotalL = (ETotalL + PL + PMagL)/rhoL
	
	local rhoR, vxR, vyR, vzR, bxR, byR, bzR, PR = self:calcPrimFromCons(unpack(qR))
	local ETotalR = qR[8]
	local sqrtRhoR = math.sqrt(rhoR)
	local PMagR = .5 * (bxR * bxR + byR * byR + bzR * bzR)
	local hTotalR = (ETotalR + PR + PMagR)/rhoR
	
	local invDenom = 1 / (sqrtRhoL + sqrtRhoR)
	local rho  = sqrtRhoL * sqrtRhoR
	local vx = (sqrtRhoL * vxL + sqrtRhoR * vxR) * invDenom
	local vy = (sqrtRhoL * vyL + sqrtRhoR * vyR) * invDenom
	local vz = (sqrtRhoL * vzL + sqrtRhoR * vzR) * invDenom
	local bx = (sqrtRhoL * bxL + sqrtRhoR * bxR) * invDenom
	
	-- why does athena switch the weights of the by and bz components?
	local by = (sqrtRhoR * byL + sqrtRhoL * byR) * invDenom
	local bz = (sqrtRhoR * bzL + sqrtRhoL * bzR) * invDenom
	
	local hTotal = (sqrtRhoL * hTotalL + sqrtRhoR * hTotalR) * invDenom
	local X = .5*((byL - byR)^2 + (bzL - bzR)^2)/((sqrtRhoL + sqrtRhoR)^2)
	local Y = .5*(rhoL + rhoR)/rho
	
	return rho, vx, vy, vz, bx, by, bz, hTotal, X, Y
end

function MHD:calcEigenBasis(lambda, evR, evL, dF_dU, rho, vx, vy, vz, bx, by, bz, hTotal, X, Y)
	local gamma = self.gamma
	local gamma_1 = gamma - 1
	local gamma_2 = gamma - 2
	local gamma_3 = gamma - 3

	local _1_rho = 1 / rho
	local vSq = vx*vx + vy*vy + vz*vz
	local bDotV = bx*vx + by*vy + bz*vz
	local bPerpSq = by*by + bz*bz
	local bStarPerpSq = (gamma_1 - gamma_2 * Y) * bPerpSq
	local CAxSq = bx*bx*_1_rho
	local CASq = CAxSq + bPerpSq * _1_rho
	local hHydro = hTotal - CASq
	-- hTotal = (EHydro + EMag + P)/rho
	-- hHydro = hTotal - CASq, CASq = EMag/rho 
	-- hHydro = eHydro + P/rho = eKin + eInt + P/rho
	-- hHydro - eKin = eInt + P/rho = (1/(gamma-1) + 1) P/rho = gamma/(gamma-1) P/rho
	-- a^2 = (gamma-1)(hHydro - eKin) = gamma P / rho
	local aTildeSq = math.max((gamma_1 * (hHydro - .5 * vSq) - gamma_2 * X), 1e-20)

	local bStarPerpSq_rho = bStarPerpSq * _1_rho
	local CATildeSq = CAxSq + bStarPerpSq_rho
	local CStarSq = .5 * (CATildeSq + aTildeSq)
	local CA_a_TildeSqDiff = .5 * (CATildeSq - aTildeSq)
	local sqrtDiscr = math.sqrt(CA_a_TildeSqDiff * CA_a_TildeSqDiff + aTildeSq * bStarPerpSq_rho)
	
	local CfSq = CStarSq + sqrtDiscr
	local Cf = math.sqrt(CfSq)

	local CsSq = aTildeSq * CAxSq / CfSq
	local Cs = math.sqrt(CsSq)

	local bPerpLen = math.sqrt(bPerpSq)
	local bStarPerpLen = math.sqrt(bStarPerpSq)
	local betaY, betaZ
	if bPerpLen == 0 then
		betaY = 1
		betaZ = 0
	else
		betaY = by / bPerpLen
		betaZ = bz / bPerpLen
	end
	local betaStarY = betaY / math.sqrt(gamma_1 - gamma_2*Y)
	local betaStarZ = betaZ / math.sqrt(gamma_1 - gamma_2*Y)
	local betaStarSq = betaStarY*betaStarY + betaStarZ*betaStarZ
	local vDotBeta = vy*betaStarY + vz*betaStarZ

	local alphaF, alphaS
	if CfSq - CsSq == 0 then
		alphaF = 1
		alphaS = 0
	elseif aTildeSq - CsSq <= 0 then
		alphaF = 0
		alphaS = 1
	elseif CfSq - aTildeSq <= 0 then
		alphaF = 1
		alphaS = 0
	else
		alphaF = math.sqrt((aTildeSq - CsSq) / (CfSq - CsSq))
		alphaS = math.sqrt((CfSq - aTildeSq) / (CfSq - CsSq))
	end

	local sqrtRho = math.sqrt(rho)
	local _1_sqrtRho = 1/sqrtRho
	local sbx = bx >= 0 and 1 or -1
	local aTilde = math.sqrt(aTildeSq)
	local Qf = Cf*alphaF*sbx
	local Qs = Cs*alphaS*sbx
	local Af = aTilde*alphaF*_1_sqrtRho
	local As = aTilde*alphaS*_1_sqrtRho
	local Afpbb = Af*bStarPerpLen*betaStarSq
	local Aspbb = As*bStarPerpLen*betaStarSq

	local CAx = math.sqrt(CAxSq)
	fill(lambda, vx-Cf, vx-CAx, vx-Cs, vx, vx+Cs, vx+CAx, vx+Cf)

	-- dF/dU 
	if dF_dU then
		fill(dF_dU[1], 0, 											1,										0,							0,							0,			0,							0							)
		fill(dF_dU[2], -vx*vx + .5*gamma_1*vSq - gamma_2*X,			-gamma_3*vx,							-gamma_1*vy,				-gamma_1*vz,				gamma_1,	-gamma_2*Y*by,				-gamma_2*Y*bz				)
		fill(dF_dU[3], -vx*vy,										vy,										vx,							0, 							0,			-bx,						0							)
		fill(dF_dU[4], -vx*vz,										vz,										0,							vx, 						0,			0,							-bx							)
		fill(dF_dU[5], vx*(.5*gamma_1*vSq - hTotal) + bx*bDotV/rho,	-gamma_1*vx*vx + hTotal - bx*bx/rho,	-gamma_1*vx*vy - bx*by/rho,	-gamma_1*vx*vz - bx*bz/rho,	gamma*vx,	-gamma_2*Y*by*vx - bx*vy,	-gamma_2*Y*bz*vx - bx*vz	)
		fill(dF_dU[6], (bx*vy - by*vx)/rho,							by/rho,									-bx/rho,					0, 							0,			vx,							0							)
		fill(dF_dU[7], (bx*vz - bz*vx)/rho,							bz/rho,									0,							-bx/rho, 					0,			0,							vx							)
	end

	-- right eigenvectors
	local qa3 = alphaF*vy
	local qb3 = alphaS*vy
	local qc3 = Qs*betaStarY
	local qd3 = Qf*betaStarY
	local qa4 = alphaF*vz
	local qb4 = alphaS*vz
	local qc4 = Qs*betaStarZ
	local qd4 = Qf*betaStarZ
	local r52 = -(vy*betaZ - vz*betaY)
	local r61 = As*betaStarY
	local r62 = -betaZ*sbx*_1_sqrtRho
	local r63 = -Af*betaStarY
	local r71 = As*betaStarZ
	local r72 = betaY*sbx*_1_sqrtRho
	local r73 = -Af*betaStarZ
	if evR then
		fill(evR[1], alphaF, 0, alphaS, 1, alphaS, 0, alphaF)
		fill(evR[2], alphaF*lambda[1], 0, alphaS*lambda[3], vx, alphaS*lambda[5], 0, alphaF*lambda[7])
		fill(evR[3], qa3 + qc3, -betaZ, qb3 - qd3, vy, qb3 + qd3, betaZ, qa3 - qc3)
		fill(evR[4], qa4 + qc4, betaY, qb4 - qd4, vz, qb4 + qd4, -betaY, qa4 - qc4)
		fill(evR[5], alphaF*(hHydro - vx*Cf) + Qs*vDotBeta + Aspbb, r52, alphaS*(hHydro - vx*Cs) - Qf*vDotBeta - Afpbb, .5*vSq + gamma_2*X/gamma_1, alphaS*(hHydro + vx*Cs) + Qf*vDotBeta - Afpbb, -r52, alphaF*(hHydro + vx*Cf) - Qs*vDotBeta + Aspbb)
		fill(evR[6], r61, r62, r63, 0, r63, r62, r61)
		fill(evR[7], r71, r72, r73, 0, r73, r72, r71)
	end

	-- left eigenvectors
	local norm = .5/aTildeSq
	local Cff = norm*alphaF*Cf
	local Css = norm*alphaS*Cs
	Qf = Qf * norm
	Qs = Qs * norm
	local AHatF = norm*Af*rho
	local AHatS = norm*As*rho
	local afpb = norm*Af*bStarPerpLen
	local aspb = norm*As*bStarPerpLen

	norm = norm * gamma_1
	alphaF = alphaF * norm
	alphaS = alphaS * norm
	local QStarY = betaStarY/betaStarSq
	local QStarZ = betaStarZ/betaStarSq
	local vqstr = (vy*QStarY + vz*QStarZ)
	norm = norm * 2
	if evL then
		fill(evL[1], alphaF*(vSq-hHydro) + Cff*(Cf+vx) - Qs*vqstr - aspb, -alphaF*vx - Cff, -alphaF*vy + Qs*QStarY, -alphaF*vz + Qs*QStarZ, alphaF, AHatS*QStarY - alphaF*by, AHatS*QStarZ - alphaF*bz)
		fill(evL[2], .5*(vy*betaZ - vz*betaY), 0, -.5*betaZ, .5*betaY, 0, -.5*sqrtRho*betaZ*sbx, .5*sqrtRho*betaY*sbx)
		fill(evL[3], alphaS*(vSq-hHydro) + Css*(Cs+vx) + Qf*vqstr + afpb, -alphaS*vx - Css, -alphaS*vy - Qf*QStarY, -alphaS*vz - Qf*QStarZ, alphaS, -AHatF*QStarY - alphaS*by, -AHatF*QStarZ - alphaS*bz)
		fill(evL[4], 1 - norm*(.5*vSq - gamma_2*X/gamma_1) , norm*vx, norm*vy, norm*vz, -norm, norm*by, norm*bz)
		fill(evL[5], alphaS*(vSq-hHydro) + Css*(Cs-vx) - Qf*vqstr + afpb, -alphaS*vx + Css, -alphaS*vy + Qf*QStarY, -alphaS*vz + Qf*QStarZ, alphaS, evL[3][6], evL[3][7])
		fill(evL[6], -evL[2][1], 0, -evL[2][3], -evL[2][4], 0, evL[2][6], evL[2][7])
		fill(evL[7], alphaF*(vSq-hHydro) + Cff*(Cf-vx) + Qs*vqstr - aspb, -alphaF*vx + Cff, -alphaF*vy - Qs*QStarY, -alphaF*vz - Qs*QStarZ, alphaF, evL[1][6], evL[1][7])
	end
end

function MHD:calcInterfaceEigenvalues(solver,i)	
	return self:calcEigenBasis(
		solver.eigenvalues[i],
		nil,
		nil,
		nil,
		self:calcInterfaceRoeValues(solver, i))
end

function MHD:calcCellCenterRoeValues(solver, i)
	local q = solver.qs[i]
	local rho, vx, vy, vz, bx, by, bz, P = self:calcPrimFromCons(table.unpack(q))
	local ETotal = q[8]
	local PMag = .5 * (bx * bx + by * by + bz * bz)
	local hTotal = (ETotal + P + PMag)/rho
	return rho, vx, vy, vz, bx, by, bz, hTotal, 0, 1
end

local function permute8to7(v1,v2,v3,v4,v5,v6,v7,v8)
	return v1,v2,v3,v4,v8,v6,v7
end
local function permute7to8(v1,v2,v3,v4,v5,v6,v7)
	return v1,v2,v3,v4,0,v6,v7,v5
end

function MHD:eigenTransform(solver, m, x, from, to)
	-- x is energy-last order, so convert to energy-5th order
	if from then
		x = {permute8to7(table.unpack(x))}
	end
	local y = {}
	for j=1,self.numWaves do
		local sum = 0
		for k=1,self.numWaves do
			sum = sum + m[j][k] * x[k]
		end
		y[j] = sum
	end
	-- convert back to energy-last order
	if to then
		y = {permute7to8(table.unpack(y))}
	end
	return y
end

return MHD
