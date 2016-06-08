-- https://arxiv.org/pdf/0804.0402v1.pdf
-- based on Athena's version of eigenvectors of derivative of adiabatic MHD flux wrt primitives
-- then dF/dW transformed to dF/dU
-- and with conservatives rearranged so ETotal is last
-- and converted to 8x8 instead of 7x7
local class = require 'ext.class'
local Equation = require 'equation'

local MHD = class(Equation)
MHD.name = 'MHD'

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
	local P = (gamma - 1) * EInt	-- hydro pressure
	local PMag = EMag	-- magnetic pressure
	local PTotal = P + PMag	-- total pressure
	local S = P / rho^gamma	-- mhd entropy the same as non-mhd?
	-- [[ full mhd
	MHD:buildGraphInfos{
		{rho=rho},
		{vx=vx}, {vy=vy}, {vz=vz},
		{bx=bx}, {by=by}, {bz=bz},
		
		{P=P}, 
		--{PTotal=PTotal},
		{S=S},
		
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
	--local b = 1
	--local b = 1/math.sqrt(4*math.pi)
	local b = 1
	local bx = b
	local by = x < 0 and b or -b
	local bz = 0
	--]]
	-- [[ some other tests
	local bx, by, bz = 0, math.sin(math.pi/2*x), 0
	--local bx, by, bz = 0, 1, 0	-- constant field works
	--]]
	--[[ Sod
	local bx, by, bz = 0, 0, 0	-- zero field works ... sort of.
	--]]
	local P = x < 0 and 1 or .1
	return {self:consFromPrim(rho, vx, vy, vz, bx, by, bz, P)}
end

function MHD:consFromPrim(rho, vx, vy, vz, bx, by, bz, P)
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

local _1_sqrt2 = 1/math.sqrt(2)
 
function MHD:calcMinMaxEigenvaluesFromCons(rho, mx, my, mz, bx, by, bz, ETotal)
	local gamma = self.gamma	
	local rho, vx, vy, vz, bx, by, bz, P = self:calcPrimFromCons(rho, mx, my, mz, bx, by, bz, ETotal)
	
	local bSq = bx*bx + by*by + bz*bz
	local _1_rho = 1/rho

	local aSq = gamma * P * _1_rho
	local CaxSq = bx * bx * _1_rho
assertfinite(CaxSq)
	local CaSq = bSq * _1_rho
assertfinite(CaSq)
	
	local CStarSq = .5 * (CaSq + aSq)
assertfinite(CStarSq)
	local sqrtCfsDiscr = math.sqrt(math.max(0, CStarSq * CStarSq - aSq * CaxSq))
assertfinite(sqrtCfsDiscr)
	
	local CfSq = CStarSq + sqrtCfsDiscr
assertfinite(CfSq)
	local CsSq = CStarSq - sqrtCfsDiscr
assertfinite(CsSq)

	local Cf = math.sqrt(CfSq)
	local Cs = math.sqrt(math.max(CsSq,0))

	return vx - Cf, vx + Cf	
end

function MHD:calc_hTotal(rho, PTotal, ETotal)
	return (ETotal + PTotal) / rho
end
		
function MHD:calcSpeedOfSoundSq(rho, vSq, bSq, hTotal)
	return (self.gamma - 1) * (hTotal - .5 * vSq - bSq / rho)
end

--[[
with arithmetic averaging, rho, v, b, and P are averaged.  hTotal and aSq are derived.
with Roe averaging, rho v b and hTotal are averaged.  aSq is derived
hTotal itself is only ever used by dF/dU, which is only used for verifying the eigenbasis reconstruction
 (which I haven't got to ever work, symbolically or numerically, using anyone's eigenvectors of MHD flux)
--]]
MHD.useArithmeticAvergingForRoeValues = true 
function MHD:calcRoeValues(qL, qR)
	local gamma = self.gamma
	local rhoL, vxL, vyL, vzL, bxL, byL, bzL, PL = self:calcPrimFromCons(unpack(qL))
	local rhoR, vxR, vyR, vzR, bxR, byR, bzR, PR = self:calcPrimFromCons(unpack(qR))

	local rho, vx, vy, vz, bx, by, bz
	local vSq, bSq, hTotal, aSq
	
	if self.useArithmeticAvergingForRoeValues then
		rho = .5 * (rhoL + rhoR)
		vx = .5 * (vxL + vxR)
		vy = .5 * (vyL + vyR)
		vz = .5 * (vzL + vzR)
		bx = .5 * (bxL + bxR)
		by = .5 * (byL + byR)
		bz = .5 * (bzL + bzR)
		local P = .5 * (PL + PR) 

		vSq = vx*vx + vy*vy + vz*vz
		bSq = bx*bx + by*by + bz*bz
		local ETotal = P / (gamma - 1) + .5 * rho * vSq + .5 * bSq
		local PTotal = P + .5 * bSq
		hTotal = self:calc_hTotal(rho, PTotal, ETotal)
		aSq = self:calcSpeedOfSoundSq(rho, vSq, bSq, hTotal)
	else
		local PMagL = .5 * (bxL*bxL + byL*byL + bzL*bzL)
		local PTotalL = PL + PMagL
		local ETotalL = qL[8]
		local hTotalL = self:calc_hTotal(rhoL, PTotalL, ETotalL)
		local sqrtRhoL = math.sqrt(rhoL)
		
		local PMagR = .5 * (bxR*bxR + byR*byR + bzR*bzR)
		local PTotalR = PR + PMagR
		local ETotalR = qR[8]
		local hTotalR = self:calc_hTotal(rhoR, PTotalR, ETotalR)
		local sqrtRhoR = math.sqrt(rhoR)
		
		local invDenom = 1 / (sqrtRhoL + sqrtRhoR)
		rho = sqrtRhoL * sqrtRhoR
		vx = (sqrtRhoL*vxL + sqrtRhoR*vxR)*invDenom
		vy = (sqrtRhoL*vyL + sqrtRhoR*vyR)*invDenom
		vz = (sqrtRhoL*vzL + sqrtRhoR*vzR)*invDenom
		bx = (sqrtRhoL*bxL + sqrtRhoR*bxR)*invDenom
		by = (sqrtRhoL*byL + sqrtRhoR*byR)*invDenom
		bz = (sqrtRhoL*bzL + sqrtRhoR*bzR)*invDenom
		hTotal = (sqrtRhoL*hTotalL + sqrtRhoR*hTotalR)*invDenom
		vSq = vx*vx + vy*vy + vz*vz
		bSq = bx*bx + by*by + bz*bz
		aSq = self:calcSpeedOfSoundSq(rho, vSq, bSq, hTotal)
	end
assertfinite(rho)
assertfinite(vx)
assertfinite(vy)
assertfinite(vz)
assertfinite(bx)
assertfinite(by)
assertfinite(bz)
assertfinite(hTotal)
assertfinite(aSq)
assertfinite(vSq)
assertfinite(bSq)

	return rho, vx, vy, vz, bx, by, bz, hTotal, aSq, vSq, bSq
end

function MHD:calcInterfaceEigenBasis(sim,i,qL,qR)
	local gamma = self.gamma	

	local rho, vx, vy, vz, bx, by, bz, hTotal, aSq, vSq, bSq = self:calcRoeValues(qL,qR)

	local sqrtRho = math.sqrt(rho)
	local _1_sqrtRho = 1 / sqrtRho
	local _1_rho = _1_sqrtRho * _1_sqrtRho
assertfinite(_1_rho)
	
	local b_v = bx*vx + by*vy + bz*vz

	local CaxSq = bx*bx*_1_rho
assertfinite(CaxSq)
	local CaSq = bSq*_1_rho
assertfinite(CaSq)
	
	local CStarSq = .5 * (CaSq + aSq)
assertfinite(CStarSq)
	local sqrtCfsDiscr = math.sqrt(math.max(0, CStarSq * CStarSq - aSq * CaxSq))
assertfinite(sqrtCfsDiscr)
	
	local CfSq = CStarSq + sqrtCfsDiscr
assertfinite(CfSq)
	local CsSq = CStarSq - sqrtCfsDiscr
assertfinite(CsSq)
	-- should negative slow speeds be allowed to exist?
	-- even influence alpha calculations?
	-- why clamp Cs to zero a few lines below but not do so here?

	local invDenom = 1 / (CfSq - CsSq)
assertfinite(invDenom)
	local alphaS = math.sqrt((CfSq - aSq) * invDenom)
	local alphaF = math.sqrt((aSq - CsSq) * invDenom)
	if alphaS < 1e-7 then alphaS = 0 end
	if alphaF < 1e-7 then alphaF = 0 end
assertfinite(alphaS)
assertfinite(alphaF)
	
	local sbx = bx >= 0 and 1 or -1
	
	local betaPerpLen = math.sqrt(by*by + bz*bz)
assertfinite(betaPerpLen)
	local betaY, betaZ
	if betaPerpLen < 1e-12 then
		betaY = _1_sqrt2
		betaZ = _1_sqrt2
	else
		local _1_betaPerpLen = 1/betaPerpLen
		betaY = by * _1_betaPerpLen
		betaZ = bz * _1_betaPerpLen
	end
assertfinite(betaY)
assertfinite(betaZ)
	
	local Cf = math.sqrt(CfSq)
assertfinite(Cf)
	local Cs = math.sqrt(math.max(CsSq, 0))
assertfinite(Cs)
	local Cax = math.sqrt(CaxSq)
assertfinite(Cax)
	local a = math.sqrt(aSq)
	
	local _1_a = 1 / a
	local _1_aSq = _1_a * _1_a 

	local gamma_1 = gamma - 1
	local gamma_2 = gamma - 2
	local gamma_3 = gamma - 3
	local dF_dU = sim.fluxMatrix[i]
	--[[ identity for bx dimension 
	local Au25 = 0
	local Au35 = 0
	local Au45 = 0
	local Au55 = 1
	local Au65 = 0
	local Au75 = 0
	local Au85 = 0
	--]]
	--[[ analytical calculations for dF/bx column which coincide with the 8-variable flux, which itself coincides with the dF/by and dF/bz flux derivative columns used in the Athena paper 
	local Au25 = -gamma * bx
	local Au35 = -by
	local Au45 = -bz
	local Au55 = 0
	local Au65 = -vy
	local Au75 = -vz
	local Au85 = -gamma_1 * vx * bx - b_v
	--]]
	-- [[ Athena's 7x7, rearranged to put P last
	local Au25 = -gamma_1*bx
	local Au35 = 0
	local Au45 = 0
	local Au55 = vx
	local Au65 = 0
	local Au75 = 0
	local Au85 = -gamma_1*bx*vx
	--]]
	fill(dF_dU[1], 0, 											1,										0,							0,							0,		0,						0,						0		)
	fill(dF_dU[2], -vx*vx + .5*gamma_1*vSq,						-gamma_3*vx,							-gamma_1*vy,				-gamma_1*vz,				Au25,	-gamma_2*by,			-gamma_2*bz,			gamma_1	)
	fill(dF_dU[3], -vx*vy,										vy,										vx,							0, 							Au35,	-bx,					0,						0		)
	fill(dF_dU[4], -vx*vz,										vz,										0,							vx, 						Au45,	0,						-bx,					0		)
	fill(dF_dU[5], 0,											0,										0,							0, 							Au55,	0,						0,						0		)
	fill(dF_dU[6], (bx*vy - by*vx)/rho,							by/rho,									-bx/rho,					0, 							Au65,	vx,						0,						0		)
	fill(dF_dU[7], (bx*vz - bz*vx)/rho,							bz/rho,									0,							-bx/rho, 					Au75,	0,						vx,						0		)
	fill(dF_dU[8], vx*(.5*gamma_1*vSq - hTotal) + bx*b_v/rho,	-gamma_1*vx*vx + hTotal - bx*bx/rho,	-gamma_1*vx*vy - bx*by/rho,	-gamma_1*vx*vz - bx*bz/rho,	Au85,	-gamma_2*by*vx - bx*vy,	-gamma_2*bz*vx - bx*vz,	gamma*vx)

	local eigenvalues = sim.eigenvalues[i]
	eigenvalues[1] = vx - Cf
	eigenvalues[2] = vx - Cax
	eigenvalues[3] = vx - Cs
	eigenvalues[4] = vx
	eigenvalues[5] = vx
	eigenvalues[6] = vx + Cax
	eigenvalues[7] = vx + Cs
	eigenvalues[8] = vx + Cf

for j=1,8 do
	assertfinite(eigenvalues[j])
end

	--[[
	these are eigenvectors of the flux wrt the primitive variables
	to convert them to eigenvectors of flux wrt conservatives, use this:
		conservative:	
	du/dt + df/dx = 0
	du/dt + df/du du/dx = 0
	du/dt + Au du/dx = 0
	...for Au = df/du
		primitive:
	du/dw dw/dt + df/du du/dw dw/dx = 0
	dw/dt + dw/du df/du du/dw dw/dx = 0
	dw/dt + Aw dw/dx = 0
	...for Aw = dw/du df/du du/dw = dw/du Au du/dw
	so Au = du/dw Aw dw/du 
	and if Aw = Rw Lambda Lw by eigendecomposition
	then Au = (du/dw Rw) Lambda (Lw dw/du)
	and if Au = Ru Lambda Lu
	then Ru = du/dw Rw
	and Lu = Lw dw/du
	--]]

	local evrw = {
		{rho*alphaF,				0,					rho*alphaS,					1,	0,	rho*alphaS,					0,					rho*alphaF				},
		{-alphaF*Cf,				0,					-alphaS*Cs,					0,	0,	alphaS*Cs,					0,					alphaF*Cf				},
		{alphaS*Cs*sbx*betaY,		-betaZ,				-alphaF*Cf*sbx*betaY,		0,	0,	alphaF*Cf*sbx*betaY,		betaZ,				-alphaS*Cs*sbx*betaY	},
		{alphaS*Cs*sbx*betaZ,		betaY,				-alphaF*Cf*sbx*betaZ,		0,	0,	alphaF*Cf*sbx*betaZ,		-betaY,				-alphaS*Cs*sbx*betaZ	},
		{0,							0,					0,							0,	1,	0,							0,					0						},
		{sqrtRho*alphaS*a*betaY,	-sqrtRho*sbx*betaZ,	-sqrtRho*alphaF*a*betaY,	0,	0,	-sqrtRho*alphaF*a*betaY,	-sqrtRho*sbx*betaZ,	sqrtRho*alphaS*a*betaY	},
		{sqrtRho*alphaS*a*betaZ,	sqrtRho*sbx*betaY,	-sqrtRho*alphaF*a*betaZ,	0,	0,	-sqrtRho*alphaF*a*betaZ,	sqrtRho*sbx*betaY,	sqrtRho*alphaS*a*betaZ	},
		{rho*alphaF*aSq,			0,					rho*alphaS*aSq,				0,	0,	rho*alphaS*aSq,				0,					rho*alphaF*aSq			},
	}

for j=1,8 do
	for k=1,8 do
		if not math.isfinite(evrw[j][k]) then
			error(tolua({
				tmp=tmp,
				betaZ=betaZ,
				--evrw=evrw,
			},{indent=true}))
		end
	end
end

	local evlw = {
		{0,-.5*alphaF*Cf*_1_aSq,.5*alphaS*Cs*sbx*betaY*_1_aSq,.5*alphaS*Cs*sbx*betaZ*_1_aSq,0,.5*alphaS*betaY*_1_sqrtRho*_1_a,.5*alphaS*betaZ*_1_sqrtRho*_1_a,.5*alphaF*_1_rho*_1_aSq},
		{0,0,-.5*betaZ,.5*betaY,0,-.5*sbx*betaZ*_1_sqrtRho,.5*sbx*betaY*_1_sqrtRho,0},
		{0,-.5*alphaS*Cs*_1_aSq,-.5*alphaF*Cf*sbx*betaY*_1_aSq,-.5*alphaF*Cf*sbx*betaZ*_1_aSq,0,-.5*alphaF*betaY*_1_sqrtRho*_1_a,-.5*alphaF*betaZ*_1_sqrtRho*_1_a,.5*alphaS*_1_rho*_1_aSq},
		{1,0,0,0,0,0,0,-_1_aSq},
		{0,0,0,0,1,0,0,0},
		{0,.5*alphaS*Cs*_1_aSq,.5*alphaF*Cf*sbx*betaY*_1_aSq,.5*alphaF*Cf*sbx*betaZ*_1_aSq,0,-.5*alphaF*betaY*_1_sqrtRho*_1_a,-.5*alphaF*betaZ*_1_sqrtRho*_1_a,.5*alphaS*_1_rho*_1_aSq},
		{0,0,.5*betaZ,-.5*betaY,0,-.5*sbx*betaZ*_1_sqrtRho,.5*sbx*betaY*_1_sqrtRho,0},
		{0,.5*alphaF*Cf*_1_aSq,-.5*alphaS*Cs*sbx*betaY*_1_aSq,-.5*alphaS*Cs*sbx*betaZ*_1_aSq,0,.5*alphaS*betaY*_1_sqrtRho*_1_a,.5*alphaS*betaZ*_1_sqrtRho*_1_a,.5*alphaF*_1_rho*_1_aSq}
	}

for j=1,8 do
	for k=1,8 do
		assertfinite(evlw[j][k])
	end
end

	-- for conservatives ordered: rho, mx, my, mz, bx, by, bz, ETotal 
	-- and primitives ordered: rho, vx, vy, vz, bx, by, bz, P
	
	local du_dw = {
		{1,0,0,0,0,0,0,0},
		{vx,rho,0,0,0,0,0,0},
		{vy,0,rho,0,0,0,0,0},
		{vz,0,0,rho,0,0,0,0},
		{0,0,0,0,1,0,0,0},
		{0,0,0,0,0,1,0,0},
		{0,0,0,0,0,0,1,0},
		{.5*vSq,rho*vx,rho*vy,rho*vz,bx,by,bz,1/(gamma-1)},
	}

	local dw_du = {
		{1,0,0,0,0,0,0,0},
		{-vx/rho,1/rho,0,0,0,0,0,0},
		{-vy/rho,0,1/rho,0,0,0,0,0},
		{-vz/rho,0,0,1/rho,0,0,0,0},
		{0,0,0,0,1,0,0,0},
		{0,0,0,0,0,1,0,0},
		{0,0,0,0,0,0,1,0},
		{.5*gamma_1*vSq,-gamma_1*vx,-gamma_1*vy,-gamma_1*vz,-gamma_1*bx,-gamma_1*by,-gamma_1*bz,gamma_1},
	}

	-- Ru = du/dw Rw
	local evru = sim.eigenvectors[i]
	for j=1,8 do
		for k=1,8 do
			local sum = 0
			for l=1,8 do
				sum = sum + du_dw[j][l] * evrw[l][k]
			end
			evru[j][k] = sum
		end
	end
	
	-- Lu = Lw dw/du
	local evlu = sim.eigenvectorsInverse[i]
	for j=1,8 do
		for k=1,8 do
			local sum = 0
			for l=1,8 do
				sum = sum + evlw[j][l] * dw_du[l][k]
			end
			evlu[j][k] = sum
		end
	end

--[[
	print()
	print('error of eigenbasis '..i)
	local errtotal = 0
	local everr = {}
	for j=1,#evrw do
		everr[j] = table()
		for k=1,#evrw do
			local sum = 0
			for l=1,#evrw do
				sum = sum + evlw[j][l] * evrw[l][k]
			end
			everr[j][k] = sum
			errtotal = errtotal + math.abs(sum - (j == k and 1 or 0))
		end
		print(everr[j]:concat', ')
	end
	print('error total',errtotal)
--]]
end

--[[ should I prevent unphysical densities and pressures here or in calcPrimFromCons?
-- turns out certain proportionality differences between P and ETotal makes it so even constraining things here
-- then recomputing them into conservative still creates conservative that give rise to nonphysical primitives
-- so I'll put the correction code in calcPrimFromCons
function MHD:postIterate(sim)
	if MHD.super.postIterate then MHD.super.postIterate(self, sim) end
	for i=1,sim.gridsize do
		local rho, vx, vy, vz, bx, by, bz, P = self:calcPrimFromCons(table.unpack(sim.qs[i]))
		rho = math.max(rho, 1e-7)
		P = math.max(P, 1e-7)
		fill(sim.qs[i], self:consFromPrim(rho, vx, vy, vz, bx, by, bz, P))
		-- hmm, what to do when the conservative variables don't want to stay conservative ...
		local rho, vx, vy, vz, bx, by, bz, P = self:calcPrimFromCons(table.unpack(sim.qs[i]))
		assert(P >= 1e-7)
	end
end
--]]

return MHD
