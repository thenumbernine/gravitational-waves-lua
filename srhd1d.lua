local class = require 'ext.class'
local Equation = require 'equation'

local SRHD1D = class(Equation)
SRHD1D.name = 'SRHD 1D'
SRHD1D.numStates = 3

local gamma = 5/3

do
	local q = function(self,i) return self.qs[i] end
	local prim = function(self,i) return self.primitives[i] end
	local rho = prim:_'rho'	-- rest-mass rho
	local vx = prim:_'vx'
	local P = prim:_'P'
	local eInt = prim:_'eInt'	-- hmm this looks like densitized internal energy
	local h = prim:_'h'			-- and this looks densitized too ...
	local D = q:_(1)	-- rho * W
	local W = D / rho	-- Lorentz factor
	local Sx = q:_(2)	-- rho h W^2 vx
	local tau = q:_(3)	-- rho h W^2 - P - D
	SRHD1D:buildGraphInfos{
		{rho = rho},
		{vx = vx},
		{P = P},
--		{EInt = EInt},
--		{H = H},
		{D = D},
		{Sx = Sx},
		{tau = tau},
		{W = W},
		{['log eigenbasis error'] = function(self,i) return self.eigenbasisErrors and math.log(self.eigenbasisErrors[i]) end},
		{['log reconstruction error'] = function(self,i) return self.fluxMatrixErrors and math.log(self.fluxMatrixErrors[i]) end},
	}
end

local function pressureFunc(rho, eInt)
	return (gamma - 1) * rho * eInt
end
local function energyInternalFromPressure(rho, P)
	return P / ((gamma - 1) * rho)
end

function SRHD1D:initCell(sim,i)
	--primitives:
	-- [[ Sod
	-- density 
	local rho = sim.xs[i] < 0 and 1 or .1
	-- local velocity
	local vx = 0
	-- local pressure
	local P = 1
	--]]
	-- specific internal energy ... unspecified by the paper, but I'll use ideal gas heat capacity ratio
	local eInt = energyInternalFromPressure(rho, P)
	--aux variables:
	local vSq = vx^2 -- + vy^2 + vz^2
	local W = 1/math.sqrt(1 - vSq)
	-- specific enthalpy
	local h = 1 + eInt + P/rho
	--conservative variables:
	-- rest-mass density
	local D = rho * W
	-- momentum density
	local Sx = rho * h * W^2 * vx
	-- energy density
	local tau = rho * h * W^2 - P - D
	-- store primitive variables for later use
	print('assigning primitives to cell',i)
	sim.primitives[i] = {rho=rho, vx=vx, h=h, eInt=eInt, P=P}
	-- use conservative variables with the Roe scheme
	return {D, Sx, tau}
end

function SRHD1D:calcInterfaceEigenBasis(sim,i)
	-- but the Roe section says to use these variables:
	local sqrtAbsGL = 1	-- this will change once you start using a proper metric
	local primL = sim.primitives[i-1]
	local qL = sim.qs[i-1]
	local rhoL = assert(primL.rho, "failed to find rho for index "..(i-1))
	local hL = primL.h
	local eIntL = primL.eInt
	local vxL = primL.vx
	local PL = pressureFunc(rhoL, eIntL)
	local kL = math.sqrt(sqrtAbsGL * rhoL * hL)
	local DL = qL[1]
	local WL = DL / rhoL
	local u0L = WL
	local u1L = vxL * u0L
	local wL = {u0L, u1L, PL / (rhoL * hL)}

	local sqrtAbsGR = 1
	local primR = sim.primitives[i]
	local qR = sim.qs[i]
	local rhoR = primR.rho
	local hR = primR.h
	local eIntR = primR.eInt
	local vxR = primR.vx
	local PR = pressureFunc(rhoR, eIntR)
	local kR = math.sqrt(sqrtAbsGR * rhoR * hR)
	local DR = qR[1]
	local WR = DR / rhoR
	local u0R = WR
	local u1R = vxR * u0R
	local wR = {u0R, u1R, PR / (rhoR * hR)}

--[[ roe averaging.  i need to find a paper that describes how to average all variables, and not just these ...
	local w = {}
	for i=1,3 do
		w[i] = (wL[i] * kL + wR[i] * kR) / (kL + kR)
	end
	local W = w[1]	-- = u0
	local u1 = w[2]	-- = vx * W
	-- TODO now recover the prim values from the Roe averaged 4-velocities
	local vx = u1 / W
--]]
-- [[ arithmetic averaging
	local W = (WL + WR) / 2
	local rho = (rhoL + rhoR) / 2
	local eInt = (eIntL + eIntR) / 2
	local vx = (vxL + vxR) / 2
--]]
	
	local vSq = vx^2 -- + vy^2 + vz^2
	-- w[3] = P / (rho * h) <=> h = P(rho,eInt) / (rho * w[3])
	-- so how do you recover eInt, rho, and h from w[3]?
	local P = pressureFunc(rho, eInt)
	local h = 1 + eInt + P/rho
	local csSq = gamma * P / (rho * h)
	local cs = math.sqrt(csSq)

	local discr = math.sqrt((1 - vSq) * (1 - vx * vx - (vSq - vx * vx) * csSq))
	
	local lambda = sim.eigenvalues[i]
	-- slow
	lambda[1] = (vx * (1 - csSq) - cs * discr) / (1 - vSq * csSq)
	-- med
	lambda[2] = vx
	-- fast
	lambda[3] = (vx * (1 - csSq) + cs * discr) / (1 - vSq * csSq)

	local Aminus = (1 - vx * vx) / (1 - vx * lambda[1])
	local Aplus  = (1 - vx * vx) / (1 - vx * lambda[3])
	
	--local kappa = ...
	--local kappaTilde = kappa / rho
	--local Kappa = kappaTilde / (kappaTilde - csSq)
	-- "for an ideal gas equation of state, Kappa = h, so Kappa > 1, so Delta != 0 for |vx| < 1"
	-- but doesn't vx have no limit? doesn't the lorentz transform map [0,inf) to [0,1] ?
	local Kappa = h

	local evr = sim.eigenvectors[i]
	-- r- is the slow column
	evr[1][1] = 1
	evr[2][1] = h * W * Aminus * lambda[1]
	-- h * W * vy
	-- h * W * vz
	evr[3][1] = h * W * Aminus - 1

	-- r0,1 is the velocity wave in the x direction
	evr[1][2] = Kappa / (h * W)
	evr[2][2] = vx
	-- vy
	-- vz
	evr[3][2] = 1 - Kappa / (h * W)

	-- r0,2 = (W * vy, 2 * h * W^2 * vx * vy, h * (1 + 2 * W^2 * vy^2),    2 * h * W^2 * vy * vz, (2 * h * W - 1) * W * vy)
	-- r0,3 = (W * vz, 2 * h * W^2 * vx * vz, h * (1 + 2 * W^2 * vy * vz), 2 * h * W^2 * vz^2,    (2 * h * W - 1) * W * vz)

	-- r+ is the fast column
	evr[1][3] = 1
	evr[2][3] = h * W * Aplus * lambda[3]
	-- h * W * vy
	-- h * W * vz
	evr[3][3] = h * W * Aplus - 1

	local Delta = h^3 * W * (Kappa - 1) * (1 - vx^2) * (Aplus * lambda[3] - Aminus * lambda[1])
	local evl = sim.eigenvectorsInverse[i]
	-- Notice the left slow/fast vectors has a def wrt -+ and then uses 20x over +- ... I think the -+ might be a typo ... why would anything be defined in terms of -+?  usually -+ is used to denote the negative of a +- definition
	-- if not then just swap the fast and the slow waves
	-- l- is the slow row
	evl[1][1] = -h^2 / Delta * (h * W * Aminus * (vx - lambda[1]) - vx - W^2 * (vSq - vx^2) * (2 * Kappa - 1) * (vx - Aminus * lambda[1]) + Kappa * Aminus * lambda[1])
	evl[1][2] = -h^2 / Delta * (1 + W^2 * (vSq - vx^2) * (2 * Kappa - 1) * (1 - Aminus) - Kappa * Aminus)
	-- vy (-h^2 / Delta * W^2 (2 Kappa - 1) Aminus (vx - lambda[1]))
	-- vz (-h^2 / Delta * W^2 (2 Kappa - 1) Aminus (vx - lambda[1]))
	evl[1][3] = -h^2 / Delta * (-vx - W^2 * (vSq - vx^2) * (2 * Kappa - 1) * (vx - Aminus * lambda[1]) + Kappa * Aminus * lambda[1])
	
	-- l0,1 is the velocity wave in the x direction
	evl[2][1] = W / (Kappa - 1) * (h - W)
	evl[2][2] = W / (Kappa - 1) * W * vx
	-- W / (Kappa - 1) * W * vy
	-- W / (Kappa - 1) * W * vz
	evl[2][3] = W / (Kappa - 1) * (-W)

	-- l0,2 = (-vy, vx vy, 1 - vx^2, 0, -vy) / (h(1 - vx^2))
	-- l0,3 = (-vz, vx vz, 0, 1 - vx^2, -vz) / (h(1 - vx^2))

	-- l+ is the fast row
	evl[3][1] = h^2 / Delta * (h * W * Aplus * (vx - lambda[3]) - vx - W^2 * (vSq - vx^2) * (2 * Kappa - 1) * (vx - Aplus * lambda[3]) + Kappa * Aplus * lambda[3])
	evl[3][2] = h^2 / Delta * (1 + W^2 * (vSq - vx^2) * (2 * Kappa - 1) * (1 - Aplus) - Kappa * Aplus)
	-- vy (h^2 / Delta * W^2 (2 Kappa - 1) Aplus (vx - lambda[3]))
	-- vz (h^2 / Delta * W^2 (2 Kappa - 1) Aplus (vx - lambda[3]))
	evl[3][3] = h^2 / Delta * (-vx - W^2 * (vSq - vx^2) * (2 * Kappa - 1) * (vx - Aplus * lambda[3]) + Kappa * Aplus * lambda[3])

	-- has to be defined, even if it's not used ...
	-- TODO fix that
	local F = sim.fluxMatrix[i]
	for i=1,3 do
		for j=1,3 do
			F[i][j] = 0
		end
	end
end

function SRHD1D:calcPrims(sim,i,prims, qs)
	-- start with last iteration's pressure prim as our initial guess
	local PStar = pressureFunc(prims.rho, prims.eInt)

	-- use the new state variables to newton converge
	local D, Sx, tau = table.unpack(qs)
	for iter=1,100 do
		local vxStar = Sx / (tau + D + PStar)
		local vSqStar = vxStar^2	-- + vyStar^2 + vzStar^2
		local WStar = 1 / math.sqrt(1 - vSqStar)
		local eIntStar = (tau + D * (1 - WStar) + PStar * (1 - WStar^2)) / (D * WStar)
		local rhoStar = D / WStar
		local csSqStar = gamma * PStar / (rhoStar * PStar)	-- Alcubierre, 7.6.52
		-- for f = pressureFunc = (gamma - 1) * rho * eInt
		local f = pressureFunc(rhoStar, eIntStar) - PStar
		local df_dPStar = vSqStar * csSqStar - 1
		-- hmm, seems to be diverging ...
		local deltaPStar = -f / df_dPStar
		print(i,'deltaPStar',deltaPStar)
		if math.abs(deltaPStar) < 1e-7 then
			prims.rho = rhoStar
			prims.vx = vxStar
			prims.eInt = eIntStar
			prims.h = 1 + eIntStar + PStar/rhoStar
			prims.P = PStar
			break
		end
		PStar = PStar + deltaPStar
		if iter==100 then error("failed to converge") end
	end
end

return SRHD1D
