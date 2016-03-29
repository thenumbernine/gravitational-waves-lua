local class = require 'ext.class'
local Equation = require 'equation'
local mat33 = require 'mat33'

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
	local eInt = prim:_'eInt'
	local h = prim:_'h'
	local H = h * rho
	local EInt = eInt * rho
	local D = q:_(1)	-- rho * W
	local Sx = q:_(2)	-- rho h W^2 vx
	local tau = q:_(3)	-- rho h W^2 - P - D
	local W = D / rho	-- Lorentz factor
	SRHD1D:buildGraphInfos{
		{rho = rho},
		{vx = vx},
		{P = P},
		{EInt = EInt},
		{H = H},
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
	local x = sim.xs[i]
	--primitives:
	-- [[ Sod
	-- density 
	local rho = x < 0 and 1 or .125
	-- local velocity
	local vx = 0
	-- local pressure
	local P = 1 -- x < 0 and 1 or .1
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
	local roeWeightL = math.sqrt(sqrtAbsGL * rhoL * hL)
	local DL = qL[1]
	local WL = DL / rhoL
	local u0L = WL
	local u1L = vxL * u0L
	local roeVarsL = {u0L, u1L, PL / (rhoL * hL)}

	local sqrtAbsGR = 1
	local primR = sim.primitives[i]
	local qR = sim.qs[i]
	local rhoR = primR.rho
	local hR = primR.h
	local eIntR = primR.eInt
	local vxR = primR.vx
	local PR = pressureFunc(rhoR, eIntR)
	local roeWeightR = math.sqrt(sqrtAbsGR * rhoR * hR)
	local DR = qR[1]
	local WR = DR / rhoR
	local u0R = WR
	local u1R = vxR * u0R
	local roeVarsR = {u0R, u1R, PR / (rhoR * hR)}

-- [[ roe averaging.  i need to find a paper that describes how to average all variables, and not just these ...
	local normalize = 1 / (roeWeightL + roeWeightR)
	local roeVars = {}
	for i=1,3 do
		roeVars[i] = (roeVarsL[i] * roeWeightL + roeVarsR[i] * roeWeightR) * normalize
	end
	local W = roeVars[1]	-- = u0
	local u1 = roeVars[2]	-- = vx * W
	local vx = u1 / W
	-- eInt / h = (p / (rho h)) / (gamma - 1) 
	local P_over_rho_h = roeVars[3]

	-- so how do you figure 'h' for the eigenvector calculations?
	local h = (hL * roeWeightL + hR * roeWeightR) * normalize
	local rho = (rhoL * roeWeightL + rhoR * roeWeightR) * normalize
	local P = P_over_rho_h * rho * h
--]]
--[[ arithmetic averaging
	local rho = (rhoL + rhoR) / 2
	local eInt = (eIntL + eIntR) / 2
	local vx = (vxL + vxR) / 2
	local W = 1/math.sqrt(1-vx^2)
	local P = pressureFunc(rho, eInt)
	local h = 1 + eInt + P/rho
	local P_over_rho_h = P / (rho  * h)
--]]
	
	local vSq = vx^2 -- + vy^2 + vz^2
	local csSq = gamma * P_over_rho_h
	local cs = math.sqrt(csSq)

	-- Marti 1998 eqn 19
	-- also Marti & Muller 2008 eqn 68
	-- also Font 2008 eqn 106
	local lambda = sim.eigenvalues[i]
	local discr = math.sqrt((1 - vSq) * ((1 - vSq * csSq) - vx * vx * (1 - csSq)))
	-- slow
	lambda[1] = (vx * (1 - csSq) - cs * discr) / (1 - vSq * csSq)
	-- med
	lambda[2] = vx
	-- fast
	lambda[3] = (vx * (1 - csSq) + cs * discr) / (1 - vSq * csSq)

	-- Font 2008 eqn 113
	local VXMinus = (vx - lambda[3]) / (1 - vx * lambda[1])
	local VXPlus = (vx - lambda[1]) / (1 - vx * lambda[3])
	local AXTildeMinus = (1 - vx * vx) / (1 - vx * lambda[1])
	local AXTildePlus = (1 - vx * vx) / (1 - vx * lambda[3])

	--Font 2008 finally says that kappa = dp/deInt
	-- Marti & Muller 2008, Alcubierre 2008, and a few others I'm betting all have the eigenvalues, vectors, all the same ... but forget kappa
	local kappa = (gamma - 1) * rho
	
	-- Font 2008 eqn 112
	local CXMinus = vx - VXMinus
	local CXPlus  = vx - VXPlus
	local kappaTilde = kappa / rho
	local Kappa = kappaTilde / (kappaTilde - csSq)

	local evr = sim.eigenvectors[i]
	-- r- is the slow column
	evr[1][1] = 1
	evr[2][1] = h * W * CXMinus
	evr[3][1] = h * W * CXMinus - 1

	-- r0,1 is the velocity wave in the x direction
	evr[1][2] = Kappa / (h * W)
	evr[2][2] = vx
	evr[3][2] = 1 - Kappa / (h * W)

	-- r+ is the fast column
	evr[1][3] = 1
	evr[2][3] = h * W * CXPlus
	evr[3][3] = h * W * CXPlus - 1

	-- Font 2008 eqn 121
	local xi = 1 - vx * vx
	local Delta = h^3 * W * (Kappa - 1) * (CXPlus - CXMinus) * xi

	-- Font 2008 eqn 118
	-- Notice the left slow/fast vectors has a def wrt -+ and then uses 20x over +- ... I think the -+ might be a typo ... why would anything be defined in terms of -+?  usually -+ is used to denote the negative of a +- definition
	local evl = sim.eigenvectorsInverse[i]
	-- l- is the slow row
	evl[1][2] = -h^2 / Delta * ((1 - Kappa * AXTildePlus) + (2 * Kappa - 1) * VXPlus * (W * W * vx * xi - vx))
	evl[1][3] = -h^2 / Delta * ((1 - Kappa) * (-vx + VXPlus * (W * W * xi - 1)) - Kappa * W * W * VXPlus * xi)
	evl[1][1] = -h^2 / Delta * (h * W * VXPlus * xi + Delta / h^2 * evl[1][3])
	
	-- l0,1 is the velocity wave in the x direction
	evl[2][1] = W / (Kappa - 1) * (h - W)
	evl[2][2] = W / (Kappa - 1) * (W * vx)
	evl[2][3] = W / (Kappa - 1) * (-W)

	-- l+ is the fast row
	evl[3][2] = h^2 / Delta * ((1 - Kappa * AXTildeMinus) + (2 * Kappa - 1) * VXMinus * (W * W * vx * xi - vx))
	evl[3][3] = h^2 / Delta * ((1 - Kappa) * (-vx + VXMinus * (W * W * xi - 1)) - Kappa * W * W * VXMinus * xi)
	evl[3][1] = h^2 / Delta * (h * W * VXMinus * xi + Delta / h^2 * evl[3][3])

	-- define the flux matrix to see how accurate our 
	local dF_dW = {{},{},{}}
	dF_dW[1] = {}
	dF_dW[1][1] = W * vx
	dF_dW[2][1] = h * W * W * vx * vx - P / rho 
	dF_dW[3][1] = (h * W - 1) * W * vx
	dF_dW[1][2] = rho * W * (W * W * vx * vx + 1)
	dF_dW[2][2] = rho * h * W * W * (2 * W * W * vx * vx + 1) * vx + rho * h * W * W * vx
	dF_dW[3][2] = rho * h * W * W * (2 * W * W * vx * vx + 1) - rho * W * (W * W * vx * vx + 1)
	dF_dW[1][3] = 0
	dF_dW[2][3] = gamma * rho * W * W * vx * vx - (gamma - 1) * rho
	dF_dW[3][3] = gamma * rho * W * W * vx
	
	local dU_dW = {{},{},{}}
	dU_dW[1][1] = W
	dU_dW[2][1] = h * W * W * vx
	dU_dW[3][1] = W * (h * W - 1) * (gamma - 1)
	dU_dW[1][2] = rho * W * W * W * vx
	dU_dW[2][2] = rho * h * W * W * (2 * W * W * vx * vx + 1)
	dU_dW[3][2] = rho * W * W * W * vx * (2 * h * W - 1)
	dU_dW[1][3] = 0
	dU_dW[2][3] = gamma * rho * W * W * vx
	dU_dW[3][3] = rho * (gamma * W * W - (gamma - 1))
	local dW_dU = mat33.inv(dU_dW)

	for j=1,3 do
		for k=1,3 do
			local sum = 0
			for l=1,3 do
				sum = sum + dF_dW[j][l] * dW_dU[l][k]
			end
			sim.fluxMatrix[i][j][k] = sum
		end
	end
end

local function checknan(x, msg)
	assert(x==x, msg or "nan")
end

function SRHD1D:calcPrims(sim,i,prims, qs)
	local rho = prims.rho
	local eInt = prims.eInt

	-- use the new state variables to newton converge
	local D, Sx, tau = table.unpack(qs)
	
	-- code with Marti & Muller 2008:
	D = math.max(D, 1e-10)
	tau = math.max(tau, 1e-10)

	-- superluminal velocity will occur if |S| > tau + D + P
	-- i.e. if P < |S| - tau - D
	local absSx = math.abs(Sx)
	local velEpsilon = 1e-10
	local PMin = math.max(absSx - tau - D + absSx * velEpsilon, 1e-16)

	-- this is in the Marti & Muller 2008 code ...
	-- where do they get this from?
	local PMax = (gamma - 1) * tau
	-- why is there an upper bound again?
	--assert(PMax >= PMin, "pmax="..PMax.." < pmin="..PMin)
	PMax = math.max(PMax, PMin)

	-- start with a guess between PMin and PMax 
	local P = .5 * (PMin + PMax)
checknan(P)

	local maxiter = 1000
	for iter=1,maxiter do
		local vx = Sx / (tau + D + P)
-- PMin should prevent this from occuring
if math.abs(vx) > 1 then
	error("superluminal vx="..vx.." Sx="..Sx.." tau="..tau.." D="..D.." P="..P)
end
		-- in the case that this happens, the delta P can still be huge, but the min() keeps us in place ... 
		-- if we run up against the boundary then we should exit as well
checknan(vx)
		local vSq = vx^2	-- + vy^2 + vz^2
checknan(vSq)
		local W = 1 / math.sqrt(1 - vSq)
checknan(W, "W nan, vSq="..vSq.." vx="..vx)
		eInt = (tau + D * (1 - W) + P * (1 - W^2)) / (D * W)
checknan(eInt, "eInt nan, tau="..tau.." D="..D.." W="..W.." P="..P)
		rho = D / W
checknan(rho)

		-- for f = pressureFunc = (gamma - 1) * rho * eInt
		-- Marti & Muller 2008, eqn 61
		local f = pressureFunc(rho, eInt) - P
checknan(f)
		
		-- Alcubierre, 7.6.52
		--local csSq = gamma * P / (rho * P)	
		-- however the code with Marti & Muller uses a different equation for the cs^2 term of df/dp:
		local csSq = (gamma - 1) * (tau + D * (1 - W) + P) / (tau + D + P)
checknan(csSq)
		local df_dP = vSq * csSq - 1
checknan(df_dP)
		
		-- code with Marti & Muller says to keep things above PMin ...
		local newP = P - f / df_dP
		newP = math.max(newP, PMin)
checknan(newP)
		
		-- Marti & Muller 2008 code ...
		local PError = math.abs(1 - newP / P)
		P = newP
checknan(P)
		if math.abs(PError) < 1e-8 or iter == maxiter then
			prims.P = P
			prims.rho = rho
			prims.vx = Sx / (tau + D + P)
			if prims.vx < velEpsilon^2 then prims.vx = 0 end
			-- Marti & Muller 2008 recompute eInt here
			prims.eInt = P / (rho * (gamma - 1))
			prims.h = 1 + eInt + P/rho
			break
		end
	end
end

return SRHD1D
