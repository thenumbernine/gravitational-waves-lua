local class = require 'ext.class'
local Equation = require 'equation'
local mat33 = require 'mat33'

local function checknan(x, msg)
	assert(x==x, msg or 'nan')
end

local SRHD1D = class(Equation)
SRHD1D.name = 'SRHD 1D'
SRHD1D.numStates = 3

-- pressure functions for ideal gas
function SRHD1D:calcP(rho, eInt)
	return (self.gamma - 1) * rho * eInt
end
-- chi in most papers
function SRHD1D:calc_dP_drho(rho, eInt)
	return (self.gamma - 1) * eInt
end
-- kappa in most papers
function SRHD1D:calc_dP_deInt(rho, eInt)
	return (self.gamma - 1) * rho
end
function SRHD1D:calc_eInt_from_P(rho, P)
	return P / ((self.gamma - 1) * rho)
end

function SRHD1D:init(...)
	assert(not SRHD1D.super.init)

	self.gamma = 5/3
	
	local q = function(self,i) return self.qs[i] end
	local prim = function(self,i) return self.ws[i] end
	local rho = prim:_(1)
	local vx = prim:_(2)
	local eInt = prim:_(3)
	
	-- a function to return a function computed via arithmetic operations on functions ...
	-- the problem is that calcP will multiply self.gamma once, which will create a lambda closure with the number value
	-- so if self.gamma is changed later, it won't reflect here
	-- solution: make gamma a function (for the time begin) while the closure is built
	local originalGamma = self.gamma
	self.gamma = function(self,i) return originalGamma end
	local P = self:calcP(rho, eInt)
	-- ...and restore it
	self.gamma = originalGamma
	
	local h = 1 + eInt + P/rho	
	local eInt = eInt * rho
	local D = q:_(1)	-- rho * W
	local Sx = q:_(2)	-- rho h W^2 vx
	local tau = q:_(3)	-- rho h W^2 - P - D
	local W = D / rho	-- Lorentz factor
	self:buildGraphInfos{
		{rho = rho},
		{vx = vx},
		{eInt = eInt},
		{P = P},
		{h = h},
		{D = D},
		{Sx = Sx},
		{tau = tau},
		{W = W},
		{['log eigenbasis error'] = function(self,i) return self.eigenbasisErrors and math.log(self.eigenbasisErrors[i]) end},
		{['log reconstruction error'] = function(self,i) return self.fluxMatrixErrors and math.log(self.fluxMatrixErrors[i]) end},
		{['primitive reconstruction error'] = function(self,i) return self.primitiveReconstructionErrors and self.primitiveReconstructionErrors[i] and math.log(self.primitiveReconstructionErrors[i]) end},
	}
end

function SRHD1D:consFromPrims(rho, vx, eInt)
	local vSq = vx*vx
	local WSq = 1/(1-vSq)
	local W = math.sqrt(WSq)
	local P = self:calcP(rho, eInt)
	local h = 1 + eInt + P/rho
	-- rest-mass density
	local D = rho * W
	-- momentum density
	local Sx = rho * h * WSq * vx
	-- energy density
	local tau = rho * h * WSq - P - D
	return D, Sx, tau
end

function SRHD1D:initCell(sim,i)
	local x = sim.xs[i]
	--primitives:
	--[[ Sod
	self.gamma = 7/5
	local rho = x < 0 and 1 or .125
	local vx = 0
	local P = x < 0 and 1 or .1
	--]]
	--[[ Marti & Muller 2008 rppm's Schneider et al
	self.gamma = 5/3
	local rho = x < 0 and 10 or 1
	local vx = 0
	local P = (self.gamma-1) * rho * (x < 0 and 2 or 1e-6)
	--]]
	--[[ Marti & Muller 2008 rppm relativistic blast wave
	self.gamma = 5/3
	local rho = 1
	local vx = 0
	local P = x < 0 and 1e+3 or 1e-2
	--]]
	--[[ Marti & Muller 2008 rppm relativistic shock reflection
	-- not working with my sim...
	self.gamma = 4/3
	local rho = 1
	local vx = .99999
	local P = (self.gamma - 1) * rho * (1e-7 / math.sqrt(1 - vx*vx))
	--]]
	--[[ Marti & Muller 2008 rppm relativistic blast wave interaction
	-- gets nans in the left-eigenvectors when the shockwaves collide
	-- under both analytical and numerical calculuation
	self.gamma = 7/5
	local lhs = .9 * sim.domain.xmin + .1 * sim.domain.xmax
	local rhs = .1 * sim.domain.xmin + .9 * sim.domain.xmax
	local rho = 1
	local vx = 0
	local P = x < lhs and 1e+3 or (x > rhs and 1e+2 or 1e-2)
	--]]
	-- [[ relativistic blast wave test problem 1, Marti & Muller 2008, table 5
	self.gamma = 5/3
	local rho = x < 0 and 10 or 1
	local vx = 0
	local P = x < .01 and 40/3 or 1e-3	-- paper says 0 for rhs but I'm putting .01 and hoping for a typo 
	--]]
	--[[ relativistic blast wave test problem 2, Marti & Muller 2008, table 5
	self.gamma = 5/3
	local rho = 1
	local vx = 0
	local P = x < 0 and 1e+3 or 1e-2
	--]]
	local eInt = self:calc_eInt_from_P(rho, P)
	local vSq = vx*vx -- + vy*vy + vz*vz
	local W = 1/math.sqrt(1 - vSq)
	local h = 1 + eInt + P/rho
	sim.ws[i] = {rho, vx, eInt}
	return {self:consFromPrims(rho, vx, eInt)}
end

function SRHD1D:calcInterfaceEigenBasis(sim,i)
	-- but the Roe section says to use these variables:
	local sqrtAbsGL = 1	-- this will change once you start using a proper metric
	local DL, _, _ = table.unpack(sim.qs[i-1])
	local rhoL, vxL, eIntL = table.unpack(sim.ws[i-1])
	local WL = DL / rhoL
	local PL = self:calcP(rhoL, eIntL)
	local hL = 1 + eIntL + PL/rhoL
	local roeWeightL = math.sqrt(sqrtAbsGL * rhoL * hL)

	local sqrtAbsGR = 1
	local DR, _, _ = table.unpack(sim.qs[i])
	local rhoR, vxR, eIntR = table.unpack(sim.ws[i])
	local WR = DR / rhoR
	local PR = self:calcP(rhoR, eIntR)
	local hR = 1 + eIntR + PR/rhoR
	local roeWeightR = math.sqrt(sqrtAbsGR * rhoR * hR)

--[[ roe averaging.  i need to find a paper that describes how to average all variables, and not just these ...
	local normalize = 1 / (roeWeightL + roeWeightR)
	local u1L = vxL * WL 
	local roeVarsR = {WR, u1R, PR / (rhoR * hR)}
	local u1R = vxR * WR 
	local roeVarsL = {WL, u1L, PL / (rhoL * hL)}
	local roeVars = {}
	for i=1,3 do
		roeVars[i] = (roeVarsL[i] * roeWeightL + roeVarsR[i] * roeWeightR) * normalize
	end
	local W = roeVars[1]	-- = u0
	local u1 = roeVars[2]	-- = vx * W
	local vx = u1 / W
	-- eInt / h = (p / (rho h)) / (self.gamma - 1) 
	local P_over_rho_h = roeVars[3]

	-- so how do you figure 'h' for the eigenvector calculations?
	local h = (hL * roeWeightL + hR * roeWeightR) * normalize
	local rho = (rhoL * roeWeightL + rhoR * roeWeightR) * normalize
	local P = P_over_rho_h * rho * h
--]]
-- [[ arithmetic averaging
	local rho = (rhoL + rhoR) / 2
	local eInt = (eIntL + eIntR) / 2
	local vx = (vxL + vxR) / 2
	local vSq = vx * vx
	local W = 1/math.sqrt(1-vSq)
	local P = self:calcP(rho, eInt)
	local h = 1 + eInt + P/rho
	local P_over_rho_h = P / (rho  * h)
--]]

	local hW = h * W
	local W2 = W * W
	local W3 = W2 * W
	local W4 = W2 * W2
	local hSq = h * h

	local vxSq = vx * vx	-- this is where it's just the flux direction squared
	local csSq = self.gamma * P_over_rho_h
	local cs = math.sqrt(csSq)

	-- Marti 1998 eqn 19
	-- also Marti & Muller 2008 eqn 68
	-- also Font 2008 eqn 106
	local lambda = sim.eigenvalues[i]
	local discr = math.sqrt((1 - vSq) * ((1 - vSq * csSq) - vxSq * (1 - csSq)))
	-- slow
	lambda[1] = (vx * (1 - csSq) - cs * discr) / (1 - vSq * csSq)
	-- med
	lambda[2] = vx
	-- fast
	lambda[3] = (vx * (1 - csSq) + cs * discr) / (1 - vSq * csSq)
	
	--Font 2008 finally says that kappa = dp/deInt
	-- Marti & Muller 2008, Alcubierre 2008, and a few others I'm betting all have the eigenvalues, vectors, all the same ... but forget kappa
	local kappa = (self.gamma - 1) * rho
	
	-- Font 2008 eqn 112, also Marti & Muller 2008 eqn 74
	local kappaTilde = kappa / rho
	local Kappa = kappaTilde / (kappaTilde - csSq)

-- [=[ Marti, Muller 2008
	local AMinus= (1 - vxSq) / (1 - vx * lambda[1])
	local APlus = (1 - vxSq) / (1 - vx * lambda[3])
	local evr = sim.eigenvectors[i]
	evr[1][1] = 1
	evr[2][1] = hW * AMinus * lambda[1]
	evr[3][1] = hW * AMinus - 1
	evr[1][2] = Kappa/hW
	evr[2][2] = vx
	evr[3][2] = 1 - Kappa/hW
	evr[1][3] = 1
	evr[2][3] = hW * APlus * lambda[3]
	evr[3][3] = hW * APlus - 1

	local Delta = h*h*h * W * (Kappa - 1) * (1 - vxSq) * (APlus * lambda[3] - AMinus * lambda[1])
	local evl = sim.eigenvectorsInverse[i]
	evl[1][1] = hSq / Delta * ( hW * APlus * (vx - lambda[3]) - vx - W2 * (vSq - vxSq) * (2 * Kappa - 1) * (vx - APlus * lambda[3]) + Kappa * APlus * lambda[3] )
	evl[1][2] = hSq / Delta * (1 + W2 * (vSq - vxSq) * (2 * Kappa - 1) * (1 - APlus) - Kappa * APlus )
	evl[1][3] = hSq / Delta * (-vx - W2 * (vSq - vxSq) * (2 * Kappa - 1) * (vx - APlus * lambda[3]) + Kappa * APlus * lambda[3] )
	evl[2][1] = W / (Kappa - 1) * (h - W)
	evl[2][2] = W / (Kappa - 1) * W * vx
	evl[2][3] = W / (Kappa - 1) * -W
	evl[3][1] = -hSq / Delta * (hW * AMinus * (vx - lambda[1]) - vx - W2 * (vSq - vxSq) * (2 * Kappa - 1) * (vx - AMinus * lambda[1]) + Kappa * AMinus * lambda[1] )
	evl[3][2] = -hSq / Delta * (1 + W2 * (vSq - vxSq) * (2 * Kappa - 1) * (1 - AMinus) - Kappa * AMinus )
	evl[3][3] = -hSq / Delta * (-vx - W2 * (vSq - vxSq) * (2 * Kappa - 1) * (vx - AMinus * lambda[1]) + Kappa * AMinus * lambda[1] )
--]=]
--[=[ Font 2008
	-- Font 2008 eqn 113
	local VXMinus = (vx - lambda[3]) / (1 - vx * lambda[1])
	local VXPlus = (vx - lambda[1]) / (1 - vx * lambda[3])
	local AXTildeMinus = (1 - vxSq) / (1 - vx * lambda[1])
	local AXTildePlus = (1 - vxSq) / (1 - vx * lambda[3])
	
	-- Font 2008 eqn 112
	local CXMinus = vx - VXMinus
	local CXPlus  = vx - VXPlus

	local evr = sim.eigenvectors[i]
	-- r- is the slow column
	evr[1][1] = 1
	evr[2][1] = hW * CXMinus
	evr[3][1] = hW * CXMinus - 1

	-- r0,1 is the velocity wave in the x direction
	evr[1][2] = Kappa / hW
	evr[2][2] = vx
	evr[3][2] = 1 - Kappa / hW

	-- r+ is the fast column
	evr[1][3] = 1
	evr[2][3] = hW * CXPlus
	evr[3][3] = hW * CXPlus - 1

	-- Font 2008 eqn 121
	local xi = 1 - vxSq
	local Delta = h * h * h * W * (Kappa - 1) * (CXPlus - CXMinus) * xi

	--[[
	-- Font 2008 eqn 118
	-- Notice the left slow/fast vectors has a def wrt -+ and then uses 20x over +- ... I think the -+ might be a typo ... why would anything be defined in terms of -+?  usually -+ is used to denote the negative of a +- definition
	local evl = sim.eigenvectorsInverse[i]
	-- l- is the slow row
	evl[1][2] = -hSq / Delta * ((1 - Kappa * AXTildePlus) + (2 * Kappa - 1) * VXPlus * (W2 * vx * xi - vx))
	evl[1][3] = -hSq / Delta * ((1 - Kappa) * (-vx + VXPlus * (W2 * xi - 1)) - Kappa * W2 * VXPlus * xi)
	evl[1][1] = -hSq / Delta * (h * W * VXPlus * xi + Delta / hSq * evl[1][3])
	
	-- l0,1 is the velocity wave in the x direction
	evl[2][1] = W / (Kappa - 1) * (h - W)
	evl[2][2] = W / (Kappa - 1) * (W * vx)
	evl[2][3] = W / (Kappa - 1) * (-W)

	-- l+ is the fast row
	evl[3][2] = hSq / Delta * ((1 - Kappa * AXTildeMinus) + (2 * Kappa - 1) * VXMinus * (W2 * vx * xi - vx))
	evl[3][3] = hSq / Delta * ((1 - Kappa) * (-vx + VXMinus * (W2 * xi - 1)) - Kappa * W2 * VXMinus * xi)
	evl[3][1] = hSq / Delta * (h * W * VXMinus * xi + Delta / hSq * evl[3][3])
	--]]
	--[[ or just do it numerically
	local evl = sim.eigenvectorsInverse[i]
	sim.eigenvectorsInverse[i] = mat33.inv(evr)
	--]]
--]=]

for j=1,3 do
	checknan(lambda[j])
	for k=1,3 do
		checknan(evl[j][k])
		checknan(evr[j][k])
	end
end

	-- define the flux matrix to see how accurate our 
	-- this is coming out different than the matrix reconstructed from the eigensystem
	local dF_dw = {{},{},{}}
	dF_dw[1][1] = W * vx
	dF_dw[2][1] = (-P/rho + h * W2 - h)
	dF_dw[3][1] = W * (hW - 1) * vx
	dF_dw[1][2] = rho * W3
	dF_dw[2][2] = 2 * rho * W4 * h * vx
	dF_dw[3][2] = rho * W * (1 + W4 - 3*W2 + 2*W3*h - h*W - W4*vx*vx*vx*vx)
	dF_dw[1][3] = 0
	dF_dw[2][3] = rho * (1 - 2 * self.gamma + self.gamma * W2)
	dF_dw[3][3] = self.gamma * rho * W2 * vx
	
	local dU_dw = {{},{},{}}
	dU_dw[1][1] = W
	dU_dw[2][1] = h * W2 * vx
	dU_dw[3][1] = (-P/rho - W + h * W2)
	dU_dw[1][2] = rho * W3 * vx
	dU_dw[2][2] = rho * h * W2 * (2 * W2 - 1)
	dU_dw[3][2] = rho * W3 * (2 * hW - 1) * vx
	dU_dw[1][3] = 0
	dU_dw[2][3] = self.gamma * rho * W2 * vx
	dU_dw[3][3] = rho * (1 - self.gamma + W2 * self.gamma)
	local dw_dU = mat33.inv(dU_dw)

	for j=1,3 do
		for k=1,3 do
			local sum = 0
			for l=1,3 do
				sum = sum + dF_dw[j][l] * dw_dU[l][k]
			end
			sim.fluxMatrix[i][j][k] = sum
		end
	end
end

SRHD1D.solvePrimMaxIter = 1000
SRHD1D.solvePrimStopEpsilon = 1e-7

-- this is my attempt based on the recover pressure method described in Alcubierre, Baumgarte & Shapiro, Marti & Muller, Font, and generally everywhere
function SRHD1D:calcPrimsByPressure(sim, i ,prims, qs)
	local rho, vx, eInt = table.unpack(prims)

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
	local PMax = (self.gamma - 1) * tau
	-- why is there an upper bound again?
	--assert(PMax >= PMin, 'pmax='..PMax..' < pmin='..PMin)
	PMax = math.max(PMax, PMin)

	-- start with a guess between PMin and PMax 
	local P = .5 * (PMin + PMax)
checknan(P)

	for iter=1,self.solvePrimMaxIter do
		vx = Sx / (tau + D + P)
-- PMin should prevent this from occuring
if math.abs(vx) > 1 then
	error('superluminal vx='..vx..' Sx='..Sx..' tau='..tau..' D='..D..' P='..P)
end
		-- in the case that this happens, the delta P can still be huge, but the min() keeps us in place ... 
		-- if we run up against the boundary then we should exit as well
checknan(vx)
		local vSq = vx*vx	-- + vy*vy + vz*vz
checknan(vSq)
		local W = 1 / math.sqrt(1 - vSq)
checknan(W, 'W nan, vSq='..vSq..' vx='..vx)
		eInt = (tau + D * (1 - W) + P * (1 - W*W)) / (D * W)
checknan(eInt, 'eInt nan, tau='..tau..' D='..D..' W='..W..' P='..P)
		rho = D / W
checknan(rho)

		-- for f = calcP = (self.gamma - 1) * rho * eInt
		-- Marti & Muller 2008, eqn 61
		local f = self:calcP(rho, eInt) - P
checknan(f)
		
		-- Alcubierre, 7.6.52
		--local csSq = self.gamma * P / (rho * P)	
		-- however the code with Marti & Muller uses a different equation for the cs^2 term of df/dp:
		local csSq = (self.gamma - 1) * (tau + D * (1 - W) + P) / (tau + D + P)
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
		if math.abs(PError) < self.solvePrimStopEpsilon then
			-- one last update ...
			vx = Sx / (tau + D + P)
			W = 1 / math.sqrt(1 - vx*vx)
			rho = D / W
			eInt = self:calc_eInt_from_P(rho, P)
			return {rho, vx, eInt}, nil, iter
		end
	end
	return nil, "didn't converge", self.solvePrimMaxIter
end

--[[
here's the method described in Anton & Zanotti 2006, used by Mara:
Mara converges Z=rho h W^2 and W, while the paper says to converge rho and W
either way, if you converge W then you're calculating things based on a W apart from your conservative W
unless you replaced the conservative W with the converged W
--]]
function SRHD1D:calcPrimsByZAndW(sim, i, prims, qs)
	local D, Sx, tau = table.unpack(qs)
	local rho, vx, eInt = table.unpack(prims)
	local vSq = vx*vx
	local W = 1/math.sqrt(1-vSq)	-- W is the other
	local P = self:calcP(rho,eInt)
	local h = 1 + eInt + P/rho
	local Z = rho * h * W*W -- Z = rho h W^2 is one var we converge
	for iter=1,self.solvePrimMaxIter do
		rho = D / W
		h = Z / (D * W)		-- h-1 = gamma eInt <=> eInt = (h-1) / gamma  <=> (h-1) = eInt + P/rho
		eInt = (h - 1) / self.gamma	-- isn't this only true for ideal gasses?
		P = self:calcP(rho,eInt)
	
		-- f = [-S^2 + Z^2 (W^2 - 1) / W^2, -tau + Z^2 - P - D]
		local f = {
			-Sx*Sx  + Z*Z * (W*W - 1) / (W*W),
			-tau +  Z - P - D,
		}
		local df_dx = {
			{2*Z * (W*W - 1) * Z*Z*Z / (W*W*Z*Z*Z), 2*Z*Z / (W*W*W)},
			{1 - (self.gamma - 1) / (self.gamma * W*W), (2*Z - D*W)/(W*W*W) * (self.gamma-1)/self.gamma},
		}
		local det = df_dx[1][1] * df_dx[2][2] - df_dx[2][1] * df_dx[1][2]
		local dx_df = {
			{df_dx[2][2] / det, -df_dx[2][1] / det},
			{-df_dx[1][2] / det, df_dx[1][1] / det},
		}
		local dZ = -(dx_df[1][1] * f[1] + dx_df[1][2] * f[2])
		local dW = -(dx_df[2][1] * f[1] + dx_df[2][2] * f[2])
		local newZ = Z + dZ
		local newW = W + dW
		newZ = math.abs(newZ)
		Z = newZ < 1e+20 and newZ or Z	-- restore old Z if we exceed 1e+20
		W = math.clamp(newW,1,1e+12)
		local err = math.abs(dZ/Z) + math.abs(dW/W)
		if err < self.solvePrimStopEpsilon then
			rho = D / W
			h = Z / (D * W)		-- h-1 = gamma eInt <=> eInt = (h-1) / gamma  <=> (h-1) = eInt + P/rho
			eInt = (h - 1) / self.gamma	-- this is ideal gas law only, right?
			P = self:calcP(rho, eInt)
			vx = math.clamp(Sx / (tau + D + P), -1, 1)
			return {rho, vx, eInt}, nil, iter
		end
	end
	return nil, "didn't converge", self.solvePrimMaxIter
end

--[[
how about 3-var newton descent?   D,Sx,tau -> rho,vx,eInt
newton update: w = w - dw/dU * U
U = the zeros of the cons eqns
-D' + rho * W
-Sx' + rho * h * W^2 * vx
-tau' + rho * h * W^2 - P - D
...for D', Sx', tau' fixed
...and all else calculated based on rho,vx,eInt
--]]
function SRHD1D:calcPrimsByPrims(sim, i, prims, qs)
	--fixed cons to converge to
	local D, Sx, tau = table.unpack(qs)
checknan(D)
checknan(Sx)
checknan(tau)
	--dynamic prims to converge
	local rho, vx, eInt = table.unpack(prims)
checknan(rho)
checknan(vx)
checknan(eInt)
	--print('cell',i)
	for iter=1,self.solvePrimMaxIter do
		--print('rho='..rho..' vx='..vx..' eInt='..eInt)
		local P = self:calcP(rho, eInt)
		local dP_drho = self:calc_dP_drho(rho, eInt)
		local dP_deInt = self:calc_dP_deInt(rho, eInt)
checknan(P)
		local h = 1 + eInt + P/rho
checknan(h)
		local W2 = 1/(1-vx*vx)
		local W = math.sqrt(W2)
checknan(W, 'W nan with vx='..vx)
		local W3 = W*W2
		local dU_dw = {
			{W, rho * W3 * vx, 0},
			{W2 * vx * (1 + eInt + dP_drho), rho * h * W2 * (2 * W2 - 1), W2 * vx * (dP_deInt + rho)},
			{-dP_drho - W + W2 * (1 + eInt + dP_drho), rho * W3 * vx * (2 * W * h - 1), -dP_deInt + dP_deInt * W2 + rho * W2}
		}
for j=1,3 do
	for k=1,3 do
		checknan(dU_dw[j][k], 'nan for dU_dw['..j..']['..k..']')
	end
end
		local det_dU_dw = mat33.det(dU_dw)
checknan(det_dU_dw, 'nan for derivative of:\n'..require'matrix'(dU_dw))
		local dw_dU = mat33.inv(dU_dw, det_dU_dw)
for j=1,3 do
	for k=1,3 do
		checknan(dw_dU[j][k])
	end
end
		local dU = {
			-D + rho * W,
			-Sx + rho * h * W2 * vx,
			-tau + rho * h * W2 - P - rho * W,
		}
checknan(dU[1], 'dU[1] nan with D='..D..' rho='..rho..' W='..W)
for j=2,3 do
	checknan(dU[j], 'nan for dU['..j..']')
end
		drho = dw_dU[1][1] * dU[1] + dw_dU[1][2] * dU[2] + dw_dU[1][3] * dU[3]
checknan(drho)
		dvx = dw_dU[2][1] * dU[1] + dw_dU[2][2] * dU[2] + dw_dU[2][3] * dU[3]
checknan(dvx)
		deInt = dw_dU[3][1] * dU[1] + dw_dU[3][2] * dU[2] + dw_dU[3][3] * dU[3]
checknan(deInt)
		local alpha = 1
		local new_rho = rho - alpha * drho
		local new_vx = vx - alpha * dvx
		local new_eInt = eInt - alpha * deInt
		-- [[ enforce boundaries
		local vx_epsilon = 1e-7		-- can't reach the speed of light or W will get a NaN
		new_rho = math.clamp(new_rho,1e-7,1e+20)
		new_vx = math.clamp(new_vx,-1+vx_epsilon,1-vx_epsilon)
		new_eInt = math.clamp(new_eInt,0,1e+20)
		--]]
		local err = math.abs(new_rho-rho) + math.abs(new_vx-vx) + math.abs(new_eInt-eInt)
		rho = new_rho
		vx = new_vx
		eInt = new_eInt
		local cons = table{self:consFromPrims(rho,vx,eInt)}
		local consErr = table{math.abs(cons[1]-qs[1]), math.abs(cons[2]-qs[2]), math.abs(cons[3]-qs[3])}
		--print('err='..err)
		if err < self.solvePrimStopEpsilon then
			return {rho, vx, eInt}, nil, iter
		end
	end
	return nil, "didn't converge", self.solvePrimMaxIter
end

return SRHD1D
