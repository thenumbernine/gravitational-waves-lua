local class = require 'ext.class'
local Equation = require 'equation'
local mat33 = require 'mat33'

local SRHD1D = class(Equation)
SRHD1D.name = 'SRHD 1D'
SRHD1D.numStates = 3


SRHD1D.solvePrimMaxIter = 1000
SRHD1D.solvePrimErrorEpsilon = 1e-7

-- used by pressure solver
-- velocity epsilon is how close we can get to the speed of light
-- set ylabel "Lorentz factor"; set xlabel "velocity epsilon -log10"; set log xy; plot [1:10] 1/sqrt(1-(1-10**(-x))**2);
--SRHD1D.velEpsilon = 1e-5	-- <=> handles up to W = 500
--SRHD1D.velEpsilon = 1e-6	-- <=> handles up to W = 600
--SRHD1D.velEpsilon = 1e-7	-- <=> handles up to W = 2,000
--SRHD1D.velEpsilon = 1e-10	-- <=> handles up to W = 100,000
SRHD1D.velEpsilon = 1e-15	-- <=> smaller than 1e-15 gnuplot x11 terminal breaks down past W = 1e+7 ...

SRHD1D.PMinEpsilon = 1e-16
SRHD1D.rhoMinEpsilon = 1e-15



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
	self.gamma = function(self,i) return self.equation.gamma end
	local P = self:calcP(rho, eInt)
	self.gamma = originalGamma

	local h = 1 + eInt + P/rho	
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
	local xmid = (sim.domain.xmin + sim.domain.xmax) / 2
	
	--primitives:
	--[[ Sod
	self.gamma = 7/5
	local rho = x < xmid and 1 or .125
	local vx = 0
	local P = x < xmid and 1 or .1
	--]]
	--[[ Marti & Muller 2008 rppm relativistic shock reflection
	-- not working with my sim...
	self.gamma = 4/3
	local rho = 1
	local vx = 1 - 1e-5
	local P = (self.gamma - 1) * rho * (1e-7 / math.sqrt(1 - vx*vx))
	-- which is approx. P ~ (gamma - 1) rho sqrt(5) 1e-5 ~ sqrt(5)/3 1e-5 ~ 7.5e-6
	--]]
	--[[ relativistic blast wave test problem 1, Marti & Muller 2008, table 5
	-- also the Marti & Muller rppm code's Schneider et al
	-- the paper says P=0.00 for rhs, but looking at the rppm code it should probably be 1.66e-6 
	-- also Odyck problem #1
	self.gamma = 5/3
	local rho = x < xmid and 10 or 1
	local vx = 0
	local P = (self.gamma - 1) * rho * (x < xmid and 2 or 1e-6)
	--]]
	--[[ relativistic blast wave test problem 2, Marti & Muller 2008, table 5
	-- also the relativistic blast wave initial conditions in the provided rppm code
	-- also Odyck problem #2
	self.gamma = 5/3
	local rho = 1
	local vx = 0
	local P = x < xmid and 1000 or .01
	--]]
	-- [[ Marti & Muller 2008 rppm relativistic blast wave interaction
	self.gamma = 7/5
	local lhs = .9 * sim.domain.xmin + .1 * sim.domain.xmax
	local rhs = .1 * sim.domain.xmin + .9 * sim.domain.xmax
	local rho = 1
	local vx = 0
	local P = x < lhs and 1000 or (x > rhs and 100 or .01)
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
	local rhoL, vxL, eIntL = table.unpack(sim.ws[i-1])
	local rhoR, vxR, eIntR = table.unpack(sim.ws[i])
--print('cell',i,'prims L=',rhoL,vxL,eIntL,'R=',rhoR,vxR,eIntR)

--[[ roe averaging.  i need to find a paper that describes how to average all variables, and not just these ...
	local sqrtAbsGL = 1	-- this will change once you start using a proper metric
	local DL, _, _ = table.unpack(sim.qs[i-1])
	local WL = DL / rhoL
	local PL = self:calcP(rhoL, eIntL)
	local hL = 1 + eIntL + PL/rhoL
	local roeWeightL = math.sqrt(sqrtAbsGL * rhoL * hL)
	local u1L = vxL * WL 
	local roeVarsL = {WL, u1L, PL / (rhoL * hL)}
	
	local sqrtAbsGR = 1
	local DR, _, _ = table.unpack(sim.qs[i])
	local WR = DR / rhoR
	local PR = self:calcP(rhoR, eIntR)
	local hR = 1 + eIntR + PR/rhoR
	local roeWeightR = math.sqrt(sqrtAbsGR * rhoR * hR)
	local u1R = vxR * WR 
	local roeVarsR = {WR, u1R, PR / (rhoR * hR)}
	
	local normalize = 1 / (roeWeightL + roeWeightR)
	local roeVars = {}
	for i=1,3 do
		roeVars[i] = (roeVarsL[i] * roeWeightL + roeVarsR[i] * roeWeightR) * normalize
	end
	local W = roeVars[1]	-- = u0
	local u1 = roeVars[2]	-- = vx * W
	local vx = u1 / W
	vx = math.clamp(vx, -1+self.velEpsilon,1-self.velEpsilon)
	-- eInt / h = (p / (rho h)) / (self.gamma - 1) 
	local P_over_rho_h = roeVars[3]

	-- so how do you figure 'h' for the eigenvector calculations?
	local rho_h = roeWeightL * roeWeightR	-- .. divided by sqrt(sqrtAbsGL * sqrtAbsGR)
	local h = (hL * roeWeightL + hR * roeWeightR) * normalize
	local rho = rho_h / h
	--local eInt = (h - 1)/self.gamma	-- ideal gas law
	local P = P_over_rho_h * rho * h
	--local P = self:calcP(rho, eInt)
	--local P_over_rho_h = P / rho_h	-- recalculate
	-- this is giving a high eigenbasis error ... which makes me think some variables are out of sync

	local W2 = W * W
	local vSq = vx * vx
	local oneOverW = 1/W
--]]
-- [[ arithmetic averaging
	local rho = (rhoL + rhoR) / 2
	local eInt = (eIntL + eIntR) / 2
	local vx = (vxL + vxR) / 2
	local vSq = vx * vx
	local oneOverW2 = 1 - vSq
	local oneOverW = math.sqrt(oneOverW2)
	local W = 1/oneOverW
	local W2 = 1/oneOverW2
	local P = self:calcP(rho, eInt)
	local h = 1 + eInt + P/rho
	local P_over_rho_h = P / (rho  * h)
--]]

	local hW = h * W
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
--print('cell',i,'lambda',table.unpack(lambda))

	--[[ general equations
	-- these explode with the two-shockwave test
	--Font 2008 finally says that kappa = dp/deInt
	-- Marti & Muller 2008, Alcubierre 2008, and a few others I'm betting all have the eigenvalues, vectors, all the same ... but forget kappa
	local kappa = (self.gamma - 1) * rho
	-- Font 2008 eqn 112, also Marti & Muller 2008 eqn 74
	-- what do we do when kappaTilde approaches csSq?
	local kappaTilde = kappa / rho
	local Kappa = kappaTilde / (kappaTilde - csSq)
	--]]
	-- [[ for ideal gas
	-- these survive the test for a few more frames, but then the velocity goes next
	local Kappa = h
	--]]

-- [=[ Marti, Muller 2008
	local AMinus= (1 - vxSq) / (1 - vx * lambda[1])
	local APlus = (1 - vxSq) / (1 - vx * lambda[3])
--print('cell',i,'A+',APlus,'A-',AMinus)
--print('cell',i,'h',h,'W',W,'hW',hW)
	local evr = sim.eigenvectors[i]
	evr[1][1] = 1
	evr[2][1] = hW * AMinus * lambda[1]
	evr[3][1] = hW * AMinus - 1
	evr[1][2] = oneOverW	-- Kappa/hW in general
	evr[2][2] = vx
	evr[3][2] = 1 - oneOverW	-- 1 - Kappa/hW in general
	evr[1][3] = 1
	evr[2][3] = hW * APlus * lambda[3]
	evr[3][3] = hW * APlus - 1

	-- [[
	local Delta = h*h*hW * (Kappa - 1) * (1 - vxSq) * (APlus * lambda[3] - AMinus * lambda[1])
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
	--]]
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
--]=]
	--[[ or just do it numerically
	local evl = sim.eigenvectorsInverse[i]
	sim.eigenvectorsInverse[i] = mat33.inv(evr)
	--]]

--[[
for j=1,3 do
	for k=1,3 do
		print('cell',i,'right eigenvectors',j,k,evr[j][k])
	end
end
for j=1,3 do
	for k=1,3 do
		print('cell',i,'left eigenvectors',j,k,evl[j][k])
	end
end
--]]

for j=1,3 do
	assertfinite(lambda[j])
	for k=1,3 do
		if not math.isfinite(evl[j][k])
		or not math.isfinite(evr[j][k])
		then 
			error(tolua({
				evl=evl,
				evr=evr,
				Kappa=Kappa,
				hW=hW,
				h=h,
				W=W,
				csSq=csSq,
				rho=rho,
			}, {indent=true}))
		end
	end
end

	-- define the flux matrix to see how accurate our 
	-- this is coming out different than the matrix reconstructed from the eigensystem
	local W3 = W2 * W
	local W4 = W2 * W2
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
	local _, dw_dU = xpcall(function()
		return mat33.inv(dU_dw)
	end, function(err)
		print(err..'\n'..tolua({
			rho=rho,
			vx=vx,
			P=P,
			h=h,
			W=W,
			W2=W2,
			W3=W3,
		},{indent=true})..'\n'..debug.traceback())
	end)

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

-- this is my attempt based on the recover pressure method described in Alcubierre, Baumgarte & Shapiro, Marti & Muller, Font, and generally everywhere
function SRHD1D:calcPrimsByPressure(sim, i ,prims, qs)
	local rho, vx, eInt = table.unpack(prims)

	-- use the new state variables to newton converge
	local D, Sx, tau = table.unpack(qs)

	-- superluminal velocity will occur if |S| > tau + D + P
	-- i.e. if P < |S| - tau - D
	local absSx = math.abs(Sx)
	local PMin = math.max(absSx - tau - D + absSx * self.velEpsilon, self.PMinEpsilon)

	-- this is in the Marti & Muller 2008 code ...
	-- where do they get this from?
	local PMax = (self.gamma - 1) * tau
	-- why is there an upper bound again?
	--assert(PMax >= PMin, 'pmax='..PMax..' < pmin='..PMin)
	PMax = math.max(PMax, PMin)

	-- start with a guess between PMin and PMax 
	local P = .5 * (PMin + PMax)
assertfinite(P)

	for iter=1,self.solvePrimMaxIter do
		vx = Sx / (tau + D + P)
-- PMin should prevent this from occuring
if math.abs(vx) > 1 then
	error('superluminal vx='..vx..' Sx='..Sx..' tau='..tau..' D='..D..' P='..P)
end
		-- in the case that this happens, the delta P can still be huge, but the min() keeps us in place ... 
		-- if we run up against the boundary then we should exit as well
assertfinite(vx)
		local vSq = vx*vx	-- + vy*vy + vz*vz
assertfinite(vSq)
		local W = 1 / math.sqrt(1 - vSq)
assertfinite(W, 'W nan, vSq='..vSq..' vx='..vx)
		eInt = (tau + D * (1 - W) + P * (1 - W*W)) / (D * W)
assertfinite(eInt, 'eInt nan, tau='..tau..' D='..D..' W='..W..' P='..P)
		rho = D / W
assertfinite(rho)

		-- for f = calcP = (self.gamma - 1) * rho * eInt
		-- Marti & Muller 2008, eqn 61
		local f = self:calcP(rho, eInt) - P
assertfinite(f)
		
		-- Alcubierre, 7.6.52
		--local csSq = self.gamma * P / (rho * P)	
		-- however the code with Marti & Muller uses a different equation for the cs^2 term of df/dp:
		local csSq = (self.gamma - 1) * (tau + D * (1 - W) + P) / (tau + D + P)
assertfinite(csSq)
		local df_dP = vSq * csSq - 1
assertfinite(df_dP)
		
		-- code with Marti & Muller says to keep things above PMin ...
		local newP = P - f / df_dP
		newP = math.max(newP, PMin)
assertfinite(newP)
		
		-- Marti & Muller 2008 code ...
		local PError = math.abs(1 - newP / P)
		P = newP
assertfinite(P)
		if PError < self.solvePrimErrorEpsilon then
			-- one last update ...
			vx = Sx / (tau + D + P)
			W = 1 / math.sqrt(1 - vx*vx)
			rho = D / W
			rho = math.max(rho, self.rhoMinEpsilon)
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
		if err < self.solvePrimErrorEpsilon then
			rho = D / W
			rho = math.max(rho, self.rhoMinEpsilon)
			h = Z / (D * W)		-- h-1 = gamma eInt <=> eInt = (h-1) / gamma  <=> (h-1) = eInt + P/rho
			h = math.max(h, 1 + self.rhoMinEpsilon)
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
assertfinite(D)
assertfinite(Sx)
assertfinite(tau)
	
	--dynamic prims to converge
	local rho, vx, eInt = table.unpack(prims)
assertfinite(rho)
assertfinite(vx)
assertfinite(eInt)
--print('cell',i)
	for iter=1,self.solvePrimMaxIter do
--print('rho='..rho..' vx='..vx..' eInt='..eInt)
		local P = self:calcP(rho, eInt)
		local dP_drho = self:calc_dP_drho(rho, eInt)
		local dP_deInt = self:calc_dP_deInt(rho, eInt)
assertfinite(P)
		local h = 1 + eInt + P/rho
assertfinite(h)
		local W2 = 1/(1-vx*vx)
		local W = math.sqrt(W2)
assertfinite(W, 'W nan with vx='..vx)
		local W3 = W*W2
		local dU_dw = {
			{W, rho * W3 * vx, 0},
			{W2 * vx * (1 + eInt + dP_drho), rho * h * W2 * (2 * W2 - 1), W2 * vx * (dP_deInt + rho)},
			{-dP_drho - W + W2 * (1 + eInt + dP_drho), rho * W3 * vx * (2 * W * h - 1), -dP_deInt + dP_deInt * W2 + rho * W2}
		}
for j=1,3 do
	for k=1,3 do
		assertfinite(dU_dw[j][k], 'nan for dU_dw['..j..']['..k..']')
	end
end
		local det_dU_dw = mat33.det(dU_dw)
assertfinite(det_dU_dw, 'nan for derivative of:\n'..require'matrix'(dU_dw))
		local dw_dU = mat33.inv(dU_dw, det_dU_dw)
for j=1,3 do
	for k=1,3 do
		assertfinite(dw_dU[j][k])
	end
end
		local dU = {
			-D + rho * W,
			-Sx + rho * h * W2 * vx,
			-tau + rho * h * W2 - P - rho * W,
		}
assertfinite(dU[1], 'dU[1] nan with D='..D..' rho='..rho..' W='..W)
for j=2,3 do
	assertfinite(dU[j], 'nan for dU['..j..']')
end
		drho = dw_dU[1][1] * dU[1] + dw_dU[1][2] * dU[2] + dw_dU[1][3] * dU[3]
assertfinite(drho)
		dvx = dw_dU[2][1] * dU[1] + dw_dU[2][2] * dU[2] + dw_dU[2][3] * dU[3]
assertfinite(dvx)
		deInt = dw_dU[3][1] * dU[1] + dw_dU[3][2] * dU[2] + dw_dU[3][3] * dU[3]
assertfinite(deInt)
		local alpha = 1
		local new_rho = rho - alpha * drho
		local new_vx = vx - alpha * dvx
		local new_eInt = eInt - alpha * deInt
		-- [[ enforce boundaries
		new_rho = math.clamp(new_rho,1e-10,1e+20)
		new_vx = math.clamp(new_vx,-1+self.velEpsilon,1-self.velEpsilon)
		new_eInt = math.clamp(new_eInt,0,1e+20)
		--]]
		local err = math.abs(new_rho-rho) + math.abs(new_vx-vx) + math.abs(new_eInt-eInt)
		rho = new_rho
		vx = new_vx
		eInt = new_eInt
		local cons = table{self:consFromPrims(rho,vx,eInt)}
		local consErr = table{math.abs(cons[1]-qs[1]), math.abs(cons[2]-qs[2]), math.abs(cons[3]-qs[3])}
--print('err='..err)
		if err < self.solvePrimErrorEpsilon then
			return {rho, vx, eInt}, nil, iter
		end
	end
	return nil, "didn't converge", self.solvePrimMaxIter
end

return SRHD1D
