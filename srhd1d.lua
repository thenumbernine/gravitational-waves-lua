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
	local eInt = eInt * rho
	local D = q:_(1)	-- rho * W
	local Sx = q:_(2)	-- rho h W^2 vx
	local tau = q:_(3)	-- rho h W^2 - P - D
	local W = D / rho	-- Lorentz factor
	SRHD1D:buildGraphInfos{
		{rho = rho},
		{vx = vx},
		{P = P},
		{eInt = eInt},
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

local function consFromPrim(rho, vx, eInt)
	local WSq = 1/(1-vx^2)
	local W = math.sqrt(WSq)
	local P = pressureFunc(rho, eInt)
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
	-- [[ Sod
	-- density 
	local rho = x < 0 and 1 or .125
	-- local velocity
	local vx = 0
	-- local pressure
	local P = x < 0 and 1 or .1
	--]]
	-- specific internal energy ... unspecified by the paper, but I'll use ideal gas heat capacity ratio
	local eInt = energyInternalFromPressure(rho, P)
	--aux variables:
	local vSq = vx^2 -- + vy^2 + vz^2
	local W = 1/math.sqrt(1 - vSq)
	-- specific enthalpy
	local h = 1 + eInt + P/rho
	-- store primitive variables for later use
	sim.primitives[i] = {rho=rho, vx=vx, h=h, eInt=eInt, P=P}
	--conservative variables:
	return {consFromPrim(rho, vx, eInt)}
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

--[[ roe averaging.  i need to find a paper that describes how to average all variables, and not just these ...
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
-- [[ arithmetic averaging
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

--[[ this is my attempt based on the recover pressure method described in Alcubierre, Baumgarte & Shapiro, Marti & Muller, Font, and generally everywhere
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
		if math.abs(PError) < 1e-8 then
			-- one last update ...
			vx = Sx / (tau + D + P)
			W = 1 / math.sqrt(1 - vSq)
			rho = D / W
			eInt = P / (rho * (gamma - 1))
			h = 1 + eInt + P/rho
			-- assign prims
			prims.P = P
			prims.rho = rho
			prims.vx = vx 
			--if prims.vx < velEpsilon^2 then prims.vx = 0 end
			-- Marti & Muller 2008 recompute eInt here
			prims.eInt = eInt 
			prims.h = h
			
			-- check reconstruction error:
			local cons = {consFromPrim(rho, vx, eInt)}
			for j=1,3 do
				local err = math.abs(cons[j]-qs[j])
				if err > 1e-7 then
					print('got bad prim reconstruction of cell '..i..' with '..({'D','Sx','tau'})[j]..' error at '..err)
					print('original qs:',table.unpack(qs))
					print('new reconstruction:',table.unpack(cons))
					error'here'
				end
			end
			break
		end
		assert(iter ~= maxiter, "didn't converge")
	end
end
--]]

--[[ here's the method described in Anton & Zanotti 2006, used by Mara:
-- Mara converges Z=rho h W^2 and W, while the paper says to converge rho and W
-- either way, if you converge W then you're calculating things based on a W apart from your conservative W
-- unless you replaced the conservative W with the converged W
function SRHD1D:calcPrims(sim,i,prims,qs)
	local D, Sx, tau = table.unpack(qs)
	local vSq = prims.vx*prims.vx
	local W = 1/math.sqrt(1-vSq)-- W is the other
	local P = prims.P
	local rho = prims.rho
	local eInt = prims.eInt
	local h = prims.h
	local Z = rho * h * W*W -- Z = rho h W^2 is one var we converge
	local maxiter = 100
	for iter=1,maxiter do
		rho = D / W
		h = Z / (D * W)		-- h-1 = gamma eInt <=> eInt = (h-1) / gamma  <=> (h-1) = eInt + P/rho
		eInt = (h - 1) / gamma
		P = (gamma - 1) * rho * eInt
	
		-- f = [-S^2 + Z^2 (W^2 - 1) / W^2, -tau + Z^2 - P - D]
		local f = {
			-Sx*Sx  + Z*Z * (W*W - 1) / (W*W),
			-tau +  Z - P - D,
		}
		local df_dx = {
			{2*Z * (W*W - 1) * Z*Z*Z / (W*W*Z*Z*Z), 2*Z*Z / (W*W*W)},
			{1 - (gamma - 1) / (gamma * W*W), (2*Z - D*W)/(W*W*W) * (gamma-1)/gamma},
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
		if err < 1e-12 then
			rho = D / W
			h = Z / (D * W)		-- h-1 = gamma eInt <=> eInt = (h-1) / gamma  <=> (h-1) = eInt + P/rho
			eInt = (h - 1) / gamma
			P = (gamma - 1) * rho * eInt		
			
			prims.P = P
			prims.rho = rho
			prims.vx = math.clamp(Sx / (tau + D + P), 0, 1)
			prims.eInt = P / (rho * (gamma - 1))
			prims.h = 1 + eInt + P/rho
			break
		end
		if iter == maxiter then
			print('hit maxiter')
		end
	end
end
--]]

--[[
how about 3-var newton descent?   D,Sx,tau -> rho,vx,eInt
newton update: w = w - dw/dU * U
U = the zeros of the cons eqns
-D' + rho * W
-Sx' + rho * h * W^2 * vx
-tau' + rho * h * W^2 - P - D
...for D', Sx', tau' fixed
...and all else calculated based on rho,vx,eInt
dU/dw = 
--]]
function SRHD1D:calcPrims(sim,i,prims,qs)
	--fixed cons to converge to
	local D, Sx, tau = table.unpack(qs)
	--dynamic prims to converge
	local rho, vx, eInt = prims.rho, prims.vx, prims.eInt
	local maxiter = 100
	--print('cell',i)
	for iter=1,maxiter do
		--print('rho='..rho..' vx='..vx..' eInt='..eInt)
		local P = pressureFunc(rho,eInt)
		local h = 1 + eInt + P/rho
		local W = 1/math.sqrt(1-vx*vx)
		local dU_dw = {
			{W, rho*W*W*W*vx, 0},
			{h*W*W*vx, rho*h*W*W*(2*W*W-1), gamma*rho*W*W*vx},
			{1*((h-1)/gamma - W + W*W*h + 1 - h), rho*W*W*W*vx*(2*W*h-1), rho*(gamma*W^2 - (gamma-1))}
		}
		local dw_dU = mat33.inv(dU_dw)
		local dU = {
			-D + rho*W,
			-Sx + rho*h*W*W*vx,
			-tau + rho*h*W*W - P - rho*W,
		}
		drho = dw_dU[1][1] * dU[1] + dw_dU[1][2] * dU[2] + dw_dU[1][3] * dU[3]
		dvx = dw_dU[2][1] * dU[1] + dw_dU[2][2] * dU[2] + dw_dU[2][3] * dU[3]
		deInt = dw_dU[3][1] * dU[1] + dw_dU[3][2] * dU[2] + dw_dU[3][3] * dU[3]
		local alpha = 1
		local new_rho = rho - alpha * drho
		local new_vx = vx - alpha * dvx
		local new_eInt = eInt - alpha * deInt
		-- [[ enforce boundaries
		new_rho = math.max(new_rho,1e-7)
		new_vx = math.clamp(new_vx,-1+1e-7,1-1e-7)
		new_eInt = math.max(new_eInt,0)
		--]]
		local err = math.abs(new_rho-rho) + math.abs(new_vx-vx) + math.abs(new_eInt-eInt)
		rho = new_rho
		vx = new_vx
		eInt = new_eInt
		--print('err='..err)
		if err < 1e-7 then
			prims.rho = rho
			prims.vx = vx
			prims.eInt = eInt
			prims.P = pressureFunc(rho,eInt)
			prims.h = 1 + eInt + P/rho
			--[[ make sure the reconstructed conservative variables match
			local cons = {consFromPrim(rho,vx,eInt)}
			for j=1,3 do
				local err = math.abs(cons[j]-qs[j])
				if err > 1e-5 then
					print('got bad prim reconstruction of cell '..i..' with '..({'D','Sx','tau'})[j]..' error at '..err)
					print('original qs:',table.unpack(qs))
					print('new reconstruction:',table.unpack(cons))
					error('reconstruction error='..err)
				end
			end		
			--]]
			-- [[ just replace them
			sim.qs[i] = {consFromPrim(rho,vx,eInt)}
			--]]
			break
		end
		if iter == maxiter then
			error("didn't converge")
		end
	end
end

return SRHD1D
