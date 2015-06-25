local class = require 'ext.class'
local fluxLimiters = require 'limiter' 

-- determine timestep based on eigenvalue interfaces
local function calcDtWithEigenvalues(self, sim)
	if sim.fixed_dt then
		return sim.fixed_dt
	else
		local result = huge
		for i=1,sim.gridsize do
			local eigenvaluesL = sim.eigenvalues[i]
			local eigenvaluesR = sim.eigenvalues[i+1]
			local maxLambda = max(0, unpack(eigenvaluesL))
			local minLambda = min(0, unpack(eigenvaluesR))
			local dx = sim.ixs[i+1] - sim.ixs[i]
			local dum = dx / (abs(maxLambda - minLambda) + 1e-9)
			result = min(result, dum)
		end
		return result * sim.cfl
	end
end

local Roe = class()

-- calculates timestep and eigenbasis
function Roe:calcDT(sim, getLeft, getRight)
	-- Roe solver:
	-- 1) calculate eigenbasis at interface with Roe-weighted average of states
	for i=2,sim.gridsize do
		local qL = getLeft and getLeft(i) or sim.qs[i-1]
		local qR = getRight and getRight(i) or sim.qs[i]
		
		sim:calcInterfaceEigenBasis(i,qL,qR)

		-- collect error for reporting
		local eigenbasisError = 0
		local fluxMatrixError = 0
		-- for the i'th cell:
		-- Q = eigenvectors matrix
		-- L = eigenvalues matrix
		-- A_jk = Q_jl L_l invQ_lk
		-- eigenbasisError = Q_jl invQ_lk - delta_jk
		-- fluxMatrixError = Q_jl L_l invQ_lk - A_jk
		for j=1,sim.numStates do
			-- local basis_j = 0's everywhere except a 1 at the j'th entry
			-- local eigencoords_j = {k,eigenfield[i][k](basis_j)}			<- dot input vector with eigenvector inverse row k
			-- local eigenscaled_j = eigencoords_j * lambda_j
			-- local newbasis_j = {k,eigenfieldInverse[i][k](eigencoords_j)}	<- dot input vector with eigenvector row k
			-- local newtransformed_j = {k,eigenfieldInverse[i][k](eigenscaled_j)}
			-- sum up abs error between basis_j and newbasis_j 
			-- sum up abs error between A_jk basis_k and newtransformed_j
		
			-- basis_k = delta_jk
			local basis = {}
			for k=1,sim.numStates do
				basis[k] = k == j and 1 or 0
			end

			-- eigenCoords_k = invQ_kl basis_l
			local eigenCoords = sim:eigenfields(i, basis)
			for k=1,sim.numStates do
				assert(type(eigenCoords[k])=='number', "failed for coord "..k.." got type "..type(eigenCoords[k]))
			end

			-- eigenScaled_k = lambda_k * eigenCoords_k
			local eigenScaled = {}
			for k=1,sim.numStates do
				eigenScaled[k] = sim.eigenvalues[i][k] * eigenCoords[k]
			end
			
			-- newbasis_k = Q_kl eigenCoords_l
			local newbasis = sim:eigenfieldsInverse(i, eigenCoords)
			
			-- newtransformed_k = Q_kl eigenScaled_l = Q_kl lambda_l eigenCoords_k
			local newtransformed = sim:eigenfieldsInverse(i, eigenScaled)

			-- transformed_k = A_kl basis_l
			local transformed = sim:fluxTransform(i, basis)

			for k=1,sim.numStates do
				eigenbasisError = eigenbasisError + math.abs(basis[k] - newbasis[k])
				fluxMatrixError = fluxMatrixError + math.abs(transformed[k] - newtransformed[k])
			end
		end
		sim.eigenbasisErrors[i] = eigenbasisError
		sim.fluxMatrixErrors[i] = fluxMatrixError
	end

	return calcDtWithEigenvalues(self, sim)
end

function Roe:calcFlux(sim, dt, getLeft, getRight, getLeft2, getRight2)
	
	-- 2) calculate interface state difference in eigenbasis coordinates
	for i=2,sim.gridsize do
		local qL = getLeft and getLeft(i) or sim.qs[i-1]
		local qR = getRight and getRight(i) or sim.qs[i]
		
		local dq = {}
		for j=1,sim.numStates do
			dq[j] = qR[j] - qL[j]
		end
		sim.deltaQTildes[i] = sim:eigenfields(i, dq)
	end

	local useFluxMatrix = false

	-- 3) slope limit on interface difference
	-- 4) transform back
	for i=2,sim.gridsize do
		local qL = getLeft and getLeft(i) or sim.qs[i-1]
		local qR = getRight and getRight(i) or sim.qs[i]

		local qAvg = {}
		for j=1,sim.numStates do
			qAvg[j] = .5 * (qR[j] + qL[j])
		end
			
		local rTildes = {}
		for j=1,sim.numStates do
			if sim.deltaQTildes[i][j] == 0 then
				rTildes[j] = 0
			else
				if sim.eigenvalues[i][j] >= 0 then
					rTildes[j] = sim.deltaQTildes[i-1][j] / sim.deltaQTildes[i][j]
				else
					rTildes[j] = sim.deltaQTildes[i+1][j] / sim.deltaQTildes[i][j]
				end
			end
		end

		local fluxTilde = {}
		for j=1,sim.numStates do
			local phi = sim.fluxLimiter(rTildes[j])
			local theta = sim.eigenvalues[i][j] >= 0 and 1 or -1
			local dx = sim.xs[i] - sim.xs[i-1]
			local epsilon = sim.eigenvalues[i][j] * dt / dx
			local deltaFluxTilde = sim.eigenvalues[i][j] * sim.deltaQTildes[i][j]
			fluxTilde[j] = -.5 * deltaFluxTilde * (theta + phi * (epsilon - theta))
		end
		
		if not useFluxMatrix then
			local qAvgTildes = sim:eigenfields(i, qAvg)
			for j=1,sim.numStates do
				fluxTilde[j] = fluxTilde[j] + sim.eigenvalues[i][j] * qAvgTildes[j]
			end
		end
		
		sim.fluxes[i] = sim:eigenfieldsInverse(i, fluxTilde)
		
		-- using the flux matrix itself allows for complete reconstruction even in the presence of zero-sim.eigenvalues
		if useFluxMatrix then
			local fluxQs = sim:fluxTransform(i, qAvg)
			for j=1,sim.numStates do
				sim.fluxes[i][j] = sim.fluxes[i][j] + fluxQs[j]
			end
		end
	end

	local dq_dts = sim:newState()
	for i=1,sim.gridsize do
		local dx = sim.ixs[i+1] - sim.ixs[i]
		for j=1,sim.numStates do
			dq_dts[i][j] = dq_dts[i][j] - (sim.fluxes[i+1][j] - sim.fluxes[i][j]) / dx
		end
	end
	return dq_dts
end

--[[

  d [ rho ]    d    [ rho ]     d    [ 0 ]
  - [rho u] +  - (u [rho u]) +  - (P [ 1 ]) = 0
 dt [rho e]   dx    [rho e]    dx    [ u ]

the 1st and 2nd terms are integrated via the flux integration
the 1st and 3rd terms are integrated via the pressure integration
	that is split into first the momentum and then the work diffusion 

--]]
local EulerBurgers = class()

function EulerBurgers:calcDT(sim, getLeft, getRight)
	local gamma = sim.gamma
	
	-- determine timestep based on cell velocity 
	local dt
	if sim.fixed_dt then
		dt = sim.fixed_dt
	else
		local result = huge
		for i=1,sim.gridsize do
			local rho = sim.qs[i][1]
			local u = sim.qs[i][2] / rho
			local eTotal = sim.qs[i][3] / rho
			local eInt = eTotal - .5 * u * u
			local Cs = sqrt(gamma * (gamma - 1) * eInt)

			local dx = sim.ixs[i+1] - sim.ixs[i]
			local dum = dx / (Cs + abs(u))
			result = min(result, dum)
		end
		dt = result * sim.cfl
	end

	return dt
end


function EulerBurgers:calcFlux(sim, dt, getLeft, getRight, getLeft2, getRight2)
	local dq_dts = sim:newState()
	
	local gamma = sim.gamma

	for i=3,sim.gridsize-1 do
		local qL2 = getLeft2 and getLeft2(i) or sim.qs[i-2]
		local qL = getLeft and getLeft(i) or sim.qs[i-1]
		local qR = getRight and getRight(i) or sim.qs[i]
		local qR2 = getRight2 and getRight2(i) or sim.qs[i+1]
		local uL = qL[2] / qL[1]
		local uR = qR[2] / qR[1]
		local iu = .5 * (uL + uR)
		local dx = sim.xs[i] - sim.xs[i-1]
		for j=1,sim.numStates do
			local dq = qR[j] - qL[j]
			local r = dq == 0 and 0 or 
				(iu >= 0 
					and ((qL[j] - qL2[j]) / dq) 
					or ((qR2[j] - qR[j]) / dq))
			local phi = sim.fluxLimiter(r)
			sim.fluxes[i][j] = 
				iu * (iu >= 0 and sim.qs[i-1][j] or sim.qs[i][j])
				-- why does a + make things worse, and a - make things *much* smoother?  (original: http://www.mpia.de/homes/dullemon/lectures/fluiddynamics/Chapter_4.pdf says + )
				- dq * phi * .5 * abs(iu) * (1 - abs(iu * dt/dx))
		end
	end

	for i=1,sim.gridsize do
		local dx = sim.ixs[i+1] - sim.ixs[i]
		for j=1,sim.numStates do
			dq_dts[i][j] = dq_dts[i][j] - (sim.fluxes[i+1][j] - sim.fluxes[i][j]) / dx
		end
	end

	return dq_dts
end
	
function EulerBurgers:postIterate(sim, dt)
	sim.Ps = sim.Ps or {}
	
	local gamma = sim.gamma

	for i=1,sim.gridsize do
		local rho = sim.qs[i][1]
		local u = sim.qs[i][2] / rho
		local eTotal = sim.qs[i][3] / rho
		local eInt = eTotal - .5 * u * u
		sim.Ps[i] = (gamma - 1) * rho * eInt
	end

	--[[ artificial viscosity
	for i=2,sim.gridsize-1 do
		local rho = sim.qs[i][1]
		local uL = sim.qs[i-1][2] / sim.qs[i-1][1]
		local uR = sim.qs[i+1][2] / sim.qs[i+1][1]
		local zeta = 2
		local du = zeta * .5 * (uR - uL)
		sim.Ps[i] = sim.Ps[i] + du * du * rho
	end
	--]]

	-- diffuse momentum
	for i=2,sim.gridsize-1 do
		local dP = sim.Ps[i+1] - sim.Ps[i-1]
		local dx = sim.xs[i+1] - sim.xs[i-1]
		sim.qs[i][2] = sim.qs[i][2] - dP * dt / dx
	end

	-- diffuse work
	for i=2,sim.gridsize-1 do
		local WR = sim.Ps[i+1] * sim.qs[i+1][2] / sim.qs[i+1][1]
		local WL = sim.Ps[i-1] * sim.qs[i-1][2] / sim.qs[i-1][1]
		local dW = WR - WL
		local dx = sim.xs[i+1] - sim.xs[i-1]
		sim.qs[i][3] = sim.qs[i][3] - dW * dt / dx
	end
end

local HLL = class()

function HLL:calcDT(sim, getLeft, getRight)
	-- matches Roe, except without eigenvectors
	for i=2,sim.gridsize do
		local qL = getLeft and getLeft(i) or sim.qs[i-1]
		local qR = getRight and getRight(i) or sim.qs[i]
		sim:calcInterfaceEigenvalues(qL, qR, sim.eigenvalues[i])
	end

	return calcDtWithEigenvalues(self, sim)
end
	
function HLL:calcFlux(sim, dt, getLeft, getRight, getLeft2, getRight2)
	local gamma = sim.gamma
	
	local iqs = sim:newState()
	
	for i=2,sim.gridsize do
		-- TODO use qL and qR to allow compatability with MUSCL
		local qL = getLeft and getLeft(i) or sim.qs[i-1]
		local qR = getRight and getRight(i) or sim.qs[i]
		
		local sL = sim.eigenvalues[i][1]
		local sR = sim.eigenvalues[i][sim.numStates]

		local fluxL = sim:calcFluxForState(qL)
		local fluxR = sim:calcFluxForState(qR)

		local flux = sim.fluxes[i]
		for i=1,sim.numStates do
			if 0 <= sL then
				flux[i] = fluxL[i]
			elseif sL <= 0 and 0 <= sR then
				flux[i] = (sR * fluxL[i] - sL * fluxR[i] + sL * sR * (qR[i] - qL[i])) / (sR - sL)
			elseif sR <= 0 then
				flux[i] = fluxR[i]
			end
		end
	end
	
	local dq_dts = sim:newState()
	for i=1,sim.gridsize do
		local dx = sim.ixs[i+1] - sim.ixs[i]
		for j=1,sim.numStates do
			dq_dts[i][j] = dq_dts[i][j] - (sim.fluxes[i+1][j] - sim.fluxes[i][j]) / dx
		end
	end

	return dq_dts
end

local EulerMUSCL = class()

function EulerMUSCL:init(args)
	-- baseScheme can be Roe or HLL
	-- TODO should MUSCL-Roe be using eigenbasis-transformed states somewhere in there? 
	self.baseScheme = args.baseScheme or require 'scheme'.Roe()

	-- limiter of ratio
	-- popular limiters: Fromm, Beam-Warming, Lax-Wendroff, minmod
	-- notice that winded-ness needs to vary per-scheme
	self.slopeLimiter = args.slopeLimiter or fluxLimiters.Fromm
end

-- first calculate new state interface left and right
-- then call the old calcDT which also calculates eigen basis stuff
function EulerMUSCL:calcDT(sim, getLeft, getRight, getLeft2, getRight2)
	
	local sigma = sim:newState()
	for i=3,sim.gridsize-1 do
		local qL2 = getLeft2 and getLeft2(i) or sim.qs[i-2]
		local qL = getLeft and getLeft(i) or sim.qs[i-1]
		local qR = getRight and getRight(i) or sim.qs[i]
		local qR2 = getRight2 and getRight2(i) or sim.qs[i+1]
		for j=1,sim.numStates do
			local dq = qR[j] - qL[j]

			-- since we're doing flux/subgrid per-state (un-transformed)
			--  then, for determining the next-cell of advection,
			--   I'm going to use the velocity
			--   ...which means this is a Euler-only MUSCL solver
			local iu = .5*(qL[2]/qL[1] + qR[2]/qR[1])
			
			-- ratio of slope to next cell slope
			local r = dq == 0 and 0 or 
				(iu >= 0 
					and ((qL[j] - qL2[j]) / dq) 
					or ((qR2[j] - qR[j]) / dq))
			
			local phi = self.slopeLimiter(r)
			
			-- slope limiter:
			-- 	sigma = minmod(dq[i],dq[i']) for i current cell and i' next advection cell
			-- flux limiter:
			--  phi = minmod(1,r) = minmod(dq[i]/dq[i], dq[i']/dq[i])
			-- ... so scaling phi by dq[i] ... or simply not dividing r by dq[i] ... should get us the slope
			-- ... but if the slope limiter is constrained to [0,2] then we want it divided at first, then multiply after 
			sigma[i][j] = phi * dq
		end
	end

	-- subgrid values
	local iqLs = sim:newState()
	local iqRs = sim:newState()
	for j=1,sim.numStates do
		iqLs[1][j] = sim.qs[1][j]
		iqRs[1][j] = sim.qs[1][j]
	end
	for i=2,sim.gridsize do
		local qL = getLeft and getLeft(i) or sim.qs[i-1]
		local qR = getRight and getRight(i) or sim.qs[i]
		local dx = sim.ixs[i] - sim.ixs[i-1]
		for j=1,sim.numStates do
			iqLs[i][j] = qL[j] + .5 * dx * sigma[i-1][j]
			iqRs[i][j] = qR[j] - .5 * dx * sigma[i][j]
		end
	end

	-- flux based on subgrid values
	local ifLs = sim:newState()
	local ifRs = sim:newState()
	for i=1,sim.gridsize do
		ifLs[i] = sim:calcFluxForState(iqLs[i])
		ifRs[i] = sim:calcFluxForState(iqRs[i])
	end

	-- how should we get this dt?
	-- based on Roe eigenvector without MUSCL?  or based on Burgers?  or fixed + implicit (optimistically)
	-- should this be the dt that we consistently use even after using MUSCL to adjust the state?
	local dt = self.baseScheme.calcDT(self, sim)

	-- half step in time
	self.iqhLs = sim:newState()
	self.iqhRs = sim:newState()
	for j=1,sim.numStates do
		self.iqhLs[1][j] = iqLs[1][j]
		self.iqhRs[1][j] = iqRs[1][j]
		self.iqhLs[sim.gridsize][j] = iqLs[sim.gridsize][j]
		self.iqhRs[sim.gridsize][j] = iqRs[sim.gridsize][j]
	end
	for i=2,sim.gridsize-1 do
		local dx = sim.ixs[i] - sim.ixs[i-1]
		for j=1,sim.numStates do
			self.iqhLs[i][j] = iqLs[i][j] + .5 * dt/dx * (ifLs[i][j] - ifRs[i-1][j])
			self.iqhRs[i][j] = iqRs[i][j] + .5 * dt/dx * (ifLs[i+1][j] - ifRs[i][j])
		end
	end
	
	-- once we have *this* collection of subgrid L & R states,
	--  we use them for whatever method you want ...
	-- ... be it HLL or Roe, etc 

	self.getLeft = function(i) return self.iqhLs[i] end
	self.getRight = function(i) return self.iqhRs[i] end
	self.getLeft2 = function(i) return self.iqhLs[i-1] end
	self.getRight2 = function(i) return self.iqhRs[i+1] end

	--[[
	TODO now the paper has officially gave me a circular dependency:
	to do the MUSCL subgrid half-step in time it needs to know dt
	however, to calculate dt, it needs the eigenvalues
		eigenvalues need Roe matrix at interface
		interface needs the left and right half-step values
		... which need the dt
	--]]
	local dt = self.baseScheme.calcDT(self, sim, self.getLeft, self.getRight)

	return dt
end

function EulerMUSCL:calcFlux(sim, dt)
	return self.baseScheme:calcFlux(sim, dt, self.getLeft, self.getRight, self.getLeft2, self.getRight2)
end

function EulerMUSCL:postIterate(...)
	if self.baseScheme.postIterate then
		return self.baseScheme:postIterate(...)
	end
end

function EulerMUSCL:addSourceToDeriv(...)
	if self.baseScheme.addSourceToDeriv then
		self.baseScheme:addSourceToDeriv(...)
	end
end

return {
	Roe = Roe,
	HLL = HLL,
	-- only works for Eulers
	EulerBurgers = EulerBurgers,
	EulerMUSCL = EulerMUSCL,
}

