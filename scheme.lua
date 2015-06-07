-- determine timestep based on eigenvalue interfaces
local function calcDtWithEigenvalues(self)
	if self.fixed_dt then
		return self.fixed_dt
	else
		local result = huge
		for i=1,self.gridsize do
			local eigenvaluesL = self.eigenvalues[i]
			local eigenvaluesR = self.eigenvalues[i+1]
			local maxLambda = max(0, unpack(eigenvaluesL))
			local minLambda = min(0, unpack(eigenvaluesR))
			local dx = self.ixs[i+1] - self.ixs[i]
			local dum = dx / (abs(maxLambda - minLambda) + 1e-9)
			result = min(result, dum)
		end
		return result * self.cfl
	end
end

local Roe = {
	-- calculates timestep and eigenbasis
	calcDT = function(self, getLeft, getRight)
		-- Roe solver:
		-- 1) calculate eigenbasis at interface with Roe-weighted average of states
		for i=2,self.gridsize do
			local qL = getLeft and getLeft(i) or self.qs[i-1]
			local qR = getRight and getRight(i) or self.qs[i]
			
			self:calcInterfaceEigenBasis(i,qL,qR)

			-- collect error for reporting
			local eigenbasisError = 0
			local fluxMatrixError = 0
			-- for the i'th cell:
			-- Q = eigenvectors matrix
			-- L = eigenvalues matrix
			-- A_jk = Q_jl L_l invQ_lk
			-- eigenbasisError = Q_jl invQ_lk - delta_jk
			-- fluxMatrixError = Q_jl L_l invQ_lk - A_jk
			for j=1,self.numStates do
				-- local basis_j = 0's everywhere except a 1 at the j'th entry
				-- local eigencoords_j = {k,eigenfield[i][k](basis_j)}			<- dot input vector with eigenvector inverse row k
				-- local eigenscaled_j = eigencoords_j * lambda_j
				-- local newbasis_j = {k,eigenfieldInverse[i][k](eigencoords_j)}	<- dot input vector with eigenvector row k
				-- local newtransformed_j = {k,eigenfieldInverse[i][k](eigenscaled_j)}
				-- sum up abs error between basis_j and newbasis_j 
				-- sum up abs error between A_jk basis_k and newtransformed_j
			
				-- basis_k = delta_jk
				local basis = {}
				for k=1,self.numStates do
					basis[k] = k == j and 1 or 0
				end

				-- eigenCoords_k = invQ_kl basis_l
				local eigenCoords = self:eigenfields(i, basis)
				for k=1,self.numStates do
					assert(type(eigenCoords[k])=='number', "failed for coord "..k.." got type "..type(eigenCoords[k]))
				end

				-- eigenScaled_k = lambda_k * eigenCoords_k
				local eigenScaled = {}
				for k=1,self.numStates do
					eigenScaled[k] = self.eigenvalues[i][k] * eigenCoords[k]
				end
				
				-- newbasis_k = Q_kl eigenCoords_l
				local newbasis = self:eigenfieldsInverse(i, eigenCoords)
				
				-- newtransformed_k = Q_kl eigenScaled_l = Q_kl lambda_l eigenCoords_k
				local newtransformed = self:eigenfieldsInverse(i, eigenScaled)

				-- transformed_k = A_kl basis_l
				local transformed = self:fluxTransform(i, basis)

				for k=1,self.numStates do
					eigenbasisError = eigenbasisError + math.abs(basis[k] - newbasis[k])
					fluxMatrixError = fluxMatrixError + math.abs(transformed[k] - newtransformed[k])
				end
			end
			self.eigenbasisErrors[i] = eigenbasisError
			self.fluxMatrixErrors[i] = fluxMatrixError
		end

		return calcDtWithEigenvalues(self)
	end,
	calcFlux = function(self, dt)
		
		-- 2) calculate interface state difference in eigenbasis coordinates
		for i=2,self.gridsize do
			local dq = {}
			for j=1,self.numStates do
				dq[j] = self.qs[i][j] - self.qs[i-1][j]
			end
			self.deltaQTildes[i] = self:eigenfields(i, dq)
		end
	
		local useFluxMatrix = false

		-- 3) slope limit on interface difference
		-- 4) transform back
		for i=2,self.gridsize do
			
			local qAvg = {}
			for j=1,self.numStates do
				qAvg[j] = .5 * (self.qs[i][j] + self.qs[i-1][j])
			end
				
			local rTildes = {}
			for j=1,self.numStates do
				if self.deltaQTildes[i][j] == 0 then
					rTildes[j] = 0
				else
					if self.eigenvalues[i][j] >= 0 then
						rTildes[j] = self.deltaQTildes[i-1][j] / self.deltaQTildes[i][j]
					else
						rTildes[j] = self.deltaQTildes[i+1][j] / self.deltaQTildes[i][j]
					end
				end
			end

			local fluxTilde = {}
			for j=1,self.numStates do
				local phi = self.fluxLimiter(rTildes[j])
				local theta = self.eigenvalues[i][j] >= 0 and 1 or -1
				local dx = self.xs[i] - self.xs[i-1]
				local epsilon = self.eigenvalues[i][j] * dt / dx
				local deltaFluxTilde = self.eigenvalues[i][j] * self.deltaQTildes[i][j]
				fluxTilde[j] = -.5 * deltaFluxTilde * (theta + phi * (epsilon - theta))
			end
			
			if not useFluxMatrix then
				local qAvgTildes = self:eigenfields(i, qAvg)
				for j=1,self.numStates do
					fluxTilde[j] = fluxTilde[j] + self.eigenvalues[i][j] * qAvgTildes[j]
				end
			end
			
			self.fluxes[i] = self:eigenfieldsInverse(i, fluxTilde)
			
			-- using the flux matrix itself allows for complete reconstruction even in the presence of zero-self.eigenvalues
			if useFluxMatrix then
				local fluxQs = self:fluxTransform(i, qAvg)
				for j=1,self.numStates do
					self.fluxes[i][j] = self.fluxes[i][j] + fluxQs[j]
				end
			end
		end

		local dq_dts = self:newState()
		for i=1,self.gridsize do
			local dx = self.ixs[i+1] - self.ixs[i]
			for j=1,self.numStates do
				dq_dts[i][j] = dq_dts[i][j] - (self.fluxes[i+1][j] - self.fluxes[i][j]) / dx
			end
		end
		return dq_dts
	end,
}

local Burgers = {
	calcDT = function(self)
		local gamma = self.gamma
		
		-- determine timestep based on cell velocity 
		local dt
		if self.fixed_dt then
			dt = self.fixed_dt
		else
			local result = huge
			for i=1,self.gridsize do
				local rho = self.qs[i][1]
				local u = self.qs[i][2] / rho
				local eTotal = self.qs[i][3] / rho
				local eInt = eTotal - .5 * u * u
				local Cs = sqrt(gamma * (gamma - 1) * eInt)

				local dx = self.ixs[i+1] - self.ixs[i]
				local dum = dx / (Cs + abs(u))
				result = min(result, dum)
			end
			dt = result * self.cfl
		end

		return dt
	end,
}

local EulerBurgers = {
	calcDT = Burgers.calcDT,
	
	calcFlux = function(self, dt)
		local dq_dts = self:newState()
		
		local gamma = self.gamma

		for i=3,self.gridsize-1 do
			local qL2 = self.qs[i-2]
			local qL = self.qs[i-1]
			local qR = self.qs[i]
			local qR2 = self.qs[i+1]
			local uL = qL[2] / qL[1]
			local uR = qR[2] / qR[1]
			local iu = .5 * (uL + uR)
			local dx = self.xs[i] - self.xs[i-1]
			for j=1,self.numStates do
				local dq = qR[j] - qL[j]
				local r = dq == 0 and 0 or 
					(iu >= 0 
						and ((qL[j] - qL2[j]) / dq) 
						or ((qR2[j] - qR[j]) / dq))
				local phi = self.fluxLimiter(r)
				self.fluxes[i][j] = 
					iu * (iu >= 0 and self.qs[i-1][j] or self.qs[i][j])
					-- why does a + make things worse, and a - make things *much* smoother?  (original: http://www.mpia.de/homes/dullemon/lectures/fluiddynamics/Chapter_4.pdf says + )
					- dq * phi * .5 * abs(iu) * (1 - abs(iu * dt/dx))
			end
		end

		for i=1,self.gridsize do
			local dx = self.ixs[i+1] - self.ixs[i]
			for j=1,self.numStates do
				dq_dts[i][j] = dq_dts[i][j] - (self.fluxes[i+1][j] - self.fluxes[i][j]) / dx
			end
		end

		return dq_dts
	end,
	
	postIterate = function(self, dt)
		self.Ps = self.Ps or {}
		
		local gamma = self.gamma

		for i=1,self.gridsize do
			local rho = self.qs[i][1]
			local u = self.qs[i][2] / rho
			local eTotal = self.qs[i][3] / rho
			local eInt = eTotal - .5 * u * u
			self.Ps[i] = (gamma - 1) * rho * eInt
		end

		--[[ artificial viscosity
		for i=2,self.gridsize-1 do
			local rho = self.qs[i][1]
			local uL = self.qs[i-1][2] / self.qs[i-1][1]
			local uR = self.qs[i+1][2] / self.qs[i+1][1]
			local zeta = 2
			local du = zeta * .5 * (uR - uL)
			self.Ps[i] = self.Ps[i] + du * du * rho
		end
		--]]

		-- diffuse momentum
		for i=2,self.gridsize-1 do
			local dP = self.Ps[i+1] - self.Ps[i-1]
			local dx = self.xs[i+1] - self.xs[i-1]
			self.qs[i][2] = self.qs[i][2] - dP * dt / dx
		end

		-- diffuse work
		for i=2,self.gridsize-1 do
			local WR = self.Ps[i+1] * self.qs[i+1][2] / self.qs[i+1][1]
			local WL = self.Ps[i-1] * self.qs[i-1][2] / self.qs[i-1][1]
			local dW = WR - WL
			local dx = self.xs[i+1] - self.xs[i-1]
			self.qs[i][3] = self.qs[i][3] - dW * dt / dx
		end
	end,
}
	
local HLL = {
	calcDT = function(self, getLeft, getRight)
		-- matches Roe, except without eigenvectors
		for i=2,self.gridsize do
			local qL = getLeft and getLeft(i) or self.qs[i-1]
			local qR = getRight and getRight(i) or self.qs[i]
			self:calcInterfaceEigenvalues(qL, qR, self.eigenvalues[i])
		end

		return calcDtWithEigenvalues(self)
	end,
	calcFlux = function(self, dt)
		local gamma = self.gamma
		
		local iqs = self:newState()
		
		for i=2,self.gridsize do
			-- TODO use qL and qR to allow compatability with MUSCL
			local qL = self.qs[i-1]
			local qR = self.qs[i]
			
			local sL = self.eigenvalues[i][1]
			local sR = self.eigenvalues[i][self.numStates]

			local fluxL = self:calcFluxForState(qL)
			local fluxR = self:calcFluxForState(qR)

			local flux = self.fluxes[i]
			for i=1,self.numStates do
				if 0 <= sL then
					flux[i] = fluxL[i]
				elseif sL <= 0 and 0 <= sR then
					flux[i] = (sR * fluxL[i] - sL * fluxR[i] + sL * sR * (qR[i] - qL[i])) / (sR - sL)
				elseif sR <= 0 then
					flux[i] = fluxR[i]
				end
			end
		end
		
		local dq_dts = self:newState()
		for i=1,self.gridsize do
			local dx = self.ixs[i+1] - self.ixs[i]
			for j=1,self.numStates do
				dq_dts[i][j] = dq_dts[i][j] - (self.fluxes[i+1][j] - self.fluxes[i][j]) / dx
			end
		end
	
		return dq_dts

	end,
}

--local EulerMUSCLParent = HLL
EulerMUSCLParent = Roe

local EulerMUSCL
EulerMUSCL = {
	-- limiter of ratio
	-- popular limiters: Fromm, Beam-Warming, Lax-Wendroff, minmod
	-- notice that winded-ness needs to vary per-scheme
	--slopeLimiter = require 'limiter'.LaxWendroff,
	slopeLimiter = require 'limiter'.Fromm,

	-- first calculate new state interface left and right
	-- then call the old calcDT which also calculates eigen basis stuff
	calcDT = function(self) 

		local sigma = self:newState()
		for i=3,self.gridsize-1 do
			local qL2 = self.qs[i-2]
			local qL = self.qs[i-1]
			local qR = self.qs[i]
			local qR2 = self.qs[i+1]
			for j=1,self.numStates do
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
				
				local phi = EulerMUSCL.slopeLimiter(r)
				
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
		local iqLs = self:newState()
		local iqRs = self:newState()
		for j=1,self.numStates do
			iqLs[1][j] = self.qs[1][j]
			iqRs[1][j] = self.qs[1][j]
		end
		for i=2,self.gridsize do
			local dx = self.ixs[i] - self.ixs[i-1]
			for j=1,self.numStates do
				iqLs[i][j] = self.qs[i-1][j] + .5 * dx * sigma[i-1][j]
				iqRs[i][j] = self.qs[i][j] - .5 * dx * sigma[i][j]
			end
		end

		-- flux based on subgrid values
		local ifLs = self:newState()
		local ifRs = self:newState()
		for i=1,self.gridsize do
			ifLs[i] = self:calcFluxForState(iqLs[i])
			ifRs[i] = self:calcFluxForState(iqRs[i])
		end

		-- how should we get this dt?
		-- based on Roe eigenvector without MUSCL?  or based on Burgers?  or fixed + implicit (optimistically)
		-- should this be the dt that we consistently use even after using MUSCL to adjust the state?
		local dt = EulerMUSCLParent.calcDT(self)

		-- half step in time
		local iqhLs = self:newState()
		local iqhRs = self:newState()
		for j=1,self.numStates do
			iqhLs[1][j] = iqLs[1][j]
			iqhRs[1][j] = iqRs[1][j]
			iqhLs[self.gridsize][j] = iqLs[self.gridsize][j]
			iqhRs[self.gridsize][j] = iqRs[self.gridsize][j]
		end
		for i=2,self.gridsize-1 do
			local dx = self.ixs[i] - self.ixs[i-1]
			for j=1,self.numStates do
				iqhLs[i][j] = iqLs[i][j] + .5 * dt/dx * (ifLs[i][j] - ifRs[i-1][j])
				iqhRs[i][j] = iqRs[i][j] + .5 * dt/dx * (ifLs[i+1][j] - ifRs[i][j])
			end
		end
		
		-- once we have *this* collection of subgrid L & R states,
		--  we use them for whatever method you want ...
		-- ... be it HLL or Roe, etc 

		local getLeft = function(i) return iqhLs[i] end
		local getRight = function(i) return iqhRs[i] end

		--[[
		TODO now the paper has officially gave me a circular dependency:
		to do the MUSCL subgrid half-step in time it needs to know dt
		however, to calculate dt, it needs the eigenvalues
			eigenvalues need Roe matrix at interface
			interface needs the left and right half-step values
			... which need the dt
		--]]
		local dt = EulerMUSCLParent.calcDT(self, getLeft, getRight)
		return dt
	end,

	calcFlux = EulerMUSCLParent.calcFlux,
}

return {
	Roe = Roe,
	-- only works for Eulers
	EulerBurgers = EulerBurgers,
	HLL = HLL,
	EulerMUSCL = EulerMUSCL,
}

