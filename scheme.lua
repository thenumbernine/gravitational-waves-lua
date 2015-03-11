return {
	Roe = function(self)
		-- Roe solver:
		-- 1) calculate eigenbasis at interface with Roe-weighted average of states
		for i=2,self.gridsize do
			self:calcInterfaceEigenBasis(i)

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

		-- determine timestep based on eigenvalue interfaces
		local dt
		if self.fixed_dt then
			dt = self.fixed_dt
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
			dt = result * self.cfl
		end

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
				local phi = self.slopeLimiter(rTildes[j])
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

		return dt
	end,
}

