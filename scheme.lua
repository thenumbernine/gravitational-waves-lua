return {
	Roe = function(self)
		-- Roe solver:
		-- 1) calculate eigenbasis at interface with Roe-weighted average of states
		for i=2,self.gridsize do
			self:calcInterfaceEigenBasis(i)

			-- collect error for reporting
			local eigenbasisError = 0
			local fluxMatrixError = 0
			for j=1,self.numStates do
				for k=1,self.numStates do
					local eigenbasisSum = 0
					local fluxMatrixSum = 0
					for l=1,self.numStates do
						eigenbasisSum = eigenbasisSum + self.eigenvectors[i][j][l] * self.eigenvectorsInverse[i][l][k]
						fluxMatrixSum = fluxMatrixSum + self.eigenvectors[i][j][l] * self.eigenvectorsInverse[i][l][k] * self.eigenvalues[i][l]
					end
					eigenbasisError = eigenbasisError + abs(eigenbasisSum - (j == k and 1 or 0))
					fluxMatrixError = fluxMatrixError + abs(fluxMatrixSum - self.fluxMatrix[i][j][k])
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
			for j=1,self.numStates do
				local s = 0
				for k=1,self.numStates do
					s = s + self.eigenvectorsInverse[i][j][k] * (self.qs[i][k] - self.qs[i-1][k])
				end
				self.deltaQTildes[i][j] = s
			end
		end
		
		-- 3) slope limit on interface difference
		-- 4) transform back
		for i=2,self.gridsize do
			local fluxTilde = {}
			for j=1,self.numStates do
				local rTilde
				if self.deltaQTildes[i][j] == 0 then
					rTilde = 0
				else
					if self.eigenvalues[i][j] >= 0 then
						rTilde = self.deltaQTildes[i-1][j] / self.deltaQTildes[i][j]
					else
						rTilde = self.deltaQTildes[i+1][j] / self.deltaQTildes[i][j]
					end
				end
				local phi = self.slopeLimiter(rTilde)
				local theta = self.eigenvalues[i][j] >= 0 and 1 or -1
				local dx = self.xs[i] - self.xs[i-1]
				local epsilon = self.eigenvalues[i][j] * dt / dx
				local deltaFluxTilde = self.eigenvalues[i][j] * self.deltaQTildes[i][j]
				fluxTilde[j] = -.5 * deltaFluxTilde * (theta + phi * (epsilon - theta))
			end
			for j=1,self.numStates do
				local s = 0
				for k=1,self.numStates do
					s = s + self.eigenvectors[i][j][k] * fluxTilde[k]
					-- using the flux matrix itself allows for complete reconstruction even in the presence of zero-self.eigenvalues
					s = s + self.fluxMatrix[i][j][k] * (self.qs[i-1][k] + self.qs[i][k]) * .5
				end
				self.fluxes[i][j] = s
			end
		end

		return dt
	end,
}

