local class = require 'ext.class'
local SolverFV = require 'solverfv'

local Roe = class(SolverFV)

function Roe:init(args)
	Roe.super.init(self, args)

	self.deltaQTildes = {}
	self.rTildes = {}
	self.Phis = {}
	self.fluxMatrix = {}
	self.eigenvectors = {}
	self.eigenvectorsInverse = {}
	self.eigenbasisErrors = {}
	self.fluxMatrixErrors = {}
end

local function buildField(matrixField)
	return function(self, i, v)
		local m = self[matrixField][i]
		local result = {}
		for j=1,self.numStates do
			local sum = 0
			for k=1,self.numStates do
				sum = sum + m[j][k] * v[k]
			end
			result[j] = sum
		end
		return result 
	end
end

--[[
default implementation will dot with j'th row of eigenvectorsInverse[i]
subclasses with sparse matrices (like ADM) will be able to override this and optimize away (those 37x37 matrices)

another note: eigenfields never have input vectors.  they are made of state vaules, and their input is state values, so there's no need to define an inner product.
...except the fact that some of the state variables are on the i'th entry, and some are of the i+1/2'th entry...
--]]
Roe.fluxTransform = buildField'fluxMatrix'
Roe.eigenfields = buildField'eigenvectorsInverse'
Roe.eigenfieldsInverse = buildField'eigenvectors'

function Roe:reset()
	Roe.super.reset(self)

	for i=1,self.gridsize do
		self.rTildes[i] = {}
		for j=1,self.numStates do
			self.rTildes[i][j] = 0
		end
	end

	-- state interfaces
	for i=1,self.gridsize+1 do
		self.deltaQTildes[i] = {}
		for j=1,self.numStates do
			self.deltaQTildes[i][j] = 0
		end
	end

	for i=1,self.gridsize do
		self.Phis[i] = {}
		for j=1,self.numStates do
			self.Phis[i][j] = 1
		end
	end

	for i=1,self.gridsize+1 do
		self.fluxMatrix[i] = {}
		self.eigenvectors[i] = {}
		self.eigenvectorsInverse[i] = {}
		for j=1,self.numStates do
			self.fluxMatrix[i][j] = {}
			self.eigenvectors[i][j] = {}
			self.eigenvectorsInverse[i][j] = {}
		end
		self.eigenbasisErrors[i] = 0
		self.fluxMatrixErrors[i] = 0
	end
end

-- calculates timestep and eigenbasis
function Roe:calcDT(getLeft, getRight)
	-- Roe solver:
	-- 1) calculate eigenbasis at interface with Roe-weighted average of states
	for i=2,self.gridsize do
		local qL = getLeft and getLeft(i) or self.qs[i-1]
		local qR = getRight and getRight(i) or self.qs[i]
		
		self.equation:calcInterfaceEigenBasis(self,i,qL,qR)

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

	return Roe.super.calcDT(self)
end

function Roe:calcDeltaQTildes(getLeft, getRight)
	for i=2,self.gridsize do
		local qL = getLeft and getLeft(i) or self.qs[i-1]
		local qR = getRight and getRight(i) or self.qs[i]
		
		local dq = {}
		for j=1,self.numStates do
			dq[j] = qR[j] - qL[j]
		end
		self.deltaQTildes[i] = self:eigenfields(i, dq)
	end
end

function Roe:calcRTildes(getLeft, getRight)
	for i=2,self.gridsize do
		local qL = getLeft and getLeft(i) or self.qs[i-1]
		local qR = getRight and getRight(i) or self.qs[i]
		for j=1,self.numStates do
			if self.deltaQTildes[i][j] == 0 then
				self.rTildes[i][j] = 0
			else
				if self.eigenvalues[i][j] >= 0 then
					self.rTildes[i][j] = self.deltaQTildes[i-1][j] / self.deltaQTildes[i][j]
				else
					self.rTildes[i][j] = self.deltaQTildes[i+1][j] / self.deltaQTildes[i][j]
				end
			end
		end
	end
end

function Roe:calcPhis(dt)
	-- compute Phi vector along interfaces
	-- Phi = diag(sign(v_i) + 1/2 phi * (dt/dx v_i - sign(v_i)))
	for i=2,self.gridsize do
		local dx = self.xs[i] - self.xs[i-1]
		for j=1,self.numStates do
			local phi = self.fluxLimiter(self.rTildes[i][j])
			local theta = self.eigenvalues[i][j] >= 0 and 1 or -1
			local epsilon = self.eigenvalues[i][j] * dt / dx
			self.Phis[i][j] = theta + phi * (epsilon - theta)
		end
	end
end

function Roe:calcFlux(dt, getLeft, getRight, getLeft2, getRight2)

	-- 2) calculate interface state difference in eigenbasis coordinates
	self:calcDeltaQTildes(getLeft, getRight)

	-- 3) slope limit on interface difference
	self:calcRTildes(getLeft, getRight)

	-- 4) phis based on eigenvalues and flux limiter of rTildes
	self:calcPhis(dt)
	-- TODO make use of Phis below (take from roe_backwardeuler_linear)

	local useFluxMatrix = false
	
	-- 5) transform back
	for i=2,self.gridsize do
		local qL = getLeft and getLeft(i) or self.qs[i-1]
		local qR = getRight and getRight(i) or self.qs[i]

		local qAvg = {}
		for j=1,self.numStates do
			qAvg[j] = .5 * (qR[j] + qL[j])
		end
			
		local fluxTilde = {}
		for j=1,self.numStates do
			local phi = self.fluxLimiter(self.rTildes[i][j])
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
end

return Roe

