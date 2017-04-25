local class = require 'ext.class'
local SolverFV = require 'solverfv'

local Roe = class(SolverFV)

function Roe:init(args)
	self.equation = assert(args.equation or self.equation)

	-- add graphs of eigenbasis orthogonality and flux derivative reconstruction errors 
	local graphInfo = {
		name = 'log eigenbasis error',
		getter = function(self,i)
			return self.eigenbasisErrors and math.log(self.eigenbasisErrors[i], 10)
		end,
	}
	self.equation.graphInfos:insert(graphInfo)
	self.equation.graphInfoForNames[graphInfo.name] = graphInfo
	
	local graphInfo = {
		name = 'log reconstruction error',
		getter = function(self,i)
			return self.fluxMatrixErrors and math.log(self.fluxMatrixErrors[i], 10)
		end,
	}
	self.equation.graphInfos:insert(graphInfo)
	self.equation.graphInfoForNames[graphInfo.name] = graphInfo

	self.name = self.equation.name .. ' Roe'
	
	Roe.super.init(self, args)

	self.deltaQTildes = {}
	self.rTildes = {}
	self.fluxMatrix = {}
	self.eigenvalues = {}
	self.rightEigenvectors = {}
	self.leftEigenvectors = {}
	self.eigenbasisErrors = {}
	self.fluxMatrixErrors = {}
end

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

	for i=1,self.gridsize+1 do
		self.eigenvalues[i] = {}
		for j=1,self.numWaves do
			self.eigenvalues[i][j] = 0
		end

		-- #states x #states
		self.fluxMatrix[i] = {}
		for j=1,self.numStates do
			self.fluxMatrix[i][j] = {}
			for k=1,self.numStates do
				self.fluxMatrix[i][j][k] = 0
			end
		end
		-- #waves x #states
		self.leftEigenvectors[i] = {}
		for j=1,self.numWaves do
			self.leftEigenvectors[i][j] = {}
			for k=1,self.numStates do
				self.leftEigenvectors[i][j][k] = 0
			end
		end
		-- #states x #waves
		self.rightEigenvectors[i] = {}
		for j=1,self.numStates do
			self.rightEigenvectors[i][j] = {}
			for k=1,self.numWaves do
				self.rightEigenvectors[i][j][k] = 0
			end
		end
		
		self.eigenbasisErrors[i] = 0
		self.fluxMatrixErrors[i] = 0
	end
end

-- calculates timestep and eigenbasis
-- these are used by calcFluxAtInterface as well
-- so this method counts on the fact that calcDT is called first, before calcFluxAtInterface
-- however, for RK4 integration, shouldn't eigenbasis be recomputed at each intermediate state? 
function Roe:calcInterfaceEigenBasis()
	-- Roe solver:
	-- 1) calculate eigenbasis at interface with Roe-weighted average of states
	for i=2,self.gridsize do
		local qL = self:get_qL(i)
		local qR = self:get_qR(i)
		
		self.equation:calcEigenBasis(
			self.eigenvalues[i],
			self.rightEigenvectors[i],
			self.leftEigenvectors[i],
			self.fluxMatrix[i],
			self.equation:calcInterfaceRoeValues(self, i))

		-- collect error for reporting
		local eigenbasisError = 0
		local fluxMatrixError = 0
		-- for the i'th cell:
		-- Q = right eigenvectors matrix
		-- L = eigenvalues matrix
		-- A_jk = Q_jl L_l invQ_lk
		-- eigenbasisError = Q_jl invQ_lk - delta_jk
		-- fluxMatrixError = Q_jl L_l invQ_lk - A_jk
		for j=1,self.numStates do
			-- local basis_j = 0's everywhere except a 1 at the j'th entry
			-- local eigencoords_j = {k,leftEigenvectors[i][k](basis_j)}			<- dot input vector with eigenvector inverse row k
			-- local eigenscaled_j = eigencoords_j * lambda_j
			-- local newbasis_j = {k,rightEigenvectors[i][k](eigencoords_j)}	<- dot input vector with eigenvector row k
			-- local newtransformed_j = {k,rightEigenvectors[i][k](eigenscaled_j)}
			-- sum up abs error between basis_j and newbasis_j 
			-- sum up abs error between A_jk basis_k and newtransformed_j
		
			-- basis_k = delta_jk
			local basis = {}
			for k=1,self.numStates do
				basis[k] = k == j and 1 or 0
			end

			-- eigenCoords_k = invQ_kl basis_l
			local eigenCoords = self.equation:eigenLeftTransform(self, self.leftEigenvectors[i], basis)
			for k=1,self.numWaves do
				assert(type(eigenCoords[k])=='number', "failed for coord "..k.." got type "..type(eigenCoords[k]))
			end

			-- eigenScaled_k = lambda_k * eigenCoords_k
			local eigenScaled = {}
			for k=1,self.numWaves do
				eigenScaled[k] = self.eigenvalues[i][k] * eigenCoords[k]
			end
			
			-- newbasis_k = Q_kl eigenCoords_l
			local newbasis = self.equation:eigenRightTransform(self, self.rightEigenvectors[i], eigenCoords)
			
			-- newtransformed_k = Q_kl eigenScaled_l = Q_kl lambda_l eigenCoords_k
			local newtransformed = self.equation:eigenRightTransform(self, self.rightEigenvectors[i], eigenScaled)

			-- transformed_k = A_kl basis_l
			local transformed = self.equation:fluxMatrixTransform(self, self.fluxMatrix[i], basis)

			for k=1,self.numStates do
				eigenbasisError = eigenbasisError + math.abs(basis[k] - newbasis[k])
				fluxMatrixError = fluxMatrixError + math.abs(transformed[k] - newtransformed[k])
			end
		end
		
		-- the eigenbasis should only reconstruct to the identity with dimension equal to numWaves
		-- however numStates basis vectors have to be tested to find this
		-- so subtract off the difference
		-- or, if you know specifically what conservative values are zeroed, exclude those from the basis test  ... but that only works if it is a specific conservative variable, and not a linear combination of some
		eigenbasisError = eigenbasisError - (self.numStates - self.numWaves)
		
		self.eigenbasisErrors[i] = eigenbasisError
		self.fluxMatrixErrors[i] = fluxMatrixError
	end
end

-- get the q at the left side of the interface
function Roe:get_qL(i)
	return self.qs[i-1]
end

-- get the q at the right side of the interface
function Roe:get_qR(i)
	return self.qs[i]
end

function Roe:calcDeltaQTildes()
	for i=2,self.gridsize do
		local qL = self:get_qL(i)
		local qR = self:get_qR(i)
		
		local dq = {}
		for j=1,self.numStates do
			dq[j] = qR[j] - qL[j]
		end
		self.deltaQTildes[i] = self.equation:eigenLeftTransform(self, self.leftEigenvectors[i], dq)
	end
end

-- depends on calcDeltaQTildes
function Roe:calcRTildes()
	for i=2,self.gridsize do
		for j=1,self.numWaves do
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

-- calcluate self.fluxes
-- depends on self:calcDeltaQTildes and self:calcRTiles
function Roe:calcFluxAtInterface(dt, i)
	-- if we can calculate the flux directly then use that
	local canCalcFlux = self.equation.calcFluxForState
	
	local qL = self:get_qL(i)
	local qR = self:get_qR(i)
	local lambdas = self.eigenvalues[i]

	--[[ shortcut if the equation has a 'calcFluxForState' function
	-- but this bypasses any influence from the flux limiter ... 
	-- so I don't think I'll use it
	if canCalcFlux then
		if lambdas[1] >= 0 then	-- smallest eigenvalue >= 0
			fill(self.fluxes[i], self.equation:calcFluxForState(qL))
			return	
		end
		if lambdas[self.numWaves] <= 0 then	-- largest eigenvalues <= 0
			fill(self.fluxes[i], self.equation:calcFluxForState(qR))
			return
		end
	end
	--]]

	local fluxTilde = {}
	for j=1,self.numWaves do
		local lambda = lambdas[j]
		local phi = self.fluxLimiter.func(self.rTildes[i][j])
		local sgnLambda = lambda >= 0 and 1 or -1
		local dx = self.xs[i] - self.xs[i-1]
		local epsilon = lambda * dt / dx
		local deltaQTilde = self.deltaQTildes[i][j]
		fluxTilde[j] = -.5 * lambda * deltaQTilde * (sgnLambda + phi * (epsilon - sgnLambda))
	end
	
	if not canCalcFlux then
		local qAvg = {}
		for j=1,self.numStates do
			qAvg[j] = .5 * (qR[j] + qL[j])
		end	
		local qAvgTildes = self.equation:eigenLeftTransform(self, self.leftEigenvectors[i], qAvg)
		for j=1,self.numWaves do
			fluxTilde[j] = fluxTilde[j] + lambdas[j] * qAvgTildes[j]
		end
	end
	
	local flux = self.equation:eigenRightTransform(self, self.rightEigenvectors[i], fluxTilde)
	
	-- using the flux itself allows for complete reconstruction even in the presence of zero-self.eigenvalues
	if canCalcFlux then
		local FL = {self.equation:calcFluxForState(qL)}
		local FR = {self.equation:calcFluxForState(qR)}
		for j=1,self.numStates do
			flux[j] = flux[j] + .5 * (FL[j] + FR[j])
		end
	end

	self.fluxes[i] = flux
end

function Roe:calcFluxes(dt)
	-- used by all following methods in calcFluxes
	self:calcInterfaceEigenBasis()

	-- calculate interface state difference in eigenbasis coordinates
	self:calcDeltaQTildes()

	-- slope limit on interface difference
	self:calcRTildes()

	for i=2,self.gridsize do
		self:calcFluxAtInterface(dt, i)
	end
	
	-- this eliminates the flux within the ghost cells
	-- which is especially useful for the implicit solvers
	-- but that wasn't really causing any trouble anyways
	for j=1,self.numStates do
		self.fluxes[1][j] = self.fluxes[2][j]
		self.fluxes[self.gridsize+1][j] = self.fluxes[self.gridsize][j]
	end
end

return Roe
