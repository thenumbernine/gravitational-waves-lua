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
			local eigenCoords = self.equation:eigenfields(self, i, basis)
			for k=1,self.numStates do
				assert(type(eigenCoords[k])=='number', "failed for coord "..k.." got type "..type(eigenCoords[k]))
			end

			-- eigenScaled_k = lambda_k * eigenCoords_k
			local eigenScaled = {}
			for k=1,self.numStates do
				eigenScaled[k] = self.eigenvalues[i][k] * eigenCoords[k]
			end
			
			-- newbasis_k = Q_kl eigenCoords_l
			local newbasis = self.equation:eigenfieldsInverse(self, i, eigenCoords)
			
			-- newtransformed_k = Q_kl eigenScaled_l = Q_kl lambda_l eigenCoords_k
			local newtransformed = self.equation:eigenfieldsInverse(self, i, eigenScaled)

			-- transformed_k = A_kl basis_l
			local transformed = self.equation:fluxTransform(self, i, basis)

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
		self.deltaQTildes[i] = self.equation:eigenfields(self, i, dq)
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

--[[
compute Phi vector along interfaces
Phi = diag(sign(v_i) + 1/2 phi * (dt/dx v_i - sign(v_i)))
based on self.rTildes and self.eigenvalues
--]]
function Roe:calcPhis(dt)
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

--[[
calc interface flux vector influence from left/right states
not used by Roe (for efficiency's sake), but used by subclasses

flux = ALeft * qs[i-1] + ARight * qs[i]
for matrixes ALeft & ARight

i = interface index
qs = left/right q vector
dir = direction sign (1 = from left, -1 = from right)
--]]
local zero
function Roe:calcCellFluxForSide(i, qs, dir)
	if not zero then
		zero = {}
		for j=1,self.numStates do
			zero[j] = 0
		end
	end
	if i == 1 or i == self.gridsize+1 then return zero end
	local qTilde = self.equation:eigenfields(self, i, qs)
	local fluxTilde = {}
	for j=1,self.numStates do
		fluxTilde[j] = .5 * self.eigenvalues[i][j] * qTilde[j] * (1 + dir * self.Phis[i][j])
	end
	return self.equation:eigenfieldsInverse(self, i, fluxTilde)
end

--[[
calc interface flux vector
getLeft and getRight would be convenient for applying this to MUSCL
but I would also like an arbitrary state vector ...
combining the two concepts means making qs a parameter of getLeft/getRight ...
not used by Roe (for efficiency's sake), but used by subclasses
--]]
function Roe:calcCellFlux(i, getLeft, getRight)
	local qL = getLeft and getLeft(i) or self.qs[i-1]
	local qR = getRight and getRight(i) or self.qs[i]
	local fluxFromL = self:calcCellFluxForSide(i, qL, 1)
	local fluxFromR = self:calcCellFluxForSide(i, qR, -1)
	local result = {}
	for i=1,self.numStates do
		result[i] = fluxFromL[i] + fluxFromR[i]
	end
	return result
end

--[[
calc derivative of state
not used by Roe (for efficiency's sake), but used by subclasses
--]]
function Roe:calcDeriv(getLeft, getRight)
	local dq_dt = self:newState()
	for i=1,self.gridsize do
		local dx = self.ixs[i+1] - self.ixs[i]
		local fluxL = self:calcCellFlux(i, getLeft, getRight)
		local fluxR = self:calcCellFlux(i+1, getLeft, getRight)
		for j=1,self.numStates do
			dq_dt[i][j] = (fluxL[j] - fluxR[j]) / dx
		end		
	end
	return dq_dt
end

function Roe:calcFlux(dt, getLeft, getRight, getLeft2, getRight2)

	-- 2) calculate interface state difference in eigenbasis coordinates
	self:calcDeltaQTildes(getLeft, getRight)

	-- 3) slope limit on interface difference
	self:calcRTildes(getLeft, getRight)

	-- 4) phis based on eigenvalues and flux limiter of rTildes
	self:calcPhis(dt)
	-- TODO make use of Phis below (take from roe_backwardeuler_linear)

--[=[ uses fluxes (and optionally the flux matrix)
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
			local qAvgTildes = self.equation:eigenfields(self, i, qAvg)
			for j=1,self.numStates do
				fluxTilde[j] = fluxTilde[j] + self.eigenvalues[i][j] * qAvgTildes[j]
			end
		end
		
		self.fluxes[i] = self.equation:eigenfieldsInverse(self, i, fluxTilde)
		
		-- using the flux matrix itself allows for complete reconstruction even in the presence of zero-self.eigenvalues
		if useFluxMatrix then
			local fluxQs = self.equation:fluxTransform(self, i, qAvg)
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
--]=]
-- [=[ doesn't use the flux vector (extra calculatiosn) but uses the same routine that the implicit solver uses
	return self:calcDeriv(getLeft, getRight)
--]=]
end

return Roe
