--[[
going by this paper:
https://arxiv.org/pdf/0804.0402v1.pdf

this is working for Roe, but not for MHD, ADM, etc
--]]
local class = require 'ext.class'

local function PLMBehavior(parentClass)

	local PLMTemplate = class(parentClass)

	local function sign(x)
		return x == 0 and 0 or (x > 0 and 1 or -1)
	end

	function PLMTemplate:init(args)
		
		-- disable flux limiter
		args = table(args)
		local limiter = require 'limiter' 
		args.fluxLimiter = limiter.donorCell
		
		PLMTemplate.super.init(self, args)
		self.name = self.name .. ' PLM'

		-- cells 
		self.cellLeftEigenvectors = {}
		self.cellRightEigenvectors = {}
		
		-- interface
		self.qLs = {}
		self.qRs = {}
	end

	function PLMTemplate:reset()
		PLMTemplate.super.reset(self)

		-- centered
		for i=1,self.gridsize do
			self.cellLeftEigenvectors[i] = {}
			self.cellRightEigenvectors[i] = {}
			for j=1,self.numStates do
				self.cellLeftEigenvectors[i][j] = {}
				self.cellRightEigenvectors[i][j] = {}
			end
				
			self.qLs[i] = {}
			self.qRs[i] = {}
			for j=1,self.numStates do
				self.qLs[i][j] = 0
				self.qRs[i][j] = 0
			end
		end
	end
	
	function PLMTemplate:calcFluxes(dt)

		for i=2,self.gridsize-1 do
			-- only need to keep one copy of cell eigenvalues, for the current cell we're operating on
			-- I could do the same with eigenvectors, but their application function (which isn't necessarily a matrix multiply)
			--  is currently tied up with storage, passing the index
			local lambdas = {}
			local rightEigenvectors = {}
			local leftEigenvectors = {}
			for j=1,self.numStates do
				rightEigenvectors[j] = {}
				leftEigenvectors[j] = {}
			end
			
			-- calc eigenvalues and vectors at cell center
			self.equation:calcEigenBasis(
				lambdas,
				rightEigenvectors,
				leftEigenvectors,
				nil,	-- not computing dF/dU at the cell center.  I only ever compute this for testing eigenbasis reconstruction accuracy.
				self.equation:calcCellCenterRoeValues(self, i))
		
			-- (36) calc delta q's
			local deltaQL = {}
			local deltaQR = {}
			local deltaQC = {}
			for j=1,self.numStates do
				deltaQL[j] = self.qs[i][j] - self.qs[i-1][j]
				deltaQR[j] = self.qs[i+1][j] - self.qs[i][j]
				deltaQC[j] = self.qs[i+1][j] - self.qs[i-1][j]
			end

			-- (37) calc 'delta qtildes' by 'hydrodynamics ii' / 'delta a's by the paper
			local deltaQTildeL = self.equation:eigenLeftTransform(self, leftEigenvectors, deltaQL)
			local deltaQTildeR = self.equation:eigenLeftTransform(self, leftEigenvectors, deltaQR)
			local deltaQTildeC = self.equation:eigenLeftTransform(self, leftEigenvectors, deltaQC)
		
			-- (38) calc delta a^m TVD reconstruction
			local deltaQTildeM = {}
			for j=1,self.numWaves do
				deltaQTildeM[j] = 
					math.min(
						2 * math.abs(deltaQTildeL[j]), 
						2 * math.abs(deltaQTildeR[j]), 
						math.abs(deltaQTildeC[j]))
					* sign(deltaQTildeC[j]) 
					-- I didn't see this in the paper ... but it's in the code 
					-- this looks like a 2nd deriv constraint to me:
					* math.max(sign(deltaQTildeL[j] * deltaQTildeR[j]), 0)
			end

			local dx = self.ixs[i+1] - self.ixs[i]
			local dt_dx = dt / dx

			local pl = {}
			local pr = {}
			for j=1,self.numWaves do
				pl[j] = deltaQTildeM[j] * .5 * (lambdas[j] >= 0 and 1 - dt_dx * lambdas[j] or 0)
				pr[j] = deltaQTildeM[j] * .5 * (lambdas[j] <= 0 and 1 + dt_dx * lambdas[j] or 0)
			end

			local qTilde = self.equation:eigenLeftTransform(self, leftEigenvectors, self.qs[i])
			local qp = {}
			local qm = {}
			for j=1,self.numWaves do
				qp[j] = qTilde[j] + pl[j]
				qm[j] = qTilde[j] - pr[j]
			end
			self.qLs[i+1] = self.equation:eigenRightTransform(self, rightEigenvectors, qp)
			self.qRs[i] = self.equation:eigenRightTransform(self, rightEigenvectors, qm)
--[[
print('q', tolua(self.qs[i]))
local qTilde = self.equation:eigenLeftTransform(self, leftEigenvectors, self.qs[i])
print('l', tolua(leftEigenvectors))
print('l q', tolua(qTilde))
local q = self.equation:eigenRightTransform(self, rightEigenvectors, qTilde)
print('r', tolua(rightEigenvectors))
print('r l q', tolua(q))
self.qLs[i+1] = q 
self.qRs[i] = q 
--]]
		end
		
		-- now qLs and qRs can be used
		PLMTemplate.super.calcFluxes(self, dt)
	end

	function PLMTemplate:get_qL(i)
		for j=1,3 do
			assert(self.qLs[i][j], "failed for i,j="..i..','..j)
		end
		return self.qLs[i]
	end

	function PLMTemplate:get_qR(i)
		for j=1,3 do
			assert(self.qRs[i][j], "failed for i,j="..i..','..j)
		end
		return self.qRs[i]
	end

	return PLMTemplate
end

return PLMBehavior
