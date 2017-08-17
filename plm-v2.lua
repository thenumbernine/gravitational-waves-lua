-- going by this paper:
-- https://arxiv.org/pdf/0804.0402v1.pdf
local class = require 'ext.class'

local function PLMBehavior(parentClass)

	local PLMTemplate = class(parentClass)

	local function sign(x)
		return x == 0 and 0 or (x > 0 and 1 or -1)
	end

	function PLMTemplate:init(args)
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
			local deltaQG = {}
			for j=1,self.numStates do
				deltaQL[j] = self.qs[i][j] - self.qs[i-1][j]
				deltaQR[j] = self.qs[i+1][j] - self.qs[i][j]
				deltaQC[j] = self.qs[i+1][j] - self.qs[i-1][j]
				deltaQG[j] = math.max(0, 2 * deltaQL[j] * deltaQR[j] / (deltaQL[j] + deltaQR[j]))
			end
		
			-- (37) calc 'delta qtildes' by 'hydrodynamics ii' / 'delta a's by the paper
			local deltaQTildeL = self.equation:eigenLeftTransform(self, leftEigenvectors, deltaQL)
			local deltaQTildeR = self.equation:eigenLeftTransform(self, leftEigenvectors, deltaQR)
			local deltaQTildeC = self.equation:eigenLeftTransform(self, leftEigenvectors, deltaQC)
			local deltaQTildeG = self.equation:eigenLeftTransform(self, leftEigenvectors, deltaQG)
		
			-- (38) calc delta a^m TVD reconstruction
			local deltaQTildeM = {}
			for j=1,self.numWaves do
				--[[ athena paper and ravi's code
				deltaQTildeM[j] = 
					sign(deltaQTildeC[j]) 
					* math.min(
						2 * math.abs(deltaQTildeL[j]), 
						2 * math.abs(deltaQTildeR[j]), 
						math.abs(deltaQTildeC[j]))
				--]]
				-- [[ athena's code
				deltaQTildeM[j] = deltaQTildeL[j]*deltaQTildeR[j] <= 0 and 0 or
					sign(deltaQTildeC[j])
					* math.min(
						2 * math.abs(deltaQTildeL[j]), 
						2 * math.abs(deltaQTildeR[j]),
						.5 * math.abs(deltaQTildeC[j]),
						math.abs(deltaQTildeG[j]))
				--]]
			end
		
			-- (39) char -> prim
			local deltaQM = self.equation:eigenRightTransform(self, rightEigenvectors, deltaQTildeM)
			
			local dx = self.ixs[i+1] - self.ixs[i]
			local dt_dx = dt / dx
			local lambdaMin = lambdas[1]	-- jacobian of of flux at cell center, min eigenvalue
			local lambdaMax = lambdas[self.numWaves]	-- ... and max eigenvalue
			
			-- [[ the paper
			-- (40, 41)
			local omg = 1 - math.max(lambdaMax, 0) * dt_dx
			local opg = 1 - math.min(lambdaMin, 0) * dt_dx
			--]]
			--[[ the code
			local opg = 1
			local omg = 1
			--]]
		
			-- naming of these matches the paper, opposite the code
			local qHatLs = {}	-- right interface left side
			local qHatRs = {}	-- left interface right side 
			for j=1,self.numStates do
				qHatLs[j] = self.qs[i][j] + .5 * deltaQM[j] * omg
				qHatRs[j] = self.qs[i][j] - .5 * deltaQM[j] * opg
		
				--[[ the code
				local C = qHatLs[j] + omg/opg * qHatRs[j]
				qHatRs[j] = math.max(math.min(self.qs[i][j], self.qs[i-1][j]), qHatRs[j])
				qHatRs[j] = math.min(math.max(self.qs[i][j], self.qs[i-1][j]), qHatRs[j])
				qHatLs[j] = C - qHatRs[j] * omg/opg
				
				qHatLs[j] = math.max(math.min(self.qs[i][j], self.qs[i+1][j]), qHatLs[j])
				qHatLs[j] = math.min(math.max(self.qs[i][j], self.qs[i+1][j]), qHatLs[j])
				qHatRs[j] = C * opg/omg - qHatLs[j] * opg/omg
				--]]
			end


			--[[ the code says to ignore the rest unless it is the ctu integrator ... which I don't think the paper specified was specific to the ctu integrator 
			self.qLs[i+1] = qHatRs
			self.qRs[i] = qHatLs
			--]]
			-- [[ the code
			local deltaQ = {}
			for j=1,self.numStates do
				deltaQ[j] = qHatRs[j] - qHatLs[j]
			end
			local deltaQChar = self.equation:eigenLeftTransform(self, leftEigenvectors, deltaQ)

			local qx = .5 * math.max(lambdaMax, 0) * dt_dx
			for j=1,self.numStates do
				self.qLs[i+1][j] = qHatRs[j] - qx * deltaQ[j]
			end

			local qx = -.5 * math.min(lambdaMin, 0) * dt_dx
			for j=1,self.numStates do
				self.qRs[i][j] = qHatLs[j] + qx * deltaQ[j]
			end
		
			--[=[
			local dQxs = {}
			local deltaQCharLs = {}
			local deltaQCharRs = {}
			for j=1,self.numWaves do
				local lambda = lambdas[j]
				deltaQCharLs[j] = 0
				deltaQCharRs[j] = 0
				if lambda >= 0 then
					deltaQCharLs[j] = deltaQCharLs[j] + (lambdaMax - lambda) * deltaQChar[j]
					deltaQCharRs[j] = deltaQCharRs[j] + (lambdaMin - lambda) * deltaQChar[j]
				end
				if lambda <= 0 then
					deltaQCharLs[j] = deltaQCharLs[j] + (lambdaMax - lambda) * deltaQChar[j]
					deltaQCharRs[j] = deltaQCharRs[j] + (lambdaMin - lambda) * deltaQChar[j]
				end
			end
			local qaqLs = self.equation:eigenRightTransform(self, rightEigenvectors, deltaQCharLs)
			local qaqRs = self.equation:eigenRightTransform(self, rightEigenvectors, deltaQCharRs)
			for j=1,self.numWaves do
				self.qLs[i+1][j] = self.qLs[i+1][j] + .5 * dt_dx * qaqLs[j]
				self.qRs[i][j] = self.qRs[i][j] + .5 * dt_dx * qaqRs[j]
			end
			--]=]
			--]]
		
			-- [[
			-- (42, 43)
			-- delta w^m, in char space (original times left eigenvectors), with all eigenvalues <=0 times zero, and the others times the max eigenvalue minus its eigenvalue
			local deltaQTildeM_pos = {}
			local deltaQTildeM_neg = {}
			for j=1,self.numWaves do
				-- paper says in eqn 44 & 45 that, for HLL only (not Roe), we should be adding *all* nonzero waves to both, and not just positive or negative ones to each
				deltaQTildeM_pos[j] = (lambdas[j] <= 0 and 0 or 1) * (lambdaMax - lambdas[j]) * deltaQChar[j]
				deltaQTildeM_neg[j] = (lambdas[j] >= 0 and 0 or 1) * (lambdaMin - lambdas[j]) * deltaQChar[j]
			end
			
			local deltaQM_pos = self.equation:eigenRightTransform(self, rightEigenvectors, deltaQTildeM_pos)
			local deltaQM_neg = self.equation:eigenRightTransform(self, rightEigenvectors, deltaQTildeM_neg)
			
			for j=1,self.numStates do
				self.qLs[i+1][j] = qHatLs[j] + .5 * dt_dx * deltaQM_pos[j]
				self.qRs[i][j] = qHatRs[j] + .5 * dt_dx * deltaQM_neg[j]
			end
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
