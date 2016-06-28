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

		-- should I be modifying the equation?
		self.equation.applyLeftEigenvectorsAtCenter = self.equation.createTransformFunc('cellLeftEigenvectors', true, false)
		self.equation.applyRightEigenvectorsAtCenter = self.equation.createTransformFunc('cellRightEigenvectors', false, true)

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
			
			-- calc eigenvalues and vectors at cell center
			self.equation:calcEigenBasis(
				lambdas,
				self.cellRightEigenvectors[i],
				self.cellLeftEigenvectors[i],
				nil,	-- not computing dF/dU at the cell center.  I only ever compute this for testing eigenbasis reconstruction accuracy.
				self.equation:calcRoeValuesAtCellCenter(self.qs[i]))
		
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
			local deltaQTildeL = self.equation:applyLeftEigenvectorsAtCenter(self, i, deltaQL)
			local deltaQTildeR = self.equation:applyLeftEigenvectorsAtCenter(self, i, deltaQR)
			local deltaQTildeC = self.equation:applyLeftEigenvectorsAtCenter(self, i, deltaQC)
		
			-- (38) calc delta a^m TVD reconstruction
			local deltaQTildeM = {}
			for j=1,self.numWaves do
				deltaQTildeM[j] = 
					sign(deltaQTildeC[j]) 
					* math.min(
						2 * math.abs(deltaQTildeL[j]), 
						2 * math.abs(deltaQTildeR[j]), 
						math.abs(deltaQTildeC[j]))
			end
		
			-- (39) char -> prim
			local deltaQM = self.equation:applyRightEigenvectorsAtCenter(self, i, deltaQTildeM)
		
			-- (40, 41)
			local dx = self.ixs[i+1] - self.ixs[i]
			local dt_dx = dt / dx
			local lambdaMin = lambdas[1]	-- jacobian of of flux at cell center, min eigenvalue
			local lambdaMax = lambdas[self.numWaves]	-- ... and max eigenvalue
		
			local qHatLs = {}	-- right interface left side
			local qHatRs = {}	-- left interface right side
			for j=1,self.numStates do
				-- \hat{w}_{L,i+1/2} left of interface = right of interface's left cell's center
				qHatLs[j] = self.qs[i][j] + .5 * (1 - math.max(lambdaMax, 0) * dt_dx) * deltaQM[j]
				-- \hat{w}_{R,i-1/2} right of interface = left of interface's right cell's center
				qHatRs[j] = self.qs[i][j] - .5 * (1 - math.min(lambdaMin, 0) * dt_dx) * deltaQM[j]
			end

			-- (42, 43)

			-- delta w^m in characteristic space
			-- wait, if deltaQM are the right apply of deltaQTidleMs, and nothing modified them since then,
			-- then how is delta_wm_char anything other than deltaQTildeMs?
			local delta_wm_char = self.equation:applyLeftEigenvectorsAtCenter(self, i, deltaQM)
			
			-- delta w^m, in char space (original times left eigenvectors), with all eigenvalues <=0 times zero, and the others times the max eigenvalue minus its eigenvalue
			local deltaQTildeM_pos = {}
			-- same as above, but only negative eigenvalues, and times min eigenvalue minus this eigenvalue
			local deltaQTildeM_neg = {}
			for j=1,self.numWaves do
				-- paper says in eqn 44 & 45 that, for HLL only (not Roe), we should be adding *all* nonzero waves to both, and not just positive or negative ones to each
				deltaQTildeM_pos[j] = (lambdas[j] <= 0 and 0 or 1) * (lambdaMax - lambdas[j]) * delta_wm_char[j]
				deltaQTildeM_neg[j] = (lambdas[j] >= 0 and 0 or 1) * (lambdaMin - lambdas[j]) * delta_wm_char[j]
			end
			
			local deltaQM_pos = self.equation:applyRightEigenvectorsAtCenter(self, i, deltaQTildeM_pos)
			local deltaQM_neg = self.equation:applyRightEigenvectorsAtCenter(self, i, deltaQTildeM_neg)
			
			for j=1,self.numStates do
				self.qLs[i+1][j] = qHatLs[j] + .5 * dt_dx * deltaQM_pos[j]
				self.qRs[i][j] = qHatRs[j] + .5 * dt_dx * deltaQM_neg[j]
			end
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
