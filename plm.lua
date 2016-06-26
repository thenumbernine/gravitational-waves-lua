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
		self.cellEigenvalues = {}
		self.cellLeftEigenvectors = {}
		self.cellRightEigenvectors = {}
		self.ws = {}
		self.deltaQLs = {}
		self.deltaQRs = {}
		self.deltaQCs = {}
		self.deltaQTildeLs = {}
		self.deltaQTildeRs = {}
		self.deltaQTildeCs = {}
		self.deltaQTildeMs = {}
		self.deltaQMs = {}
	
		-- interface
		self.qHatLs = {}
		self.qHatRs = {}
		self.wLs = {}
		self.wRs = {}
		self.qLs = {}
		self.qRs = {}
	end

	function PLMTemplate:reset()
		PLMTemplate.super.reset(self)

		-- centered
		for i=1,self.gridsize do
			self.cellEigenvalues[i] = {}
			self.cellLeftEigenvectors[i] = {}
			self.cellRightEigenvectors[i] = {}
			for j=1,self.numStates do
				self.cellLeftEigenvectors[i][j] = {}
				self.cellRightEigenvectors[i][j] = {}
			end
				
			self.ws[i] = {}	
			self.deltaQLs[i] = {}
			self.deltaQRs[i] = {}
			self.deltaQCs[i] = {}
			self.deltaQTildeLs[i] = {}
			self.deltaQTildeRs[i] = {}
			self.deltaQTildeCs[i] = {}
			self.deltaQTildeMs[i] = {}
			self.deltaQMs[i] = {}

			self.qLs[i] = {}
			self.qRs[i] = {}
			self.wLs[i] = {}
			self.wRs[i] = {}
			for j=1,self.numStates do
				self.qLs[i][j] = 0
				self.qRs[i][j] = 0
				self.wLs[i][j] = 0
				self.wRs[i][j] = 0
			end
		end
		for i=1,self.gridsize+1 do
			self.qHatLs[i] = {}
			self.qHatRs[i] = {}
			for j=1,self.numStates do
				self.qHatLs[i][j] = 0
				self.qHatRs[i][j] = 0
			end
		end
	end
	
	function PLMTemplate:calcFluxes(dt)

		-- calc eigenvalues and vectors at cell center
		for i=1,self.gridsize do
			self.equation:calcEigenBasis(
				self.cellEigenvalues[i],
				self.cellRightEigenvectors[i],
				self.cellLeftEigenvectors[i],
				nil,	-- not computing dF/dU at the cell center.  I only ever compute this for testing eigenbasis reconstruction accuracy.
				self.equation:calcRoeValuesAtCellCenter(self.qs[i]))
		end
		
		-- (36) calc delta q's
		for i=2,self.gridsize-1 do
			for j=1,self.numStates do
				self.deltaQLs[i][j] = self.qs[i][j] - self.qs[i-1][j]
				self.deltaQRs[i][j] = self.qs[i+1][j] - self.qs[i][j]
				self.deltaQCs[i][j] = self.qs[i+1][j] - self.qs[i-1][j]
			end
		end
		
		-- (37) calc 'delta qtildes' by 'hydrodynamics ii' / 'delta a's by the paper
		for i=2,self.gridsize-1 do
			self.deltaQTildeLs[i] = self.equation:applyLeftEigenvectorsAtCenter(self, i, self.deltaQLs[i])
			self.deltaQTildeRs[i] = self.equation:applyLeftEigenvectorsAtCenter(self, i, self.deltaQRs[i])
			self.deltaQTildeCs[i] = self.equation:applyLeftEigenvectorsAtCenter(self, i, self.deltaQCs[i])
		end
		
		-- (38) calc delta a^m TVD reconstruction
		for i=2,self.gridsize-1 do
			for j=1,self.numWaves do
				self.deltaQTildeMs[i][j] = 
					sign(self.deltaQTildeCs[i][j]) 
					* math.min(
						2 * math.abs(self.deltaQTildeLs[i][j]), 
						2 * math.abs(self.deltaQTildeRs[i][j]), 
						math.abs(self.deltaQTildeCs[i][j]))
			end
		end
		
		-- (39) char -> prim
		for i=2,self.gridsize-1 do
			self.deltaQMs[i] = self.equation:applyRightEigenvectorsAtCenter(self, i, self.deltaQTildeMs[i])
		end
		
		-- (40, 41)
		for i=2,self.gridsize-1 do
			local dx = self.ixs[i+1] - self.ixs[i]
			local dt_dx = dt / dx
			local lambda = self.cellEigenvalues[i]
			local lambdaMin = lambda[1]	-- jacobian of of flux at cell center, min eigenvalue
			local lambdaMax = lambda[self.numWaves]	-- ... and max eigenvalue
			for j=1,self.numStates do
				-- \hat{w}_{L,i+1/2} left of interface = right of interface's left cell's center
				self.qHatLs[i+1][j] = self.qs[i][j] + .5 * (1 - math.max(lambdaMax, 0) * dt_dx) * self.deltaQMs[i][j]
				-- \hat{w}_{R,i-1/2} right of interface = left of interface's right cell's center
				self.qHatRs[i][j] = self.qs[i][j] - .5 * (1 - math.min(lambdaMin, 0) * dt_dx) * self.deltaQMs[i][j]
			end
		end

		-- (42, 43)
		for i=2,self.gridsize-1 do
			local dx = self.ixs[i+1] - self.ixs[i]
			local dt_dx = dt / dx

			local lambda = self.cellEigenvalues[i]
			local lambdaMin = lambda[1]	-- jacobian of of flux at cell center, min eigenvalue
			local lambdaMax = lambda[self.numWaves]	-- ... and max eigenvalue

			-- delta w^m in characteristic space
			local delta_wm_char = self.equation:applyLeftEigenvectorsAtCenter(self, i, self.deltaQMs[i])
			
			-- delta w^m, in char space (original times left eigenvectors), with all eigenvalues <=0 times zero, and the others times the max eigenvalue minus its eigenvalue
			local deltaQTildeM_pos = {}
			-- same as above, but only negative eigenvalues, and times min eigenvalue minus this eigenvalue
			local deltaQTildeM_neg = {}
			for j=1,self.numWaves do
				-- paper says in eqn 44 & 45 that, for HLL only (not Roe), we should be adding *all* nonzero waves to both, and not just positive or negative ones to each
				deltaQTildeM_pos[j] = (lambda[j] <= 0 and 0 or 1) * (lambdaMax - lambda[j]) * delta_wm_char[j]
				deltaQTildeM_neg[j] = (lambda[j] >= 0 and 0 or 1) * (lambdaMin - lambda[j]) * delta_wm_char[j]
			end
			
			local deltaQM_pos = self.equation:applyRightEigenvectorsAtCenter(self, i, deltaQTildeM_pos)
			local deltaQM_neg = self.equation:applyRightEigenvectorsAtCenter(self, i, deltaQTildeM_neg)
			
			for j=1,self.numStates do
				self.qLs[i+1][j] = self.qHatLs[i+1][j] + .5 * dt_dx * deltaQM_pos[j]
				self.qRs[i][j] = self.qHatRs[i][j] + .5 * dt_dx * deltaQM_neg[j]
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
