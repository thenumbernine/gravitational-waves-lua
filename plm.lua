-- going by this paper:
-- https://arxiv.org/pdf/0804.0402v1.pdf
local class = require 'ext.class'

local function PLMBehavior(parentClass)

	local PLM = class(parentClass)

	function PLM:init(args)
		self.equation = assert(args.equation or self.equation) 
		self.name = self.equation.name .. ' PLM'
		PLM.super.init(self, args)

		-- cells 
		self.lambdaCs = {}
		self.evLCWrtPrims = {}
		self.evRCWrtPrims = {}
		self.ws = {}
		self.delta_wLs = {}
		self.delta_wRs = {}
		self.delta_wCs = {}
		self.delta_aLs = {}
		self.delta_aRs = {}
		self.delta_aCs = {}
		self.delta_ams = {}
		self.delta_wms = {}
	
		-- interface
		self.wHatLs = {}
		self.wHatRs = {}
		self.wLs = {}
		self.wRs = {}
		self.qLs = {}
		self.qRs = {}
	end

	function PLM:reset()
		PLM.super.reset(self)

		-- centered
		for i=1,self.gridsize do
			self.lambdaCs[i] = {}
			self.evLCWrtPrims[i] = {}
			self.evRCWrtPrims[i] = {}
			for j=1,self.numStates do
				self.evLCWrtPrims[i][j] = {}
				self.evRCWrtPrims[i][j] = {}
			end
				
			self.ws[i] = {}	
			self.delta_wLs[i] = {}
			self.delta_wRs[i] = {}
			self.delta_wCs[i] = {}
			self.delta_aLs[i] = {}
			self.delta_aRs[i] = {}
			self.delta_aCs[i] = {}
			self.delta_ams[i] = {}
			self.delta_wms[i] = {}
		end
	
		for i=1,self.gridsize+1 do
			self.qLs[i] = {}
			self.qRs[i] = {}
			self.wHatLs[i] = {}
			self.wHatRs[i] = {}
			self.wLs[i] = {}
			self.wRs[i] = {}
		end
	end

	local function sign(x)
		return x == 0 and 0 or (x > 0 and 1 or -1)
	end

	local solverLinearFunc = require 'solverlinearfunc'
	PLM.applyLeftEigenvectorsWrtPrimsAtCenter = solverLinearFunc'evLCWrtPrims'
	PLM.applyRightEigenvectorsWrtPrimsAtCenter = solverLinearFunc'evRCWrtPrims'

	function PLM:calcFlux(dt)

		-- calc w's
		for i=1,self.gridsize do
			fill(self.ws[i], self.equation:calcPrimFromCons(table.unpack(self.qs[i])))
		end

		-- calc eigenvalues and vectors at cell center
		for i=1,self.gridsize do
			local rho, vx, P = table.unpack(self.ws[i])
			local ETotal = self.qs[i][3]
			local hTotal = self.equation:calc_hTotal(rho, P, ETotal)
			local lambda = self.lambdaCs[i]
			local evRWrtPrim = self.evRCWrtPrims[i]
			local evLWrtPrim = self.evLCWrtPrims[i]
			local evR, lambda, evL = self.equation:calcEigenBasisWrtPrims(rho, vx, hTotal, nil, nil, lambda, evLWrtPrim, evRWrtPrim)
		end

		-- (36) calc delta w's
		for i=2,self.gridsize-1 do
			for j=1,self.numStates do
				self.delta_wLs[i][j] = self.ws[i][j] - self.ws[i-1][j]
				self.delta_wRs[i][j] = self.ws[i+1][j] - self.ws[i][j]
				self.delta_wCs[i][j] = self.ws[i+1][j] - self.ws[i-1][j]
			end
		end

		-- (37) calc 'delta qtildes' by 'hydrodynamics ii' / 'delta a's by the paper
		for i=2,self.gridsize-1 do
			self.delta_aLs[i] = self:applyLeftEigenvectorsWrtPrimsAtCenter(i, self.delta_wLs[i])
			self.delta_aRs[i] = self:applyLeftEigenvectorsWrtPrimsAtCenter(i, self.delta_wRs[i])
			self.delta_aCs[i] = self:applyLeftEigenvectorsWrtPrimsAtCenter(i, self.delta_wCs[i])
		end

		-- (38) calc delta a^m TVD reconstruction
		for i=2,self.gridsize-1 do
			for j=1,self.numStates do
				self.delta_ams[i][j] = sign(self.delta_aCs[i][j]) * math.min(2 * math.abs(self.delta_aLs[i][j]), 2 * math.abs(self.delta_aRs[i][j]), math.abs(self.delta_aCs[i][j]))
			end
		end

		-- (39) char -> prim
		for i=2,self.gridsize-1 do
			self.delta_wms[i] = self:applyRightEigenvectorsWrtPrimsAtCenter(i, self.delta_ams[i])
		end

		-- (40, 41)
		for i=2,self.gridsize-1 do
			local dx = self.ixs[i+1] - self.ixs[i]
			local dt_dx = dt / dx
			local lambda = self.lambdaCs[i]
			local lambdaMin = lambda[1]	-- jacobian of of flux at cell center, min eigenvalue
			local lambdaMax = lambda[self.numStates]	-- ... and max eigenvalue
			for j=1,self.numStates do
				-- \hat{w}_{L,i+1/2} left of interface = right of interface's left cell's center
				self.wHatLs[i+1][j] = self.ws[i][j] + .5 * (1 - math.max(lambdaMax, 0) * dt_dx) * self.delta_wms[i][j]
				-- \hat{w}_{R,i-1/2} right of interface = left of interface's right cell's center
				self.wHatRs[i][j] = self.ws[i][j] - .5 * (1 - math.min(lambdaMin, 0) * dt_dx) * self.delta_wms[i][j]
			end
		end

		-- (42, 43)
		for i=2,self.gridsize-1 do
			local dx = self.ixs[i+1] - self.ixs[i]
			local dt_dx = dt / dx

			local lambda = self.lambdaCs[i]
			local lambdaMin = lambda[1]	-- jacobian of of flux at cell center, min eigenvalue
			local lambdaMax = lambda[self.numStates]	-- ... and max eigenvalue

			-- delta w^m in characteristic space
			local delta_wm_char = self:applyLeftEigenvectorsWrtPrimsAtCenter(i, self.delta_wms[i])
			
			-- delta w^m, in char space (original times left eigenvectors), with all eigenvalues <=0 times zero, and the others times the max eigenvalue minus its eigenvalue
			local delta_wm_char_pos = {}
			-- same as above, but only negative eigenvalues, and times min eigenvalue minus this eigenvalue
			local delta_wm_char_neg = {}
			for j=1,self.numStates do
				-- paper says in eqn 44 & 45 that, for HLL only (not Roe), we should be adding *all* nonzero waves to both, and not just positive or negative ones to each
				delta_wm_char_pos[j] = (lambda[j] <= 0 and 0 or 1) * (lambdaMax - lambda[j]) * delta_wm_char[j]
				delta_wm_char_neg[j] = (lambda[j] >= 0 and 0 or 1) * (lambdaMin - lambda[j]) * delta_wm_char[j]
			end
			
			local delta_wm_pos = self:applyRightEigenvectorsWrtPrimsAtCenter(i, delta_wm_char_pos)
			local delta_wm_neg = self:applyRightEigenvectorsWrtPrimsAtCenter(i, delta_wm_char_neg)
			
			for j=1,self.numStates do
				self.wLs[i+1][j] = self.wHatLs[i+1][j] + .5 * dt_dx * delta_wm_pos[j]
				self.wRs[i][j] = self.wHatRs[i][j] + .5 * dt_dx * delta_wm_neg[j]
			end
		end

		-- step 8
		for i=2,self.gridsize-1 do
			self.qLs[i+1] = {self.equation:calcPrimFromCons(table.unpack(self.wLs[i+1]))}
			self.qRs[i] = {self.equation:calcPrimFromCons(table.unpack(self.wRs[i]))}
		end

		-- now qLs and qRs can be used
		PLM.super.calcFlux(self, dt)
	end

	function PLM:get_qL(i)
		return self.qLs[i]
	end

	function PLM:get_qR(i)
		return self.qRs[i]
	end

	return PLM
end

return PLMBehavior
