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
		
		-- interface
		self.qLs = {}
		self.qRs = {}
	end

	function PLMTemplate:reset()
		PLMTemplate.super.reset(self)

		-- centered
		for i=1,self.gridsize do
			self.qLs[i] = {}
			self.qRs[i] = {}
			for j=1,self.numStates do
				self.qLs[i][j] = 0
				self.qRs[i][j] = 0
			end
		end
	end
	
	function PLMTemplate:get_qL(i)
		return self.qLs[i]
	end

	function PLMTemplate:get_qR(i)
		return self.qRs[i]
	end

	function PLMTemplate:calcFluxes(dt)
		for i=2,self.gridsize-1 do
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
		
			-- calculate deltas in conserved variable space wrt grid index 
			local deltaQL = {}
			local deltaQR = {}
			local deltaQC = {}
			for j=1,self.numStates do
				deltaQL[j] = self.qs[i][j] - self.qs[i-1][j]
				deltaQR[j] = self.qs[i+1][j] - self.qs[i][j]
				deltaQC[j] = .5 * (self.qs[i+1][j] - self.qs[i-1][j])
			end

			-- calc left right and centered deltas in characteristic variable space
			local deltaQTildeL = self.equation:eigenLeftTransform(self, leftEigenvectors, deltaQL)
			local deltaQTildeR = self.equation:eigenLeftTransform(self, leftEigenvectors, deltaQR)
			local deltaQTildeC = self.equation:eigenLeftTransform(self, leftEigenvectors, deltaQC)
		
			-- apply slope limiter 
			local deltaQTildeM = {}
			for j=1,self.numWaves do
				deltaQTildeM[j] = 
					2 * math.min(
						math.abs(deltaQTildeL[j] / deltaQTildeC[j]), 
						math.abs(deltaQTildeR[j] / deltaQTildeC[j]), 
						1)
					* deltaQTildeC[j]
					-- if dQL * dQR < 0 then 0 else ...
					* math.max(sign(deltaQTildeL[j] * deltaQTildeR[j]), 0)
			end
			-- deltaQTildeM is the post-limiter slope in characteristic space
			
			local dx = self.ixs[i+1] - self.ixs[i]
			local dt_dx = dt / dx

--[=[ my original plm, which is working just fine

			-- calculate left and right slopes in characteristic space
			local pl = {}
			local pr = {}
			for j=1,self.numWaves do
				pl[j] = lambdas[j] < 0 and 0 or (deltaQTildeM[j] * (1 - dt_dx * lambdas[j]))
				pr[j] = lambdas[j] > 0 and 0 or (deltaQTildeM[j] * (1 + dt_dx * lambdas[j]))
			end

			-- transform slopes back to conserved variable space
			local ql = self.equation:eigenRightTransform(self, rightEigenvectors, pl)
			local qr = self.equation:eigenRightTransform(self, rightEigenvectors, pr)

			-- linearly extrapolate the slopes forward and backward from the cell center
			for j=1,self.numStates do
				self.qLs[i+1][j] = self.qs[i][j] + .5 * ql[j]
				self.qRs[i][j] = self.qs[i][j] - .5 * qr[j]
			end
		
--]=]
-- [=[ subtract off reference states
			local lambdaMin = math.min(0, lambdas[1])
			local lambdaMax = math.max(0, lambdas[#lambdas])

			local deltaQM = self.equation:eigenRightTransform(self, rightEigenvectors, deltaQTildeM)
			local qMinusRef = {}
			local qPlusRef = {}
			for j=1,self.numStates do
				qMinusRef[j] = self.qs[i][j] + .5 * (1 - dt_dx * lambdaMax) * deltaQM[j]	-- left side of the interface on the right of the cell
				qPlusRef[j] = self.qs[i][j] - .5 * (1 + dt_dx * lambdaMin) * deltaQM[j]	-- right side of the interface on the left of the cell
			end

			-- calculate left and right slopes in characteristic space
			local pl = {}
			local pr = {}
			for j=1,self.numWaves do
				pl[j] = lambdas[j] < 0 and 0 or (deltaQTildeM[j] * dt_dx * (lambdaMax - lambdas[j]))
				pr[j] = lambdas[j] > 0 and 0 or (deltaQTildeM[j] * dt_dx * (lambdaMin - lambdas[j]))
			end

			-- transform slopes back to conserved variable space
			local ql = self.equation:eigenRightTransform(self, rightEigenvectors, pl)
			local qr = self.equation:eigenRightTransform(self, rightEigenvectors, pr)

			-- linearly extrapolate the slopes forward and backward from the cell center
			for j=1,self.numStates do
				self.qLs[i+1][j] = qMinusRef[j] + .5 * ql[j]
				self.qRs[i][j] = qPlusRef[j] + .5 * qr[j]
			end
--]=]
		
		end
		
		-- now qLs and qRs can be used
		PLMTemplate.super.calcFluxes(self, dt)
	end

	return PLMTemplate
end

return PLMBehavior
