-- based on code at https://www.mathworks.com/matlabcentral/fileexchange/44639-weighted-essentially-non-oscillatory--weno--scheme
-- also based on Mara test code
local class = require 'ext.class'
local matrix = require 'matrix'
local weno5 = require 'weno5'

-- technically a Godunov solver and not a Roe solver
-- and that's where I should put the eigenbasis stuff
local SolverFV = require 'solverfv'
	
local WENO5FV = class(SolverFV)

-- why does this matter?
-- why does this help?
WENO5FV.useXEdgeCentered = true

function WENO5FV:init(args)
	-- disable flux limiter
	args = table(args)
	local limiter = require 'limiter' 
	args.fluxLimiter = limiter.donorCell	
	
	WENO5FV.super.init(self, args)
	self.name = 'WENO5FV'
end


WENO5FV.numGhost = 3


function WENO5FV:calcFluxes(dt)
	local numStates = self.numStates
	local numWaves = self.numWaves

	-- cell prims
	local Ws = self:newState()
	for i=1,self.gridsize do
		Ws[i] = matrix{self.equation:calcPrimFromCons(self.qs[i]:unpack())}
	end

	-- cell centered fluxes
	local Fs = self:newState()
	for i=1,self.gridsize do
		Fs[i] = matrix{self.equation:calcFluxForState(self.qs[i])}
	end

	for i=3,self.gridsize-3 do
		local Wl = weno5(Ws, i-2, 'fv', 'r', numWaves)
		local Wr = weno5(Ws, i-1, 'fv', 'l', numWaves)

		local Ul = matrix{self.equation:calcConsFromPrim(Wl:unpack())}
		local Ur = matrix{self.equation:calcConsFromPrim(Wr:unpack())}
	
		-- get the eigen system from the l and r
			-- left side
		local lambdaL = {}
		local evlL = {}
		local evrL = {}
			-- right side
		local lambdaR = {}
		local evlR = {}
		local evrR = {}
		for q=1,numWaves do
			evlL[q] = {}
			evrL[q] = {}
			evlR[q] = {}
			evrR[q] = {}
		end
		self.equation:calcEigenBasisFromCons(lambdaL, evlL, evrL, nil, Ul)
		self.equation:calcEigenBasisFromCons(lambdaR, evlR, evrR, nil, Ur)

		-- transform the state and flux into char space along the window
		local al = {}
		local ar = {}
		local afl = {}
		local afr = {}
		for j=-2,3 do
			al[j] = self.equation:eigenLeftTransform(self, evlL, self.qs[i+j])
			ar[j] = self.equation:eigenLeftTransform(self, evlR, self.qs[i+j])
			afl[j] = self.equation:eigenLeftTransform(self, evlL, Fs[i+j])
			afr[j] = self.equation:eigenLeftTransform(self, evlR, Fs[i+j])
		end

		-- do a minmod slope limiter <-> separate positive and negative waves
		--	with the flux char vars
	
		-- then reconstruct those
		-- then eeigen right mul
		-- and tada
	

		local fp = {}
		local fm = {}
		for j=-2,3 do
			fp[j] = {}
			fm[j] = {}
		end
		for q=1,numWaves do
			if lambdaL[q] > 0 and lambdaR[q] > 0 then
				-- No sign change, right-going waves only: set fm to zero and fp to f
				for j=-2,3 do
            		fp[j][q] = afl[j][q]
            		fm[j][q] = 0
          		end
			elseif lambdaL[q] < 0 and lambdaR[q] < 0 then
          		-- No sign change, left-going waves only: set fp to zero and fm to f
          		for j=-2,3 do
					fp[j][q] = 0
					fm[j][q] = afr[j][q]
          		end
			else
				-- There is a sign change in the speed of this characteristic field
				local abslambdal = math.abs(lambdaL[q])
				local abslambdar = math.abs(lambdaR[q])
				local a = math.max(abslambdal, abslambdar)
				for j=-2,3 do
					fp[j][q] = .5 * (afl[j][q] + a * al[j][q])
					fm[j][q] = .5 * (afr[j][q] - a * ar[j][q])
				end
        	end
      	end

		local fweno_p = weno5(fp, -2, 'fd', 'r', numWaves)
		local fweno_m = weno5(fm, -1, 'fd', 'l', numWaves)

		local Fweno_p = self.equation:eigenRightTransform(self, evrL, fweno_p)
		local Fweno_m = self.equation:eigenRightTransform(self, evrR, fweno_m)

		for j=1,self.numStates do
			self.fluxes[i+1][j] = Fweno_p[j] + Fweno_m[j]
		end
	end
end

return WENO5FV
