-- based on code at https://www.mathworks.com/matlabcentral/fileexchange/44639-weighted-essentially-non-oscillatory--weno--scheme
-- also based on Mara test code
local class = require 'ext.class'
local matrix = require 'matrix'
local weno5 = require 'weno5'

-- technically a Godunov solver and not a Roe solver
-- and that's where I should put the eigenbasis stuff
local SolverFV = require 'solverfv'
	
local WENO5FD = class(SolverFV)

WENO5FD.useXEdgeCentered = true

function WENO5FD:init(args)
	-- disable flux limiter
	args = table(args)
	local limiter = require 'limiter' 
	args.fluxLimiter = limiter.donorCell	

	self.weno5method = args.weno5method

	WENO5FD.super.init(self, args)
	self.name = 'WENO5FD'
end


WENO5FD.numGhost = 3


-- TODO 3.2.1 of 2013 Lou et al:
function WENO5FD:calcFluxes(dt)
	local numStates = self.numStates
	local numWaves = self.numWaves

	-- cell centered fluxes
	local Fs = self:newState()
	for i=1,self.gridsize do
		Fs[i] = matrix{self.equation:calcFluxForState(self.qs[i])}
	end

	-- cell centered eigenvalues
	local lambdas = self:newState()
	for i=1,self.gridsize do
		lambdas[i] = matrix()
		self.equation:calcEigenBasisFromCons(lambdas[i], nil, nil, nil, self.qs[i])
	end

	-- interface Roe averages
	local iRoes = self:newState()
	for i=2,self.gridsize do
		iRoes[i] = matrix{self.equation:calcRoeValues(self.qs[i-1], self.qs[i])}
	end

	-- interface eigenbasis
	local ievls = self:newState()
	local ievrs = self:newState()
	local ilambdas = self:newState()
	for i=2,self.gridsize do
		local iRoe = iRoes[i]
		ilambdas[i] = matrix.zeros{numStates}
		ievrs[i] = matrix.zeros{numStates, numStates}
		ievls[i] = matrix.zeros{numStates, numStates}
		self.equation:calcEigenBasis(ilambdas[i], ievrs[i], ievls[i], nil, iRoe:unpack())
	end

	-- now I guess we do polynomial interpolation on the neg and pos chars
	-- and then reconstruct with the right eigenvectors at the interfaces

	for i=3,self.gridsize-3 do
		-- the i'th interface is between cells i-1 and i
		-- for each offset -3..2 from i
		
		-- max lambda of all lambdas in this range of cells
		local ml = 0
		for j=-2,3 do
			for k=1,self.numWaves do
				ml = math.max(ml, math.abs(lambdas[i+j][k]))
			end
		end

		-- extrapolate cell fluxes in left and right directions
		--  according with the max wavespeed in this region
		local Fp = {}
		local Fm = {}
		for j=-2,3 do
			Fp[j] = matrix()
			Fm[j] = matrix()
			for k=1,numStates do
				Fp[j][k] = (Fs[i+j][k] + ml * self.qs[i+j][k]) * .5
				Fm[j][k] = (Fs[i+j][k] - ml * self.qs[i+j][k]) * .5
			end
		end

		-- transform it to char space
		local fp = {}
		local fm = {}
		for j=-2,2 do
			fp[j] = self.equation:eigenLeftTransform(self, ievls[i], Fp[j])
			fm[j+1] = self.equation:eigenLeftTransform(self, ievls[i], Fm[j+1])
		end

		local f = weno5(fp, 0, 'fd', 'r', self.weno5method, numWaves)
				+ weno5(fm, 1, 'fd', 'l', self.weno5method, numWaves)

		f = self.equation:eigenRightTransform(self, ievrs[i], f)

		for j=1,self.numStates do
			self.fluxes[i+1][j] = f[j]
		end
	end
end

return WENO5FD
