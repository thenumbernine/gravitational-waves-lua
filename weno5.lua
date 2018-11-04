-- based on code at https://www.mathworks.com/matlabcentral/fileexchange/44639-weighted-essentially-non-oscillatory--weno--scheme
-- also based on Mara test code
local class = require 'ext.class'
local matrix = require 'matrix'

-- technically a Godunov solver and not a Roe solver
-- and that's where I should put the eigenbasis stuff
local SolverFV = require 'solverfv'
	
local WENO5 = class(SolverFV)

function WENO5:init(args)
	-- disable flux limiter
	args = table(args)
	local limiter = require 'limiter' 
	args.fluxLimiter = limiter.donorCell	
	
	WENO5.super.init(self, args)
	self.name = 'WENO5'
end


local CeesC2L = { {11./6., -7./6.,  1./3. },
				{ 1./3.,  5./6., -1./6. },
				{-1./6.,  5./6.,  1./3. } }
local CeesC2R = { { 1./3.,  5./6., -1./6. },
				{-1./6.,  5./6.,  1./3. },
				{ 1./3., -7./6., 11./6. } }
local DeesC2L = { 0.1, 0.6, 0.3 }
local DeesC2R = { 0.3, 0.6, 0.1 }


-- TODO 3.2.1 of 2013 Lou et al:
function WENO5:calcFluxes(dt)
	local numStates = self.numStates
	local numWaves = self.numWaves

--print('WENO5:calcFluxes')

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
--print('iRoes['..i..'] = '..iRoes[i])
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
--print('ievls['..i..'] = '..tolua(ievls[i]))
--print('ievrs['..i..'] = '..tolua(ievrs[i]))
--print('ilambdas['..i..'] = '..tolua(ilambdas[i]))
	end

	-- now I guess we do polynomial interpolation on the neg and pos chars
	-- and then reconstruct with the right eigenvectors at the interfaces

	for i=5,self.gridsize-2 do
		-- the i'th interface is between cells i-1 and i
		-- for each offset -3..2 from i
		
		-- max lambda of all lambdas in this range of cells
		local ml = 0
		for j=-3,2 do
			for k=1,self.numWaves do
				ml = math.max(ml, math.abs(lambdas[i+j][k]))
			end
		end

		-- extrapolate cell fluxes in left and right directions
		--  according with the max wavespeed in this region
		local Fp = {}
		local Fm = {}
		for j=-3,2 do
			-- [[ plus and minus extrapolation by largest eigenvector
			Fp[j] = matrix()
			Fm[j] = matrix()
			for k=1,numStates do
				Fp[j][k] = (Fs[i+j][k] + ml * self.qs[i+j][k]) * .5
				Fm[j][k] = (Fs[i+j][k] - ml * self.qs[i+j][k]) * .5
			end
			--]]
			--[[ how about positive and negative eigenvalues?
			-- F = A * U = evr * lambda * evl * U
			local char = self.equation:eigenLeftTransform(self, ievls[i], self.qs[i+j])
			local charP = matrix()
			local charM = matrix()
			for k=1,self.numWaves do
				charP[k] = char[k] * math.max(0, ilambdas[i][k])
				charM[k] = char[k] * math.min(0, ilambdas[i][k])
			end
			Fp[j] = self.equation:eigenRightTransform(self, ievls[i], charP)
			Fm[j] = self.equation:eigenRightTransform(self, ievls[i], charM)
			--]]
		end

		-- transform it to char space
		local fp = {}
		local fm = {}
		for j=-3,2 do
			fp[j] = self.equation:eigenLeftTransform(self, ievls[i], Fp[j])
			fm[j] = self.equation:eigenLeftTransform(self, ievls[i], Fm[j])
		end

		local function weno5(v, ofs, c, d)
			local eps = 1e-6
			local result = matrix()
			for k=1,numWaves do
				local B = {
					(13./12.)*(  v[ofs+2][k] - 2*v[ofs+3][k] +   v[ofs+4][k])^2 +
					 ( 1./ 4.)*(3*v[ofs+2][k] - 4*v[ofs+3][k] +   v[ofs+4][k])^2,

					 (13./12.)*(  v[ofs+1][k] - 2*v[ofs+2][k] +   v[ofs+3][k])^2 +
					 ( 1./ 4.)*(  v[ofs+1][k] - 0*v[ofs+2][k] -   v[ofs+3][k])^2,
					 
					 (13./12.)*(  v[ofs+0][k] - 2*v[ofs+1][k] +   v[ofs+2][k])^2 +
					 ( 1./ 4.)*(  v[ofs+0][k] - 4*v[ofs+1][k] + 3*v[ofs+2][k])^2}

				local vs = {c[1+0][1+0]*v[ofs+2][k] + c[1+0][1+1]*v[ofs+3][k] + c[1+0][1+2]*v[ofs+4][k],
						  c[1+1][1+0]*v[ofs+1][k] + c[1+1][1+1]*v[ofs+2][k] + c[1+1][1+2]*v[ofs+3][k],
						  c[1+2][1+0]*v[ofs+0][k] + c[1+2][1+1]*v[ofs+1][k] + c[1+2][1+2]*v[ofs+2][k]}

				local w = {d[1+0] / (eps + B[1+0])^2,
						 d[1+1] / (eps + B[1+1])^2,
						 d[1+2] / (eps + B[1+2])^2}
				
				local wtot = w[1+0] + w[1+1] + w[1+2]
				result[k] = (w[1+0]*vs[1+0] + w[1+1]*vs[1+1] + w[1+2]*vs[1+2])/wtot
			end
			return result
		end

		local f = (weno5(fp, -3, CeesC2R, DeesC2R)
					+ weno5(fm, -2, CeesC2L, DeesC2L)) * .5

		f = self.equation:eigenRightTransform(self, ievrs[i], f)

		for j=1,self.numStates do
			self.fluxes[i][j] = f[j]
		end
	end
	-- [[
	local size = 5
	for j=1,self.numStates do
		for k=1,size-1 do
			self.fluxes[k][j] = self.fluxes[size][j]
			self.fluxes[self.gridsize+2-k][j] = self.fluxes[self.gridsize+2-size][j]
		end
	end
	--]]
end

-- I always hard-coded numGhost=2 so far, but WENO5 needs 3 ...
function WENO5:applyBoundary()
	local size = 5
	for j=1,self.numStates do
		for k=1,size-1 do
			self.qs[k][j] = self.qs[size][j]
			self.qs[self.gridsize+1-k][j] = self.qs[self.gridsize+1-size][j]
		end
	end
end

return WENO5
