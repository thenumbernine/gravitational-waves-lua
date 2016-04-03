local class = require 'ext.class'
local SRHD1D = require 'srhd1d'
local Roe = require 'roe'

local SRHD1DRoe = class(Roe)
SRHD1DRoe.equation = SRHD1D()

function SRHD1DRoe:init(...)
	SRHD1DRoe.super.init(self, ...)

	-- used by the equation object, so the two are tied together pretty closely
	-- ws[i][1] = rho, ws[i][2] = vel, ws[i][3] = eInt
	self.ws = {}
	self.primitiveReconstructionErrors = {}
end

function SRHD1DRoe:reset(...)
	-- make sure the entries are reset
	self.ws = {}
	self.primitiveReconstructionErrors = {}

	SRHD1DRoe.super.reset(self, ...)
end

function SRHD1DRoe:postIterate(dt)
	assert(not SRHD1DRoe.super.postIterate)
	local eqn = self.equation

	local methods = table{
		'calcPrimsByPrims',
		'calcPrimsByPressure',
		'calcPrimsByZAndW',
	}
	local len = methods:map(function(l) return #l end):sup()

	-- now that we have iterated once, our conservative variables are out of sync with our ws
--print()
	local maxBestDist = 0
	for i=1,self.gridsize do
		local prims = self.ws[i]
		local cons = self.qs[i]
--print('cell['..i..'] prims='..tolua(prims)..' cons='..tolua(cons))
		
		local best = methods:map(function(method)
			-- there's a few ways to go about this
			-- each has varying results
			local newPrims = eqn[method](eqn, self, i, prims, cons)
			local newCons = newPrims and {eqn:consFromPrims(table.unpack(newPrims))}
			-- see how well the reconstruction worked...
			local dist = math.huge
			if newCons then
				dist = 0
				for j=1,3 do
					dist = dist + math.abs(cons[j] - newCons[j])
				end
			end
--print('...'..('%-'..len..'s'):format(method)..' prims='..tolua(newPrims)..' cons='..tolua(newCons)..' dist='..dist)
			return {prims=newPrims, cons=newCons, dist=dist}
		end):sort(function(a,b)
			return a.dist < b.dist
		end)[1]
		assert(best)
		self.ws[i] = best.prims
		self.primitiveReconstructionErrors[i] = best.dist
--print('...best.dist='..best.dist)
		maxBestDist = math.max(maxBestDist, best.dist)
	end
--print('maxBestDist='..maxBestDist)
end

-- hack the boundary function to apply to primitives as well
-- TODO implement boundaries in Equations so this doesn't have to be overridden
function SRHD1DRoe:applyBoundary(...)
	SRHD1DRoe.super.applyBoundary(self, ...)

	-- freeflow on prims...
	local prims = self.ws
	for k=1,3 do
		prims[1][k] = assert(prims[3][k], 'failed to find ws[3].'..k)
		prims[2][k] = assert(prims[3][k], 'failed to find ws[3].'..k)
		prims[self.gridsize][k] = assert(prims[self.gridsize-2][k], 'failed to find ws['..(self.gridsize-2)..'].'..k)
		prims[self.gridsize-1][k] = assert(prims[self.gridsize-2][k], 'failed to find ws['..(self.gridsize-2)..'].'..k)
	end
end

return SRHD1DRoe
