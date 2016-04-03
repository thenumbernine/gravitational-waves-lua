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

local methods = table{
	'calcPrimsByPressure',		-- 1-var (pressure) ... works!
--	'calcPrimsByZAndW',			-- 2-var (Z and W)
--	'calcPrimsByPrims',			-- 3-var (rho, vx, eInt) ... is stable & accurate, but 3-var means it's slow
}
--local len = methods:map(function(l) return #l end):sup()


function SRHD1DRoe:postIterate(dt)
	assert(not SRHD1DRoe.super.postIterate)
	local eqn = self.equation

	-- now that we have iterated once, our conservative variables are out of sync with our ws
--print()
	local maxBestDist = 0
	local maxIters = 0
	for i=1,self.gridsize do
		local prims = self.ws[i]
		local cons = self.qs[i]
--print('cell['..i..'] prims='..tolua(prims)..' cons='..tolua(cons))
	
		local failSolvers
		local best
		for _,method in ipairs(methods) do
			-- there's a few ways to go about this
			-- each has varying results
			local newPrims, failReason, iters = eqn[method](eqn, self, i, prims, cons)
			local newCons = newPrims and {eqn:consFromPrims(table.unpack(newPrims))}
			-- see how well the reconstruction worked...
			local dist = math.huge
			if newCons then
				dist = 0
				for j=1,3 do
					dist = dist + math.abs(cons[j] - newCons[j])
				end
			end
			local results = {prims=newPrims, cons=newCons, dist=dist, iters=iters}
--print('...'..('%-'..len..'s'):format(method)..' '..tolua(results))
			if not best or results.dist < best.dist then
				best = results
			end
		end
		assert(best)
		self.ws[i] = best.prims
		self.primitiveReconstructionErrors[i] = best.dist
--print('...best.dist='..best.dist)
		maxBestDist = math.max(maxBestDist, best.dist)
		maxIters = math.max(maxIters, best.iters)
	end
--print('maxBestDist='..maxBestDist..' maxIters='..maxIters)
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
