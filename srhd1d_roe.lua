local class = require 'ext.class'
local SRHD1D = require 'srhd1d'
local Roe = require 'roe'

local SRHD1DRoe = class(Roe)
SRHD1DRoe.equation = SRHD1D()

function SRHD1DRoe:init(...)
	SRHD1DRoe.super.init(self, ...)

	-- used by the equation object, so the two are tied together pretty closely
	self.primitives = {}
end

function SRHD1DRoe:reset(...)
	-- make sure the entries are reset
	self.primitives = {}
	SRHD1DRoe.super.reset(self, ...)
end

function SRHD1DRoe:postIterate(dt)
	assert(not SRHD1DRoe.super.postIterate)

	-- now that we have iterated once, our conservative variables are out of sync with our primitives
	for i=1,self.gridsize do
		self.equation:calcPrims(self, i, self.primitives[i], self.qs[i])
	end
end

-- hack the boundary function to apply to primitives as well
-- TODO implement boundaries in Equations so this doesn't have to be overridden
function SRHD1DRoe:applyBoundary(...)
	SRHD1DRoe.super.applyBoundary(self, ...)

	-- freeflow on prims...
	local prims = self.primitives
	for _,k in ipairs{'rho','vx','h','eInt','P'} do
		prims[1][k] = assert(prims[3][k], "failed to find primitives[3]."..k)
		prims[2][k] = assert(prims[3][k], "failed to find primitives[3]."..k)
		prims[self.gridsize][k] = assert(prims[self.gridsize-2][k], "failed to find primitives["..(self.gridsize-2).."]."..k)
		prims[self.gridsize-1][k] = assert(prims[self.gridsize-2][k], "failed to find primitives["..(self.gridsize-2).."]."..k)
	end
end

return SRHD1DRoe
