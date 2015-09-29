local class = require 'ext.class'
local table = require 'ext.table'
local Roe = require 'roe'
local ADM2DSpherical = require 'adm2dspherical'

local ADM2DSphericalRoe = class(Roe)

function ADM2DSphericalRoe:init(args)
	args = table(args)
	args.equation = ADM2DSpherical(args)
	ADM2DSphericalRoe.super.init(self, args)
end

function ADM2DSphericalRoe:iterate(...)
	ADM2DSphericalRoe.super.iterate(self, ...)
	-- enforce constraint of V_r = 2 * D_rhh / g_hh
	-- i.e. 1 = 2 * D_rhh / (g_hh * V_r)
	for i=1,self.gridsize do
		--[[ direct assign:
		self.qs[i][9] = 2 * self.qs[i][6] / self.qs[i][3]
		--]]
		-- [[ geometric renormalization
		local c = abs(2 * self.qs[i][6] / (self.qs[i][3] * self.qs[i][9]))^(1/3)
		self.qs[i][3] = self.qs[i][3] * c
		self.qs[i][6] = self.qs[i][6] / c
		self.qs[i][9] = self.qs[i][9] * c
		-- 2 * (self.qs[i][6] / cbrt(c)) / (self.qs[i][3] * cbrt(c) * self.qs[i][9] * cbrt(c))
		-- = (2 * self.qs[i][6] / (self.qs[i][3] * self.qs[i][9])) / c
		-- = 1
	end
end

return ADM2DSphericalRoe

