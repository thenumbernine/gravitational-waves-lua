local class = require 'ext.class'
local table = require 'ext.table'

local Roe = require 'roe'

local ADM2DSphericalRoe = class(Roe)

function ADM2DSphericalRoe:init(args)
	args = table(args)
	args.equation = require 'adm2dspherical'(args)
	ADM2DSphericalRoe.super.init(self, args)
end

function ADM2DSphericalRoe:addSourceToDerivCell(dq_dts, i)
	local alpha, g_rr, g_hh, A_r, D_rrr, D_rhh, K_rr, K_hh, V_r = unpack(self.qs[i])
	local f = self.equation.calc_f(alpha)
	local tr_K = K_rr / g_rr + K_hh / g_hh
	local S_rr = K_rr * (2 * K_hh / g_hh - K_rr / g_rr) + A_r * (D_rrr / g_rr - 2 * D_rhh / g_hh)
				+ 2 * D_rhh / g_hh * (D_rrr / g_rr - D_rhh / g_hh) + 2 * A_r * V_r
	local S_hh = K_rr * K_hh / g_rr - D_rrr * D_rhh / g_rr^2 + 1
	local P_r = -2 / g_hh * (A_r * K_hh - D_rhh * (K_hh / g_hh - K_rr / g_rr))
	dq_dts[i][1] = dq_dts[i][1] - alpha * alpha * f * tr_K 
	dq_dts[i][2] = dq_dts[i][2] - 2 * alpha * K_rr
	dq_dts[i][3] = dq_dts[i][3] - 2 * alpha * K_hh
	dq_dts[i][7] = dq_dts[i][7] + alpha * S_rr
	dq_dts[i][8] = dq_dts[i][8] + alpha * S_hh
	dq_dts[i][9] = dq_dts[i][9] + alpha * P_r
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

