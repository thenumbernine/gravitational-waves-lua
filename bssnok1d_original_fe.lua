local Solver = require 'solver'
local BSSNOK1DOriginal = require 'bssnok1d_original'

local BSSNOK1DOriginalForwardEuler = class(Solver)
BSSNOK1DOriginalForwardEuler.fixed_dt = 1/100

function BSSNOK1DOriginalForwardEuler:init(args)
	args.equation = BSSNOK1DOriginal(args)
	BSSNOK1DOriginalForwardEuler.super.init(self, args)
	self.fluxes = {}
end

function BSSNOK1DOriginalForwardEuler:reset()
	BSSNOK1DOriginalForwardEuler.super.reset(self)
end

function BSSNOK1DOriginalForwardEuler:calcDT()
	return self.fixed_dt
end

function BSSNOK1DOriginalForwardEuler:calcDerivFromFluxes(dt)
	local dq_dt = self:newState()
	for i=1,self.gridsize do
		local alpha, phi, gammaTilde_xx, GammaTildeUx, K, ATilde_xx = table.unpack(self.qs[i])
		local f = self.equation.calc.f(alpha)
		
		local gammaTildeUxx = 1 / gammaTilde_xx
		local ATildeUxx = ATilde_xx * gammaTildeUxx * gammaTildeUxx
		local gamma_xx = gammaTilde_xx * phi^4
		local gammaUxx = 1 / gamma_xx

		local dx = self.ixs[i+1] - self.ixs[i]
		local ip = math.min(i+1, self.gridsize)
		local im = math.max(i-1, 1)
		
		local alpha_p, phi_p, gammaTilde_xx_p, _, K_p, _ = table.unpack(dq_dt[ip])
		local alpha_m, phi_m, gammaTilde_xx_m, _, K_m, _ = table.unpack(dq_dt[im])
		
		local gamma_xx_p = gammaTilde_xx_p * phi_p^4
		local gamma_xx_m = gammaTilde_xx_m * phi_m^4
		
		local _1_dx = 1 / dx
		local _1_dxSq = _1_dx * _1_dx
		
		local dalpha_dx = .5 * _1_dx * (alpha_p - alpha_m)
		local dphi_dx = .5 * _1_dx * (phi_p - phi_m)
		local dgammaTilde_xx_dx = .5 * _1_dx * (gammaTilde_xx_p - gammaTilde_xx_m)
		local dK_dx = .5 * _1_dx * (K_p - K_m)

		local d2alpha_dx2 = _1_dxSq * (alpha_p - 2 * alpha + alpha_m)
		local d2phi_dx2 = _1_dxSq * (phi_p - 2 * phi + phi_m)
		local d2gammaTilde_xx_dx2 = _1_dxSq * (gammaTilde_xx_p - 2 * gammaTilde_xx + gammaTilde_xx_m)

		local dgamma_xx_dx = .5 * _1_dx * (gamma_xx_p - gamma_xx_m)
		
		-- Gamma^i_jk in 1D:
		-- Gamma^x_xx = gamma^xx Gamma_xxx = 1/2 gamma^xx gamma_xx,x
		local GammaUx_xx = 1/2 * dgamma_xx_dx * gammaUxx
		-- Gamma^i = Gamma^i_jk gamma^jk = gamma^il Gamma_ljk gamma^jk
		--	for 1D: = 1/2 gamma_xx,x / gamma_xx^2
		local GammaUx = GammaUx_xx * gammaUxx
		-- or you could just use
		-- GammaTilde^i = exp(4 phi) Gamma^i + 2 gammaTilde^ij phi_,j

		-- GammaTilde^k_ij = Gamma^k_ij - 2 (delta^k_i phi_,j + delta^k_j phi_,i - gamma_ij gamma^kl phi_,l)
		-- for 1D:
		-- GammaTilde^x_xx gammaTilde^xx = GammaTilde^x
		-- so GammaTilde^x_xx = GammaTilde^x gammaTilde_xx
		local GammaTildeUx_xx = GammaTildeUx * gammaTilde_xx
		
		local RTilde_xx = -1/2 * gammaTildeUxx * d2gammaTilde_xx_dx2
			+ gammaTilde_xx * dGammaTildeUx_dx + GammaTildeUx * gammaTilde_xx * GammaTildeUx_xx
			+ 3 * gammaTildeUxx * GammaTildeUx_xx * GammaTilde_xxx 
		local Rphi_xx = -4 * d2phi_dx2 + 4 * GammaTildeUx_xx * dphi_dx
		local R_xx = RTilde_xx + Rphi_xx

		local jTildeUx = 0
		local rho = 0
		local S = 0

		-- alpha,t = -alpha^2 f K
		dq_dt[i][1] = -alpha * alpha * f * K
		-- phi,t = -1/6 alpha K
		dq_dt[i][2] = -1/6 * alpha * K
		-- gammaTilde_ij,t = -2 alpha ATilde_ij
		dq_dt[i][3] = -2 * alpha * ATilde_xx
		-- GammaTilde^i_,t = -2 ATilde^ij alpha_,j
		--	+ 2 alpha (GammaTilde^i_jk ATilde^jk
		-- 		+ 6 ATilde^ij phi_,j
		--		- 2/3 gammaTilde^ij K_,j
		--		- 8 pi jTilde^i)
		dq_dt[i][4] = -2 * ATildeUxx * dalpha_dx
			+ 2 * alpha * (GammaTildeUx_xx * ATildeUxx
				+ 6 * ATildeUxx * dphi_dx
				- 2/3 * gammaTildeUxx * dK_dx
				- 8 * math.pi * jTildeUx)
		-- K_,t = -D_i D^i alpha + alpha (ATilde_ij ATilde^ij + K^2 / 3) + 4 pi alpha (rho + S)
		-- K_,t = -gamma^ij (alpha_,ij - Gamma^k_ij alpha_,k) + alpha (ATilde_ij ATilde^ij + K^2 / 3) + 4 pi alpha (rho + S)
		-- K_,t = -gamma^ij alpha_,ij + Gamma^k alpha_,k + alpha (ATilde_ij ATilde^ij + K^2 / 3) + 4 pi alpha (rho + S)
		dq_dt[i][5] = -gammaUxx * d2alpha_dx2 + GammaUx * dalpha_dx
			+ alpha * (ATilde_xx * ATildeUxx + K * K / 3)
			+ 4 * math.pi * alpha * (rho + S)
		-- ATilde_ij,t = exp(-4 phi) (-D_i D_j alpha + alpha R_ij + 4 pi alpha (gamma_ij (S - rho) - 2 S_ij))^TF
		-- 		+ alpha (K ATilde_ij - 2 ATilde_ik gamma^kl ATilde_lj)
		-- 1D trace-free of a 3D tensor: T_ij - 1/3 gamma_ij T^k_k
		-- = T_xx - 1/3 gamma_xx T_xx / gamma_xx = T_xx - 1/3 T_xx = 2/3 T_xx
		-- ATilde_ij,t = 2/3 exp(-4 phi) (-alpha,ij + Gamma^k_ji alpha,k + alpha R_ij + 4 pi alpha (gamma_ij (S - rho) - 2 S_ij))
		-- 		+ alpha (K ATilde_ij - 2 ATilde_ik gammaTilde^kl ATilde_lj)
		dq_dt[i][6] = 2/3 * math.exp(-4 * phi) * (-d2alpha_dx2 + GammaUx_xx * dalpha_dx + alpha * R_xx + 4 * math.pi * alpha * (gamma_xx * (S - rho) - 2 * S_xx))
			+ alpha * (K * ATilde_xx - 2 * ATilde_xx * ATilde_xx * gammaTildeUxx)
	end
	return dq_dt
end

return BSSNOK1DOriginalForwardEuler
