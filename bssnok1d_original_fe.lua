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
assert(math.isfinite(alpha))
assert(math.isfinite(phi))
assert(math.isfinite(gammaTilde_xx))
assert(math.isfinite(GammaTildeUx))
assert(math.isfinite(K))
assert(math.isfinite(ATilde_xx))
		local f = self.equation.calc.f(alpha)
assert(math.isfinite(f))
		
		local gammaTildeUxx = 1 / gammaTilde_xx
assert(math.isfinite(gammaTildeUxx ))
		local ATildeUxx = ATilde_xx * gammaTildeUxx * gammaTildeUxx
assert(math.isfinite(ATildeUxx ))
		local gamma_xx = gammaTilde_xx * phi^4
assert(math.isfinite(gamma_xx ))
		local gammaUxx = 1 / gamma_xx
assert(math.isfinite(gammaUxx ))

		local dx = self.ixs[i+1] - self.ixs[i]
		local _1_dx = 1 / dx
		local _1_dxSq = _1_dx * _1_dx
		
		local ip = math.min(i+1, self.gridsize)
		local im = math.max(i-1, 1)
		
		local alpha_p, phi_p, gammaTilde_xx_p, dGammaTildeUx_p, K_p, _ = table.unpack(dq_dt[ip])
		local alpha_m, phi_m, gammaTilde_xx_m, dGammaTildeUx_m, K_m, _ = table.unpack(dq_dt[im])	
		
		local dalpha_dx = .5 * _1_dx * (alpha_p - alpha_m)
assert(math.isfinite(dalpha_dx ))
		local dphi_dx = .5 * _1_dx * (phi_p - phi_m)
assert(math.isfinite(dphi_dx ))
		local dgammaTilde_xx_dx = .5 * _1_dx * (gammaTilde_xx_p - gammaTilde_xx_m)
assert(math.isfinite(dgammaTilde_xx_dx ))
		local dK_dx = .5 * _1_dx * (K_p - K_m)
assert(math.isfinite(dK_dx ))
		local dGammaTildeUx_dx = .5 * _1_dx * (dGammaTildeUx_p - dGammaTildeUx_m)
assert(math.isfinite(dGammaTildeUx_dx ))

		local d2alpha_dx2 = _1_dxSq * (alpha_p - 2 * alpha + alpha_m)
assert(math.isfinite(d2alpha_dx2 ))
		local d2phi_dx2 = _1_dxSq * (phi_p - 2 * phi + phi_m)
assert(math.isfinite(d2phi_dx2 ))
		local d2gammaTilde_xx_dx2 = _1_dxSq * (gammaTilde_xx_p - 2 * gammaTilde_xx + gammaTilde_xx_m)
assert(math.isfinite(d2gammaTilde_xx_dx2 ))
		
		-- for 1D:
		-- GammaTilde^x_xx gammaTilde^xx = GammaTilde^x
		-- so GammaTilde^x_xx = GammaTilde^x gammaTilde_xx
		local GammaTildeUx_xx = GammaTildeUx * gammaTilde_xx
assert(math.isfinite(GammaTildeUx_xx))

		-- GammaTilde_ijk = gammaTilde_il GammaTilde^l_jk
		local GammaTilde_xxx = gammaTilde_xx * GammaTildeUx_xx
assert(math.isfinite(GammaTilde_xxx))

		-- Gamma^k_ij = GammaTilde^k_ij + 2 (delta^k_i phi_,j + delta^k_j phi_,i - gamma_ij gamma^kl phi_,l)
		-- Gamma^x_xx = GammaTilde^x_xx + 2 phi_,x
		local GammaUx_xx = GammaTildeUx_xx + 2 * dphi_dx
assert(math.isfinite(GammaUx_xx))

		-- Gamma^i = Gamma^i_jk gamma^jk = gamma^il Gamma_ljk gamma^jk
		--	for 1D: = 1/2 gamma_xx,x / gamma_xx^2
		local GammaUx = GammaUx_xx * gammaUxx
assert(math.isfinite(GammaUx))
		-- or you could just use
		-- GammaTilde^i = exp(4 phi) Gamma^i + 2 gammaTilde^ij phi_,j

		local RTilde_xx = -1/2 * gammaTildeUxx * d2gammaTilde_xx_dx2
			+ gammaTilde_xx * dGammaTildeUx_dx + GammaTildeUx * gammaTilde_xx * GammaTildeUx_xx
			+ 3 * gammaTildeUxx * GammaTildeUx_xx * GammaTilde_xxx 
assert(math.isfinite(RTilde_xx))
		local Rphi_xx = -4 * d2phi_dx2 + 4 * GammaTildeUx_xx * dphi_dx
assert(math.isfinite(Rphi_xx))
		local R_xx = RTilde_xx + Rphi_xx
assert(math.isfinite(R_xx))
		local jTildeUx = 0
		local rho = 0
		local S = 0
		local S_xx = 0

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
		-- distribute the trace-free
		-- ATilde_ij,t = exp(-4 phi) (-(D_i D_j alpha)^TF + alpha R_ij^TF - 8 pi alpha S_ij^TF)
		-- 		+ alpha (K ATilde_ij - 2 ATilde_ik gammaTilde^kl ATilde_lj)
		
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
