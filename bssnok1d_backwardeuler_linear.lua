--[[

diff_t alpha = -alpha^2 f K
diff_t gammaTilde_ij = -2 alpha ATilde_ij
diff_t phi = -1/6 alpha K

diff_t ATilde_ij =
	exp(-4 phi) (
		-diff_i diff_j alpha
		+ 1/3 gammaTilde_ij gammaTilde^kl diff_k diff_l alpha
		
		+ GammaTilde^k_ji diff_k alpha
		- 1/3 gammaTilde_ij GammaTilde^k diff_alpha
		
		+ 2 delta^k_i diff_j phi diff_k alpha
		- 2/3 gamma_ij gamma^mn delta^k_m diff_n phi diff_k alpha
		
		+ 2 delta^k_j diff_i phi diff_k alpha
		- 2/3 gamma_ij gamma^mn delta^k_n diff_m phi diff_k alpha
		
		+ alpha R_ij
		- 1/3 alpha gamma_ij R

		- 8 pi alpha S_ij
		+ 8/3 pi alpha gamma_ij gamma^kl S_kl
	)
	+ alpha K ATilde_ij
	- 2 alpha ATilde_ik gammaTilde^kl ATilde_lj

diff_t K =
	exp(-4 phi) (
		-gammaTilde^ij diff_i diff_j alpha
		+ GammaTilde^k diff_k alpha
		- 2 gammaTilde^ij diff_j phi diff_i alpha
	)
	+ alpha ATilde_ij gammaTilde^im gammaTilde^jn ATilde_mn
	+ 1/3 alpha K^2
	+ 4 pi alpha (rho + S)

diff_t GammaTilde^i =
	-2 gammaTilde^im ATilde^mn gammaTilde^nj diff_j alpha
	+ 2 alpha GammaTilde^i_jk gammaTilde^jm ATilde^mn gammaTilde^nk
	+ 12 alpha gammaTilde^im ATilde_mn gammaTilde^nj diff_j phi
	- 4/3 alpha gammaTilde^ij diff_j K
	- 16 pi alpha jTilde^i

--]]

local BSSNOK1DBackwardEulerLinear = require 'ext.class'(require 'solver')

function BSSNOK1DBackwardEulerLinear:init(args)
	args = require 'ext.table'(args)
	args.equation = require 'bssnok1d'(args)
	BSSNOK1DBackwardEulerLinear.super.init(self, args)
end

function BSSNOK1DBackwardEulerLinear:iterate()
end

return BSSNOK1DBackwardEulerLinear

