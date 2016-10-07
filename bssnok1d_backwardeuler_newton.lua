--[[

p.85 non-hyperbolic description:
partial_t gammaTilde_ij = -2 alpha ATilde_ij
partial_t phi = -1/6 alpha K
partial_t ATilde_ij = exp(-4 phi) ( -D_i D_j alpha + alpha R_ij + 4 pi alpha (g_ij (S - rho) - 2 S_ij))^TF + alpha (K ATilde_ij - 2 ATilde_ik ATilde^k_j)
partial_t K = -D_i D^i alpha + alpha (ATilde_ij ATilde^ij + 1/3 K^2) + 4 pi alpha (rho + S)
partial_t GammaTilde^i = -2 ATilde^ij partial_j alpha + 2 alpha (GammaTilde^i_jk ATilde^jk + 6 ATilde^ij partial_j phi - 2/3 gammaTilde^ij partial_j K - 8 pi jTilde^i)
partial_t alpha = -alpha^2 f K			<- Bona-Masso lapse

rewrite to convert covariant derivatives to partial derivatives...

partial_t ATilde_ij = exp(-4 phi) ( -D_i D_j alpha + alpha R_ij + 4 pi alpha (g_ij (S - rho) - 2 S_ij))^TF + alpha (K ATilde_ij - 2 ATilde_ik ATilde^k_j)
partial_t ATilde_ij = exp(-4 phi) ( -D_i (partial_j alpha) + alpha R_ij + 4 pi alpha (g_ij (S - rho) - 2 S_ij))^TF + alpha (K ATilde_ij - 2 ATilde_ik ATilde^k_j)
partial_t ATilde_ij = exp(-4 phi) ( -partial_i partial_j alpha + Gamma^k_ji partial_k alpha + alpha R_ij + 4 pi alpha (g_ij (S - rho) - 2 S_ij))^TF + alpha (K ATilde_ij - 2 ATilde_ik ATilde^k_j)
	... using Alcubierre 2.8.14:
			GammaTilde^k_ij = Gamma^k_ij - 2 (delta^k_i partial_j phi + delta^k_j partial_i phi - gamma_ij gamma^kl partial_l phi)
			Gamma^k_ij = GammaTilde^k_ij + 2 (delta^k_i partial_j phi + delta^k_j partial_i phi - gamma_ij gamma^kl partial_l phi)
			Gamma^k_ij = GammaTilde^k_ij + 2 (delta^k_i partial_j phi + delta^k_j partial_i phi - gammaTilde_ij gammaTilde^kl partial_l phi)
partial_t ATilde_ij =
	exp(-4 phi) (
		-partial_i partial_j alpha
		+ (GammaTilde^k_ji + 2 (delta^k_i partial_j phi + delta^k_j partial_i phi - gammaTilde_ij gammaTilde^kl partial_l phi)) partial_k alpha
		+ alpha R_ij + 4 pi alpha (g_ij (S - rho) - 2 S_ij)
	)^TF
	+ alpha K ATilde_ij
	- 2 alpha ATilde_ik gammaTilde^kl ATilde_lj
	... expand trace-free ...
partial_t ATilde_ij =
	exp(-4 phi) (
		-partial_i partial_j alpha
		+ 1/3 gammaTilde_ij gammaTilde^kl partial_k partial_l alpha
		
		+ GammaTilde^k_ji partial_k alpha
		- 1/3 gammaTilde_ij GammaTilde^k partial_alpha
		
		+ 2 delta^k_i partial_j phi partial_k alpha
		- 2/3 gamma_ij gamma^mn delta^k_m partial_n phi partial_k alpha
		
		+ 2 delta^k_j partial_i phi partial_k alpha
		- 2/3 gamma_ij gamma^mn delta^k_n partial_m phi partial_k alpha
		
		+ alpha R_ij^TF
		
		- 8 pi alpha S_ij
		+ 8/3 pi alpha gamma_ij gamma^kl S_kl
	)
	+ alpha K ATilde_ij
	- 2 alpha ATilde_ik gammaTilde^kl ATilde_lj

partial_t K = -D_i D^i alpha + alpha (ATilde_ij ATilde^ij + 1/3 K^2) + 4 pi alpha (rho + S)
partial_t K = -D_i (gamma^ij D_j alpha) + alpha (ATilde_ij ATilde^ij + 1/3 K^2) + 4 pi alpha (rho + S)
partial_t K = -((D_i gamma^ij) partial_j alpha + gamma^ij D_i partial_j alpha) + alpha (ATilde_ij ATilde^ij + 1/3 K^2) + 4 pi alpha (rho + S)
	... using D_i gamma^ij = 0 ...
partial_t K = -gamma^ij (partial_i partial_j alpha - Gamma^k_ji partial_k alpha) + alpha (ATilde_ij ATilde^ij + 1/3 K^2) + 4 pi alpha (rho + S)
partial_t K = -gamma^ij partial_i partial_j alpha + gamma^ij Gamma^k_ji partial_k alpha + alpha (ATilde_ij ATilde^ij + 1/3 K^2) + 4 pi alpha (rho + S)
partial_t K = -gamma^ij partial_i partial_j alpha + Gamma^k partial_k alpha + alpha (ATilde_ij ATilde^ij + 1/3 K^2) + 4 pi alpha (rho + S)
	... using Alcubierre 2.8.15:
			GammaTilde^i = e^(4 phi) Gamma^i + 2 gammaTilde^ij partial_j phi
			e^(-4 phi) GammaTilde^i = Gamma^i + 2 e^(-4 phi) gammaTilde^ij partial_j phi
			Gamma^i = e^(-4 phi) GammaTilde^i - 2 e^(-4 phi) gammaTilde^ij partial_j phi
partial_t K = -gamma^ij partial_i partial_j alpha + (e^(-4 phi) GammaTilde^k - 2 e^(-4 phi) gammaTilde^kl partial_l phi) partial_k alpha + alpha (ATilde_ij ATilde^ij + 1/3 K^2) + 4 pi alpha (rho + S)
	... using Alcubierre 2.8.1:
			gammaTilde_ij = exp(-4 phi) gamma_ij
			gammaTilde^ij = exp(4 phi) gamma^ij		<- inverse of the above equation
			gamma^ij = exp(-4 phi) gammaTilde^ij	<- solve for gamma^ij
partial_t K =
	exp(-4 phi) (
		-gammaTilde^ij partial_i partial_j alpha
		+ GammaTilde^k partial_k alpha
		- 2 gammaTilde^ij partial_j phi partial_i alpha
	)
	+ alpha ATilde_ij gammaTilde^im gammaTilde^jn ATilde_mn
	+ 1/3 alpha K^2
	+ 4 pi alpha (rho + S)

partial_t GammaTilde^i =
	-2 gammaTilde^im ATilde^mn gammaTilde^nj partial_j alpha
	+ 2 alpha GammaTilde^i_jk gammaTilde^jm ATilde^mn gammaTilde^nk
	+ 12 alpha gammaTilde^im ATilde_mn gammaTilde^nj partial_j phi
	- 4/3 alpha gammaTilde^ij partial_j K
	- 16 pi alpha jTilde^i


state variables:
alpha, gammaTilde, phi, ATilde, K, GammaTilde^i

partial_t alpha = -alpha^2 f K
	first derivatives:
	
		d/dalpha partial_t alpha = d/dalpha( -alpha^2 f K)
		d/dalpha partial_t alpha = -alpha K (2 f + alpha d/dalpha(f))

		d/dK partial_t alpha = -alpha^2 f

	second derivatives:
		
		d/dalpha d/dalpha partial_t alpha = -K (2 f + alpha d/dalpha(f)) - alpha K (2 d/dalpha(f) + df/dalpha + alpha d^2/dalpha^2(f))
		d/dalpha d/dalpha partial_t alpha = -K (2 f + alpha (4 d/dalpha(f) + alpha d^2/dalpha^2(f)))
		
		d/dalpha d/dK partial_t alpha = -2 alpha f

partial_t gammaTilde_ij = -2 alpha ATilde_ij
	first derivatives:

		d/dalpha partial_t gammaTilde_ij = -2 ATilde_ij
		
		d/dATilde_mn partial_t gammaTilde_ij = -2 alpha delta^m_i delta^n_j

	second derivatives:
		
		d/dATilde_mn d/dalpha partial_t gammaTilde_ij = -2 delta^m_i delta^n_j

partial_t phi = -1/6 alpha K
	first derivatives:

		d/dalpha partial_t phi = -1/6 K
		
		d/dK partial_t phi = -1/6 alpha
	
	second derivatives:
	
		d/dK d/dalpha partial_t phi = -1/6

partial_t ATilde_ij
	first derivatives:

partial_t K
	first derivatives:

partial_t GammaTilde^i
	first derivatives:


Backward Euler:

q_i(t+dt) = q_i(t) + dt * dq_i/dt(t+dt)
q_i(t+dt) - dt * dq_i/dt(t+dt) = q_i(t)
	f_i(q) = dq_i/dt
q_i(t+dt) - dt * f_i(q(t+dt)) = q_i(t)
g_i = q_i(t+dt) - dt * f_i(q(t+dt)) - q_i(t)
minimizing q_i(t+dt)

g_i = q_i(t+dt) - dt * f_i(q_1(t+dt), ..., q_n(t+dt)) - q_i(t)

e = ||g|| = 1/2 sum_i (g_i * g_i)

de/dq_i(t+dt) = sum_j (g_j * d/dq_i(t+dt) g_j)
de/dq_i(t+dt) = sum_j (g_j * d/dq_i(t+dt) (q_j(t+dt) - dt * f_j - q_j(t)))
de/dq_i(t+dt) = sum_j (g_j * (delta_ji - dt * d/dq_i(t+dt) f_j ))
de/dq_i(t+dt) = g_i - dt * sum_j (d/dq_i(t+dt) f_j * g_j)

d/dq_i(t+dt) d/dq_j(t+dt) e = d/dq_j(t+dt) (g_i - dt * sum_k (g_k * d/dq_i(t+dt) f_k ))
d/dq_i(t+dt) d/dq_j(t+dt) e = d/dq_j(t+dt) g_i - dt * sum_k d/dq_j(t+dt) (g_k * d/dq_i(t+dt) f_k )
d/dq_i(t+dt) d/dq_j(t+dt) e = d/dq_j(t+dt) g_i - dt * sum_k (d/dq_j(t+dt) g_k * d/dq_i(t+dt) f_k + g_k * d/dq_i(t+dt) d/dq_j(t+dt) f_k )
d/dq_i(t+dt) d/dq_j(t+dt) e = d/dq_j(t+dt) (q_i(t+dt) - dt * f_i - q_i(t)) - dt * sum_k (d/dq_j(t+dt) g_k * d/dq_i(t+dt) f_k + g_k * d/dq_i(t+dt) d/dq_j(t+dt) f_k )
d/dq_i(t+dt) d/dq_j(t+dt) e = delta_ij - dt * d/dq_j(t+dt) f_i - dt * sum_k (d/dq_j(t+dt) g_k * d/dq_i(t+dt) f_k + g_k * d/dq_i(t+dt) d/dq_j(t+dt) f_k )
d/dq_i(t+dt) d/dq_j(t+dt) e = delta_ij - dt * d/dq_j f_i - dt * sum_k (d/dq_j g_k * d/dq_i f_k + g_k * d/dq_i d/dq_j f_k )
d/dq_i(t+dt) d/dq_j(t+dt) e = delta_ij - dt * d/dq_j f_i - dt * sum_k (d/dq_j (q_k(t+dt) - dt * f_k - q_k(t)) * d/dq_i f_k + g_k * d/dq_i d/dq_j f_k )
d/dq_i(t+dt) d/dq_j(t+dt) e = delta_ij - dt * (d/dq_j f_i + d/dq_i f_j) + dt * sum_k (dt * d/dq_j f_k * d/dq_i f_k - g_k * d/dq_i d/dq_j f_k)

solve for e = 0

q_j(t+dt)[1] = q_j(t)	<- initial guess is last frame

q_j(t+dt)[n+1] = q_j(t+dt)[n] - ||d/dq_j(t+dt) d/dq_k(t+dt) e||^-1 de/dq_j(t+dt)

--]]

local class = require 'ext.class'
local table = require 'ext.table'

local Solver = require 'solver'

local BSSNOK1DBackwardEulerNewton = class(Solver)

BSSNOK1DBackwardEulerNewton.name = 'BSSNOK 1D Backward Euler Newton'

-- fixed dt is required
BSSNOK1DBackwardEulerNewton.fixed_dt = .25

function BSSNOK1DBackwardEulerNewton:init(args)
	args = table(args)
	args.equation = require 'bssnok1d'(args)
	BSSNOK1DBackwardEulerNewton.super.init(self, args)
end

function BSSNOK1DBackwardEulerNewton:iterate()
	local dt = self.fixed_dt
end

return BSSNOK1DBackwardEulerNewton

