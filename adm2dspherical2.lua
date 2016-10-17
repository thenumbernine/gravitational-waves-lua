--[[
dl^2 = A(r,t) dr^2 + r^2 B(r,t) dOmega^2
dOmega^2 = dtheta^2 + sin(theta)^2 dphi^2
DA = (ln A),r
DB = (ln B),r
Dalpha = (ln alpha),r
MA = 2 SB - SA - rho
MB = SA - rho
KA = K^r_r
KB = K^theta_theta
SA = S^r_r
SB = S^theta_theta

A,t= -2 alpha A KA
B,t= -2 alpha B KB
DA,t = -2 alpha (KA Dalpha + KA,r)
DB,t = -2 alpha (KB Dalpha + KB,r)
KA,t = -alpha / A ( (Dalpha + DB),r + Dalpha^2 - Dalpha DA / 2 + DB^2 - DA DB / 2 - A KA (KA + 2 KB) - 1/r (KA - 2 KB) ) + 4 pi alpha MA
KB,t = -alpha / (2A) (DB,r + Dalpha DB + DB^2 - DA DB / 2 - 1/r (DA - 2 Dalpha - 4 DB) - 2 (A - B) / (r^2 B) ) + alpha KB (KA + 2 KB) + 4 pi alpha MB

H = -DB,r + 1/(r^2 B) (A - B) + A KB (2 KA + KB) + 1/r (DA - 3 DB) + DA DB / 2 - 3 DB^2 / 4 - 8 pi A rho = 0
M = -KB,r + (KA - KB) (1/r + DB/2) - 4 pi jA = 0
jA = jr = momentum density of matter in radial direction

A^0 = B^0
K^0_A = K^0_B
lambda = 1/r (1 - A/B)
KB,t = -alpha / (2A) (DB,r + Dalpha DB + DB^2 - DA DB / 2 - 1/r (DA - 2 Dalpha - 4 DB) + 2 lambda / r) + alpha KB (KA + 2 KB) + 4 pi alpha MB
H = -DB,r - lambda / r + A KB (2 KA + KB) + 1/r (DA - 3 DB) + DA DB / 2 - 3 DB^2/4 - 8 pi A rho
lambda,t = 2 alpha A / B (KA - KB) / r
	... substitute momentum constraint ...
lambda,t = 2 alpha A / B (KB,r - DB / 2 ( KA - KB) + 4 pi jA)
alpha,t = -alpha^2 f K = -alpha^2 f (KA + 2 KB)
Dalpha,t = -alpha,r f (KA + 2 KB)

	system:
Dalpha,t = -(alpha f (KA + 2 KB)),r
DA,t = -2 alpha (KA Dalpha + KA,r)
DB,t = -2 alpha (KB Dalpha + KB,r)
KA,t = -alpha / A ((Dalpha + DB),r + Dalpha^2 - Dalpha DA / 2 + DB^2 - DA DB / 2 - A KA (KA + 2 KB) - 1/r (KA - 2 KB) ) + 4 pi alpha MA
KB,t = -alpha / (2A) (DB,r + Dalpha DB + DB^2 - DA DB / 2 - 1/r (DA - 2 Dalpha - 4 DB) - 2 (A - B) / (r^2 B) ) + alpha KB (KA + 2 KB) + 4 pi alpha MB
KB,t = -alpha / (2A) (DB,r + Dalpha DB + DB^2 - DA DB / 2 - 1/r (DA - 2 Dalpha - 4 DB) + 2 lambda / r) + alpha KB (KA + 2 KB) + 4 pi alpha MB

	eigenvalues:
lambda=0, omega_0 = Dalpha / f - (DA + 2 DB) / 2
lambda=+-alpha / sqrt(A), omega^l_+- = sqrt(A) KB -+ DB / 2
lambda=+-alpha sqrt(f / A), omega^f_+- = sqrt(A) (KA + 2 (f + 1) / (f - 1) KB) -+ (Dalpha / sqrt(f) + 2 DB / (f - 1))

DTilde = DA - 2 DB
K = KA + 2 KB

DTilde,t = -2 alpha (K,r + Dalpha (K - 4 KB) - 4 (K - 3 KB) (1/r + DB/2) + 16 pi jA)
K,t = -alpha / A (Dalpha,r + Dalpha^2 + 2 Dalpha/r - Dalpha DTilde / 2) + alpha (K^2 - 4 K KB + 6 KB^2) + 4 pi alpha (rho + SA + 2 SB)

	now the eigenvalues are:
lambda = 0, omega_0=Dalpha/f - DTilde/2
lambda = +-alpha/sqrt(A), omega^l_+- = sqrt(A) KB -+ DB/2
lambda = +-alpha sqrt(f/A), omega^f_+- = sqrt(A) K -+ Dalpha / sqrt(f) ...

--]]

local class = require 'ext.class'
local table = require 'ext.table'

local Roe = require 'roe'


local ADM2DSpherical2 = class()

function ADM2DSpherical2:initCell(sim,i)
	local DAlpha
	local DA
	local DB
	local KA
	local KB
	return {Dalpha, DA, DB, KA, KB}
end

function ADM2DSpherical2:calcInterfaceEigenBasis(sim,i,qL,qR)
	local avgQ = {}
	for j=1,self.numStates do 
		avgQ[j] = (qL[j] + qR[j]) / 2
	end
	local f = self.calc.f(alpha)
	local Dalpha, DA, DB, KA, KB = unpack(avgQ)
	sim.eigenvalues[i] = {
		-alpha * sqrt(f/A)
		-alpha / sqrt(A)
		0,
		alpha / sqrt(A),
		alpha * sqrt(f/A)
	}
end

local ADM2DSpherical2Roe = class(Roe)

function ADM2DSpherical2Roe:init(args)
	args = table(args)
	args.equation = ADM2DSpherical2(args)
	ADM2DSpherical2Roe.super.init(self, args)
end


