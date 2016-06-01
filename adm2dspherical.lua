--[[
variables:
alpha, gamma_rr, gamma_hh, A_r, D_rrr, D_rhh, K_rr, K_hh, V_r

definitions:
A_r = (ln alpha),r = alpha,r / alpha
D_kij = gamma_ij,k/2
V_k = D_km^m - D^m_mk

equations:
alpha,t = -alpha^2 f tr K
gamma_rr,t = -2 alpha K_rr
gamma_hh,t = -2 alpha K_hh

A_r,t + (alpha f tr K),r = 0
D_rrr,t + (alpha K_rr),r = 0
D_rhh,t + (alpha K_hh),r = 0
K_rr,t + (alpha lambda^r_rr),r = alpha S_rr
K_hh,t + (alpha lambda^r_hh),r = alpha S_hh
V_r,t = alpha P_r

for:
tr K = gamma^rr K_rr gamma^hh K_hh = K_rr / gamma_rr + K_hh / gamma_hh
lambda^r_rr = A_r + 2 V_r - 2 D_rhh / gamma_hh
lambda^r_hh = D_rhh / gamma_rr
S_rr = K_rr(2 K_hh / gamma_hh - K_rr/gamma_rr) + A_r(D_rrr/gamma_rr - 2 D_rhh/gamma_hh) + 2 D_rhh/gamma_hh (D_rrr/gamma_rr - D_rhh/gamma_hh) + 2 A_r V_r
S_hh = K_rr K_hh/gamma_rr - D_rrr D_rhh / gamma_rr^2 + 1
P_r = -2 / gamma_hh (A_r K_hh - D_rhh (K_hh / gamma_hh - K_rr / gamma_rr))

constraint:
V_r = 2 D_rhh / gamma_hh

sphrical metric: ds^2 = -dt^2 + dr^2 + r^2 (dh^2 + sin(h)^2 dp^2)
g_uv = diag(-1, 1, r^2, r^2 sin(h)^2)
g^uv = diag(-1, 1, r^-2, r^-2 sin(h)^-2)
... h = theta, p = phi
g_tt is fixed?  or is derived?  is alpha^2?
g_pp is fixed?
either way, none of the skew elements exist so the inverse metric is the diagonal of inverse elements 

alpha, gamma_rr, gamma_hh, A_r, D_rrr, D_rhh, K_rr, K_hh, V_r
flux matrix:
[0,	0,	0,	0,	0,	0,	0,	0,	0]
[0,	0,	0,	0,	0,	0,	0,	0,	0]
[0,	0,	0,	0,	0,	0,	0,	0,	0]
[(f + df/dalpha alpha) tr K, -alpha f K_rr/gamma_rr^2, -alpha f K_hh/gamma_hh^2, 0, 0, 0, alpha f/gamma_rr, alpha f/gamma_hh, 0]
[K_rr, 0, 0, 0, 0, 0, alpha, 0, 0]
[K_hh, 0, 0, 0, 0, 0, 0, alpha, 0]
[lambda^r_rr, 0, 2 alpha D_rhh/gamma_hh^2, alpha, 0, -2 alpha/gamma_hh, 0, 0, 2 alpha]
[lambda^r_hh, -alpha/gamma_rr^2, 0, 0, 0, alpha/gamma_rr, 0, 0, 0]
[0,	0,	0,	0,	0,	0,	0,	0,	0]

constraint:
V_r = 2 D_rhh / gamma_hh

variables: alpha, gamma_rr, gamma_hh, A_r, D_rrr, D_rhh, K_rr, K_hh, V_r

eigenvalues:
-alpha*sqrt(f/gamma_rr), -alpha/sqrt(gamma_rr), 0 x5, alpha/sqrt(gamma_rr), alpha*sqrt(f/gamma_rr), 

for lambda = -alpha*sqrt(f/gamma_rr)
	omega = sqrt(f*gamma_rr)*(K_rr/gamma_rr + K_hh/gamma_hh) - (A_r + 2*V_r)
	v = [0, 0, 0, -1, 0, 0, sqrt(f/gamma_rr), sqrt(f*gamma_rr)/gamma_hh, -2]

for lambda = -alpha/sqrt(gamma_rr)
	omega = sqrt(gamma_rr)*K_hh - D_rhh
	v = [0, 0, 0, 0, 0, -1, 0, sqrt(gamma_rr), 0]

for lambda = 0
	omega = alpha
	omega = gamma_rr
	omega = gamma_hh
	omega = V_r
	omega = A_r - f * (D_rrr/gamma_rr + D_rhh/gamma_hh)
	v = [1,0,0,0,0,0,0,0,0]
	v = [0,1,0,0,0,0,0,0,0]
	v = [0,0,1,0,0,0,0,0,0]
	v = [0,0,0,0,0,0,0,0,1]
	v = [0,0,0,1,-f/gamma_rr,-f/gamma_hh,0,0,0]

for lambda = alpha/sqrt(gamma_rr)
	omega = sqrt(gamma_rr)*K_hh + D_rhh
	v = [0, 0, 0, 0, 0, 1, 0, sqrt(gamma_rr), 0]

for lambda = alpha*sqrt(f/gamma_rr)
	omega = sqrt(f*gamma_rr)*(K_rr/gamma_rr + K_hh/gamma_hh) + (A_r + 2*V_r)
	v = [0, 0, 0, 1, 0, 0, sqrt(f/gamma_rr), sqrt(f*gamma_rr)/gamma_hh, 2]

--]]
local class = require 'ext.class'
local table = require 'ext.table'
local Equation = require 'equation'
local symmath = require 'symmath'

local ADM2DSpherical = class(Equation)
	
ADM2DSpherical.numStates = 9

-- initial conditions
function ADM2DSpherical:init(args, ...)
	local r = assert(args.r)

	local h = symmath.clone(assert(args.h)):simplify()
	self.calc_h = h:compile{r}

	local dr_h = h:diff(r):simplify()
	self.calc_dr_h = dr_h:compile{r}

	local d2r_h = dr_h:diff(r):simplify()
	self.calc_d2r_h = d2r_h:compile{r}

	local gamma_rr = symmath.clone(assert(args.gamma_rr)):simplify()
	self.calc_gamma_rr = gamma_rr:compile{r}

	local dr_gamma_rr = gamma_rr:diff(r):simplify()
	self.calc_dr_gamma_rr = dr_gamma_rr:compile{r}

	local gamma_hh = symmath.clone(assert(args.gamma_hh)):simplify()
	self.calc_gamma_hh = gamma_hh:compile{r}

	local dr_gamma_hh = gamma_hh:diff(r):simplify()
	self.calc_dr_gamma_hh = dr_gamma_hh:compile{r}

	local alpha = symmath.clone(assert(args.alpha)):simplify()
	self.calc_alpha = alpha:compile{r}

	local dr_alpha = alpha:diff(r):simplify()
	self.calc_dr_alpha = dr_alpha:compile{r}

	local f = symmath.clone(assert(args.alpha)):simplify()
	self.calc_f = f:compile{assert(args.f_param)}

	local dalpha_f = f:diff(args.f_param):simplify()
	self.calc_dalpha_f = dalpha_f:compile{args.f_param}

	-- alpha, gamma_rr, gamma_hh/r^2
	-- A_r, D_rrr, D_rhh/r
	-- K_rr, K_hh, r V_r
	local r = function(self,i) return self.xs[i] end
	local q = function(self,i) return self.qs[i] end
	local alpha = q:_(1)
	local gamma_rr = q:_(2)
	local gamma_hh = q:_(3)
	local A_r = q:_(4)
	local D_rrr = q:_(5)
	local D_rhh = q:_(6)
	local K_rr = q:_(7)
	local K_hh = q:_(8)
	local V_r = q:_(9)
	local volume = alpha * r * math.sqrt:o(gamma_rr * gamma_hh)
	self:buildGraphInfos{
		{alpha = alpha},
		{gamma_rr = gamma_rr},
		{['gamma_hh/r^2'] = gamma_hh / r^2},
		{A_r = A_r},
		{D_rrr = D_rrr},
		{['D_rhh/r'] = D_rhh / r},
		{K_rr = K_rr},
		{K_hh = K_hh},
		{['V_r*r'] = V_r * r},
		{volume = volume},
		{['log eigenbasis error'] = function(self,i) return math.log(self.eigenbasisErrors[i]) end},
		{['log reconstruction error'] = function(self,i) return math.log(self.fluxMatrixErrors[i]) end},
	}
end

function ADM2DSpherical:initCell(sim,i)
	local r = sim.xs[i]
	local alpha = self.calc_alpha(r)
	local gamma_rr = self.calc_gamma_rr(r)
	local gamma_hh = self.calc_gamma_hh(r)
	local A_r = self.calc_dr_alpha(r)
	local D_rrr = self.calc_dr_gamma_rr(r)/2
	local D_rhh = self.calc_dr_gamma_hh(r)/2
	local K_rr = -self.calc_d2r_h(r) / sqrt(gamma_rr)
	local K_hh = -r * self.calc_dr_h(r) / sqrt(gamma_rr)
	local V_r = D_rhh / gamma_hh
	return {alpha, gamma_rr, gamma_hh, A_r, D_rrr, D_rhh, K_rr, K_hh, V_r}
end

function ADM2DSpherical:calcInterfaceEigenBasis(sim,i,qL,qR)
	local avgQ = {}
	for j=1,self.numStates do
		avgQ[j] = (qL[j] + qR[j]) / 2
	end
	
	local alpha, gamma_rr, gamma_hh, A_r, D_rrr, D_rhh, K_rr, K_hh, V_r = unpack(avgQ)
	local x = sim.ixs[i]
	local f = self.calc_f(alpha)
	local tr_K = K_rr / gamma_rr + K_hh / gamma_hh

	local l1 = alpha / sqrt(gamma_rr)
	local l2 = alpha * sqrt(f / gamma_rr)
	sim.eigenvalues[i] = {-l2, -l1, 0, 0, 0, 0, 0, l1, l2}

	-- row-major, math-indexed
	local lambdaUr_rr = A_r + 2 * V_r - 2 * D_rhh / gamma_hh
	local lambdaUr_hh = D_rhh / gamma_rr
	sim.fluxMatrix[i] = {
		{0,0,0,0,0,0,0,0,0},
		{0,0,0,0,0,0,0,0,0},
		{0,0,0,0,0,0,0,0,0},
		{f * tr_K, -alpha * f * K_rr/gamma_rr^2, -alpha * f * K_hh/gamma_hh^2, 0, 0, 0, alpha * f / gamma_rr, alpha * f / gamma_hh, 0},
		{K_rr, 0, 0, 0, 0, 0, alpha, 0, 0},
		{K_hh, 0, 0, 0, 0, 0, 0, alpha, 0},
		{lambdaUr_rr, 0, 2 * alpha * D_rhh/gamma_hh^2, alpha, 0, -2 * alpha/gamma_hh, 0, 0, 2 * alpha},
		{lambdaUr_hh, -alpha/gamma_rr^2, 0, 0, 0, alpha/gamma_rr, 0, 0, 0},
		{0,0,0,0,0,0,0,0,0},
	}
	sim.eigenvectors[i] = {
		{0,0,1,0,0,0,0,0,0},
		{0,0,0,1,0,0,0,0,0},
		{0,0,0,0,1,0,0,0,0},
		{-1/2,0,0,0,0,-2,0,0,1/2},
		{-gamma_rr/(2*f),gamma_rr/(2*gamma_hh),0,0,0,-(2*gamma_rr)/f,-1/f,-gamma_rr/(2*gamma_hh),gamma_rr/(2*f)},
		{0,-1/2,0,0,0,0,0,1/2,0},
		{.5*sqrt(gamma_rr/f),-sqrt(gamma_rr)/(2*gamma_hh),0,0,0,0,0,-sqrt(gamma_rr)/(2*gamma_hh),.5*sqrt(gamma_rr/f)},
		{0,1/(2*sqrt(gamma_rr)),0,0,0,0,0,1/(2*sqrt(gamma_rr)),0},
		{0,0,0,0,0,1,0,0,0}
	}
	sim.eigenvectorsInverse[i] = {
		{0, 0, 0, -1, 0, 0, sqrt(f/gamma_rr), sqrt(f*gamma_rr)/gamma_hh, -2},
		{0, 0, 0, 0, 0, -1, 0, sqrt(gamma_rr), 0},
		{1, 0, 0, 0, 0, 0, 0, 0, 0}, 
		{0, 1, 0, 0, 0, 0, 0, 0, 0}, 
		{0, 0, 1, 0, 0, 0, 0, 0, 0}, 
		{0, 0, 0, 0, 0, 0, 0, 0, 1}, 
		{0, 0, 0, gamma_rr, -f, -f*gamma_rr/gamma_hh, 0, 0, 0}, 
		{0, 0, 0, 0, 0, 1, 0, sqrt(gamma_rr), 0},
		{0, 0, 0, 1, 0, 0, sqrt(f/gamma_rr), sqrt(f*gamma_rr)/gamma_hh, 2}	
	}
end

function ADM2DSpherical:calcEigenvaluesFromCons(alpha, gamma_rr, ...) 
	local f = self.calc_f(alpha)
	local l1 = alpha / math.sqrt(gamma_rr)
	local l2 = alpha * math.sqrt(f / gamma_rr)
	return -l2, -l1, 0, 0, 0, 0, 0, l1, l2

end

function ADM2DSpherical:sourceTerm(sim, qs)
	local source = sim:newState()
	for i=1,sim.gridsize do
		local alpha, gamma_rr, gamma_hh, A_r, D_rrr, D_rhh, K_rr, K_hh, V_r = unpack(qs[i])
		local f = self.calc_f(alpha)
		local tr_K = K_rr / gamma_rr + K_hh / gamma_hh
		local S_rr = K_rr * (2 * K_hh / gamma_hh - K_rr / gamma_rr) + A_r * (D_rrr / gamma_rr - 2 * D_rhh / gamma_hh)
					+ 2 * D_rhh / gamma_hh * (D_rrr / gamma_rr - D_rhh / gamma_hh) + 2 * A_r * V_r
		local S_hh = K_rr * K_hh / gamma_rr - D_rrr * D_rhh / gamma_rr^2 + 1
		local P_r = -2 / gamma_hh * (A_r * K_hh - D_rhh * (K_hh / gamma_hh - K_rr / gamma_rr))
		source[i][1] = -alpha * alpha * f * tr_K 
		source[i][2] = -2 * alpha * K_rr
		source[i][3] = -2 * alpha * K_hh
		source[i][7] = alpha * S_rr
		source[i][8] = alpha * S_hh
		source[i][9] = alpha * P_r
	end
	return source
end

function ADM2DSpherical:postIterate(sim, qs)
	-- enforce constraint of V_r = 2 * D_rhh / gamma_hh
	-- i.e. 1 = 2 * D_rhh / (gamma_hh * V_r)
	for i=1,sim.gridsize do
		--[[ direct assign:
		qs[i][9] = 2 * qs[i][6] / qs[i][3]
		--]]
		-- [[ geometric renormalization
		local c = abs(2 * qs[i][6] / (qs[i][3] * qs[i][9]))^(1/3)
		qs[i][3] = qs[i][3] * c
		qs[i][6] = qs[i][6] / c
		qs[i][9] = qs[i][9] * c
		-- 2 * (qs[i][6] / cbrt(c)) / (qs[i][3] * cbrt(c) * qs[i][9] * cbrt(c))
		-- = (2 * qs[i][6] / (qs[i][3] * qs[i][9])) / c
		-- = 1
	end
end

return ADM2DSpherical
