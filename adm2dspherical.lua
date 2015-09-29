--[[
variables:
alpha, g_rr, g_hh, A_r, D_rrr, D_rhh, K_rr, K_hh, V_r

definitions:
A_r = (ln alpha),r = alpha,r / alpha
D_kij = g_ij,k/2
V_k = D_km^m - D^m_mk

equations:
alpha,t = -alpha^2 f tr K
g_rr,t = -2 alpha K_rr
g_hh,t = -2 alpha K_hh

A_r,t + (alpha f tr K),r = 0
D_rrr,t + (alpha K_rr),r = 0
D_rhh,t + (alpha K_hh),r = 0
K_rr,t + (alpha lambda^r_rr),r = alpha S_rr
K_hh,t + (alpha lambda^r_hh),r = alpha S_hh
V_r,t = alpha P_r

for:
tr K = g^rr K_rr g^hh K_hh = K_rr / g_rr + K_hh / g_hh
lambda^r_rr = A_r + 2 V_r - 2 D_rhh / g_hh
lambda^r_hh = D_rhh / g_rr
S_rr = K_rr(2 K_hh / g_hh - K_rr/g_rr) + A_r(D_rrr/g_rr - 2 D_rhh/g_hh) + 2 D_rhh/g_hh (D_rrr/g_rr - D_rhh/g_hh) + 2 A_r V_r
S_hh = K_rr K_hh/g_rr - D_rrr D_rhh / g_rr^2 + 1
P_r = -2 / g_hh (A_r K_hh - D_rhh (K_hh / g_hh - K_rr / g_rr))

constraint:
V_r = 2 D_rhh / g_hh

sphrical metric: ds^2 = -dt^2 + dr^2 + r^2 (dh^2 + sin(h)^2 dp^2)
g_uv = diag(-1, 1, r^2, r^2 sin(h)^2)
g^uv = diag(-1, 1, r^-2, r^-2 sin(h)^-2)
... h = theta, p = phi
g_tt is fixed?  or is derived?  is alpha^2?
g_pp is fixed?
either way, none of the skew elements exist so the inverse metric is the diagonal of inverse elements 

alpha, g_rr, g_hh, A_r, D_rrr, D_rhh, K_rr, K_hh, V_r
flux matrix:
[0,	0,	0,	0,	0,	0,	0,	0,	0]
[0,	0,	0,	0,	0,	0,	0,	0,	0]
[0,	0,	0,	0,	0,	0,	0,	0,	0]
[(f + df/dalpha alpha) tr K, -alpha f K_rr/g_rr^2, -alpha f K_hh/g_hh^2, 0, 0, 0, alpha f/g_rr, alpha f/g_hh, 0]
[K_rr, 0, 0, 0, 0, 0, alpha, 0, 0]
[K_hh, 0, 0, 0, 0, 0, 0, alpha, 0]
[lambda^r_rr, 0, 2 alpha D_rhh/g_hh^2, alpha, 0, -2 alpha/g_hh, 0, 0, 2 alpha]
[lambda^r_hh, -alpha/g_rr^2, 0, 0, 0, alpha/g_rr, 0, 0, 0]
[0,	0,	0,	0,	0,	0,	0,	0,	0]

constraint:
V_r = 2 D_rhh / g_hh

variables: alpha, g_rr, g_hh, A_r, D_rrr, D_rhh, K_rr, K_hh, V_r

eigenvalues:
-alpha*sqrt(f/g_rr), -alpha/sqrt(g_rr), 0 x5, alpha/sqrt(g_rr), alpha*sqrt(f/g_rr), 

for lambda = -alpha*sqrt(f/g_rr)
	omega = sqrt(f*g_rr)*(K_rr/g_rr + K_hh/g_hh) - (A_r + 2*V_r)
	v = [0, 0, 0, -1, 0, 0, sqrt(f/g_rr), sqrt(f*g_rr)/g_hh, -2]

for lambda = -alpha/sqrt(g_rr)
	omega = sqrt(g_rr)*K_hh - D_rhh
	v = [0, 0, 0, 0, 0, -1, 0, sqrt(g_rr), 0]

for lambda = 0
	omega = alpha
	omega = g_rr
	omega = g_hh
	omega = V_r
	omega = A_r - f * (D_rrr/g_rr + D_rhh/g_hh)
	v = [1,0,0,0,0,0,0,0,0]
	v = [0,1,0,0,0,0,0,0,0]
	v = [0,0,1,0,0,0,0,0,0]
	v = [0,0,0,0,0,0,0,0,1]
	v = [0,0,0,1,-f/g_rr,-f/g_hh,0,0,0]

for lambda = alpha/sqrt(g_rr)
	omega = sqrt(g_rr)*K_hh + D_rhh
	v = [0, 0, 0, 0, 0, 1, 0, sqrt(g_rr), 0]

for lambda = alpha*sqrt(f/g_rr)
	omega = sqrt(f*g_rr)*(K_rr/g_rr + K_hh/g_hh) + (A_r + 2*V_r)
	v = [0, 0, 0, 1, 0, 0, sqrt(f/g_rr), sqrt(f*g_rr)/g_hh, 2]

--]]
local class = require 'ext.class'
local table = require 'ext.table'
local Equation = require 'equation'

local ADM2DSpherical = class(Equation)
	
ADM2DSpherical.numStates = 9

-- initial conditions
function ADM2DSpherical:init(args, ...)
	local symmath = require 'symmath'

	local r = assert(args.r)

	local h = symmath.clone(assert(args.h)):simplify()
	self.calc_h = h:compile{r}

	local dr_h = h:diff(r):simplify()
	self.calc_dr_h = dr_h:compile{r}

	local d2r_h = dr_h:diff(r):simplify()
	self.calc_d2r_h = d2r_h:compile{r}

	local g_rr = symmath.clone(assert(args.g_rr)):simplify()
	self.calc_g_rr = g_rr:compile{r}

	local dr_g_rr = g_rr:diff(r):simplify()
	self.calc_dr_g_rr = dr_g_rr:compile{r}

	local g_hh = symmath.clone(assert(args.g_hh)):simplify()
	self.calc_g_hh = g_hh:compile{r}

	local dr_g_hh = g_hh:diff(r):simplify()
	self.calc_dr_g_hh = dr_g_hh:compile{r}

	local alpha = symmath.clone(assert(args.alpha)):simplify()
	self.calc_alpha = alpha:compile{r}

	local dr_alpha = alpha:diff(r):simplify()
	self.calc_dr_alpha = dr_alpha:compile{r}

	local f = symmath.clone(assert(args.alpha)):simplify()
	self.calc_f = f:compile{assert(args.f_param)}

	local dalpha_f = f:diff(args.f_param):simplify()
	self.calc_dalpha_f = dalpha_f:compile{args.f_param}

	-- alpha, g_rr, g_hh/r^2
	-- A_r, D_rrr, D_rhh/r
	-- K_rr, K_hh, r V_r
	local get_r = function(self,i) return self.xs[i] end
	local get_state = function(self, i) return self.qs[i] end
	local get_alpha = get_state:index(1)
	local get_g_rr = get_state:index(2)
	local get_g_hh = get_state:index(3)
	local get_A_r = get_state:index(4)
	local get_D_rrr = get_state:index(5)
	local get_D_rhh = get_state:index(6)
	local get_K_rr = get_state:index(7)
	local get_K_hh = get_state:index(8)
	local get_V_r = get_state:index(9)
	self.graphInfos = table{
		{viewport = {0/4, 0/3, 1/4, 1/3}, getter = get_alpha, name = 'alpha', color = {1,0,1}},
		{viewport = {1/4, 0/3, 1/4, 1/3}, getter = get_g_rr, name = 'g_rr', color = {0,1,0}},
		{viewport = {2/4, 0/3, 1/4, 1/3}, getter = get_g_hh / get_r^2, name = 'g_hh/r^2', color = {0,1,0}},
		{viewport = {0/4, 1/3, 1/4, 1/3}, getter = get_A_r, name = 'A_r', color = {0,1,0}},
		{viewport = {1/4, 1/3, 1/4, 1/3}, getter = get_D_rrr, name = 'D_rrr', color = {.5,.5,1}},
		{viewport = {2/4, 1/3, 1/4, 1/3}, getter = get_D_rhh / get_r, name = 'D_rhh/r', color = {1,1,0}},
		{viewport = {0/4, 2/3, 1/4, 1/3}, getter = get_K_rr, name = 'K_rr', color = {0,1,1}},
		{viewport = {1/4, 2/3, 1/4, 1/3}, getter = get_K_hh, name = 'K_hh', color = {0,1,1}},
		{viewport = {2/4, 2/3, 1/4, 1/3}, getter = get_V_r * get_r, name = 'V_r*r', color = {0,1,1}},
		{viewport = {3/4, 0/3, 1/4, 1/3}, getter = function(self,i) return self.qs[i][1] * self.xs[i] * math.sqrt(self.qs[i][2] * self.qs[i][3]) end, name = 'volume', color = {0,1,1}},
		{viewport = {3/4, 1/3, 1/4, 1/3}, getter = function(self,i) return math.log(self.eigenbasisErrors[i]) end, name = 'log eigenbasis error', color = {1,0,0}, range = {-30, 30}},
		{viewport = {3/4, 2/3, 1/4, 1/3}, getter = function(self,i) return math.log(self.fluxMatrixErrors[i]) end, name = 'log reconstruction error', color = {1,0,0}, range = {-30, 30}},
	}
	self.graphInfoForNames = self.graphInfos:map(function(info,i)
		return info, info.name
	end)
end

function ADM2DSpherical:initCell(sim,i)
	local r = sim.xs[i]
	local alpha = self.calc_alpha(r)
	local g_rr = self.calc_g_rr(r)
	local g_hh = self.calc_g_hh(r)
	local A_r = self.calc_dr_alpha(r)
	local D_rrr = self.calc_dr_g_rr(r)/2
	local D_rhh = self.calc_dr_g_hh(r)/2
	local K_rr = -self.calc_d2r_h(r) / sqrt(g_rr)
	local K_hh = -r * self.calc_dr_h(r) / sqrt(g_rr)
	local V_r = D_rhh / g_hh
	return {alpha, g_rr, g_hh, A_r, D_rrr, D_rhh, K_rr, K_hh, V_r}
end

function ADM2DSpherical:calcInterfaceEigenBasis(sim,i,qL,qR)
	local avgQ = {}
	for j=1,self.numStates do
		avgQ[j] = (qL[j] + qR[j]) / 2
	end
	
	local alpha, g_rr, g_hh, A_r, D_rrr, D_rhh, K_rr, K_hh, V_r = unpack(avgQ)
	local x = sim.ixs[i]
	local f = self.calc_f(alpha)
	local tr_K = K_rr / g_rr + K_hh / g_hh

	local l1 = alpha / sqrt(g_rr)
	local l2 = alpha * sqrt(f / g_rr)
	sim.eigenvalues[i] = {-l2, -l1, 0, 0, 0, 0, 0, l1, l2}

	-- row-major, math-indexed
	local lambdaUr_rr = A_r + 2 * V_r - 2 * D_rhh / g_hh
	local lambdaUr_hh = D_rhh / g_rr
	sim.fluxMatrix[i] = {
		{0,0,0,0,0,0,0,0,0},
		{0,0,0,0,0,0,0,0,0},
		{0,0,0,0,0,0,0,0,0},
		{f * tr_K, -alpha * f * K_rr/g_rr^2, -alpha * f * K_hh/g_hh^2, 0, 0, 0, alpha * f / g_rr, alpha * f / g_hh, 0},
		{K_rr, 0, 0, 0, 0, 0, alpha, 0, 0},
		{K_hh, 0, 0, 0, 0, 0, 0, alpha, 0},
		{lambdaUr_rr, 0, 2 * alpha * D_rhh/g_hh^2, alpha, 0, -2 * alpha/g_hh, 0, 0, 2 * alpha},
		{lambdaUr_hh, -alpha/g_rr^2, 0, 0, 0, alpha/g_rr, 0, 0, 0},
		{0,0,0,0,0,0,0,0,0},
	}
	sim.eigenvectors[i] = {
		{0,0,1,0,0,0,0,0,0},
		{0,0,0,1,0,0,0,0,0},
		{0,0,0,0,1,0,0,0,0},
		{-1/2,0,0,0,0,-2,0,0,1/2},
		{-g_rr/(2*f),g_rr/(2*g_hh),0,0,0,-(2*g_rr)/f,-1/f,-g_rr/(2*g_hh),g_rr/(2*f)},
		{0,-1/2,0,0,0,0,0,1/2,0},
		{.5*sqrt(g_rr/f),-sqrt(g_rr)/(2*g_hh),0,0,0,0,0,-sqrt(g_rr)/(2*g_hh),.5*sqrt(g_rr/f)},
		{0,1/(2*sqrt(g_rr)),0,0,0,0,0,1/(2*sqrt(g_rr)),0},
		{0,0,0,0,0,1,0,0,0}
	}
	sim.eigenvectorsInverse[i] = {
		{0, 0, 0, -1, 0, 0, sqrt(f/g_rr), sqrt(f*g_rr)/g_hh, -2},
		{0, 0, 0, 0, 0, -1, 0, sqrt(g_rr), 0},
		{1, 0, 0, 0, 0, 0, 0, 0, 0}, 
		{0, 1, 0, 0, 0, 0, 0, 0, 0}, 
		{0, 0, 1, 0, 0, 0, 0, 0, 0}, 
		{0, 0, 0, 0, 0, 0, 0, 0, 1}, 
		{0, 0, 0, g_rr, -f, -f*g_rr/g_hh, 0, 0, 0}, 
		{0, 0, 0, 0, 0, 1, 0, sqrt(g_rr), 0},
		{0, 0, 0, 1, 0, 0, sqrt(f/g_rr), sqrt(f*g_rr)/g_hh, 2}	
	}
end

function ADM2DSpherical:sourceTerm(sim)
	local source = sim:newState()
	for i=1,sim.gridsize do
		local alpha, g_rr, g_hh, A_r, D_rrr, D_rhh, K_rr, K_hh, V_r = unpack(sim.qs[i])
		local f = self.calc_f(alpha)
		local tr_K = K_rr / g_rr + K_hh / g_hh
		local S_rr = K_rr * (2 * K_hh / g_hh - K_rr / g_rr) + A_r * (D_rrr / g_rr - 2 * D_rhh / g_hh)
					+ 2 * D_rhh / g_hh * (D_rrr / g_rr - D_rhh / g_hh) + 2 * A_r * V_r
		local S_hh = K_rr * K_hh / g_rr - D_rrr * D_rhh / g_rr^2 + 1
		local P_r = -2 / g_hh * (A_r * K_hh - D_rhh * (K_hh / g_hh - K_rr / g_rr))
		source[i][1] = -alpha * alpha * f * tr_K 
		source[i][2] = -2 * alpha * K_rr
		source[i][3] = -2 * alpha * K_hh
		source[i][7] = alpha * S_rr
		source[i][8] = alpha * S_hh
		source[i][9] = alpha * P_r
	end
	return source
end

return ADM2DSpherical
