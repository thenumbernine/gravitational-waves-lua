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

... evolution equations become ...
alpha,t = -alpha^2 f tr K
g_rr,t = -2 alpha K_rr
g_hh,t = -2 alpha K_hh
A_r,t 
	+ alpha,r (f tr K) 
	+ f,r (alpha tr K) 	... which equals ... alpha,r (df/dalpha alpha tr K)
	+ g_rr,r (-alpha f K_rr / g_rr^2) 
	+ g_hh,r (-alpha f K_hh / g_hh^2) 
	+ K_rr,r (alpha f / g_rr)
	+ K_hh,r (alpha f / g_hh)
	= 0
D_rrr,t + alpha,r K_rr + K_rr,r alpha = 0
D_rhh,t + alpha,r K_hh + K_hh,r alpha = 0
K_rr,t 
	+ alpha,r lambda^r_rr 
	+ g_hh,r (2 alpha D_rhh / g_hh^2)
	+ A_r,r alpha
	+ D_rhh,r (-2 alpha / g_hh)
	+ V_r,r (2 alpha)
	= alpha S_rr
K_hh,t 
	+ alpha,r lambda^r_hh 
	+ g_rr,r (-alpha / g_rr^2)
	+ D_rhh,r (alpha / g_rr)
	= alpha S_hh
V_r,t = alpha P_r

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

eigenvalues:
-alpha*sqrt(f/g_rr), -alpha/sqrt(g_rr), 0 x5, alpha/sqrt(g_rr), alpha*sqrt(f/g_rr), 

eigenvector columns:
lambda = -alpha*sqrt(f/g_rr)
v = [0,0,0,f,g_rr,0,-sqrt(f g_rr),0,0]
v is equal to sqrt(g_rr) K_rr - D_rrr g_rr / sqrt(f) - A_r sqrt(f)
Alcubierre says sqrt(g_rr) K_hh - D_rhh == [0,0,0,0,0,-1,0,sqrt(g_rr),0]
so does D_rrr g_rr / sqrt(f) + A_r sqrt(f) == D_rhh ?
	does g_rr,r g_rr / sqrt(f) + 2 alpha,r sqrt(f) / alpha == g_hh,r
	does g_rr,r (g_rr / sqrt(f) - 1) + 2 alpha,r sqrt(f) / alpha == 0?
also which of many forms should Alcubierre's eqn be reformulated into a vector?

lambda = -alpha/sqrt(g_rr)
v = [0,0,0,f,-(f-2)*g_rr,(f-1)*g_hh,(f-2)*sqrt(g_rr),-((f-1)*g_hh)/sqrt(g_rr),0]
Alcubierre says sqrt(f g_rr) tr K - (A_r + 2 V_r) = [0, 0, 0, -1, 0, 0, sqrt(f/g_rr), sqrt(f g_rr) / g_hh, -2]

lambda = 0
Alcubierre says alpha, g_rr, g_hh, V_r, A_r - f D_r^m_m

lambda = alpha/sqrt(g_rr)
v = [0,0,0,f,-(f-2)*g_rr,(f-1)*g_hh,-(f-2)*sqrt(g_rr),((f-1)*g_hh)/sqrt(g_rr),0]
Alcubierre says says sqrt(f g_rr) tr K + (A_r + 2 V_r)

lambda = alpha*sqrt(f/g_rr)
v = [0,0,0,f,g_rr,0,sqrt(f g_rr),0,0]
... which amounts to A_r sqrt(f) + D_rrr g_rr / sqrt(f) + K_rr sqrt(g_rr)
Alcubierre says sqrt(g_rr) K_hh + D_rhh == [0,0,0,0,0,1,0,sqrt(g_rr),0]

--]]
require 'ext'
local Simulation = require 'adm1d.simulation'
local ADM2DSim = class(Simulation)
	
ADM2DSim.numStates = 9

-- initial conditions
function ADM2DSim:init(...)
	Simulation.init(self, ...)

	self.graphInfos = {
		{viewport={0/3, 0/3, 1/3, 1/3}, getter=index:bind(self.qs):index(1), name='alpha', color={1,0,1}},
		{viewport={0/3, 1/3, 1/3, 1/3}, getter=index:bind(self.qs):index(3), name='A', color={0,1,0}},
		{viewport={1/3, 0/3, 1/3, 1/3}, getter=index:bind(self.qs):index(2), name='g', color={.5,.5,1}},
		{viewport={1/3, 1/3, 1/3, 1/3}, getter=index:bind(self.qs):index(4), name='D', color={1,1,0}},
		{viewport={2/3, 0/3, 1/3, 1/3}, getter=index:bind(self.qs):index(5), name='K', color={0,1,1}},
		{viewport={0/3, 2/3, 1/3, 1/3}, getter=log:compose(index:bind(self.eigenbasisErrors)), name='error', color={1,0,0}, range={-30, 30}},
		{viewport={1/3, 2/3, 1/3, 1/3}, getter=log:compose(index:bind(self.fluxMatrixErrors)), name='error', color={1,0,0}, range={-30, 30}},
	}
end

do
	local sigma = 10
	local xc = 150
	local H = 5
	-- aux var for init
	local function calc_h(r) return H * exp(-(r - xc)^2 / sigma^2) end
	local function d_calc_h(r) return -2 * (r - xc) / sigma^2 * calc_h(r) end
	local function d2_calc_h(r) return (-2 / sigma^2 + 4 * (r - xc)^2 / sigma^4) * calc_h(r) end
	local function calc_g_rr(r) return 1 - d_calc_h(r)^2 end
	local function dr_calc_g_rr(r) return -2 * d_calc_h(r) * d2_calc_h(r) end
	local function calc_g_hh(r) return r^2 end
	local function dr_calc_g_hh(r) return 2*r end
	local function calc_alpha(r) return 1 end
	local function dr_calc_alpha(r) return 0 end

	function ADM2DSim:initCell(i)
		local r = self.xs[i]
		local alpha = calc_alpha(r)
		local g_rr = calc_g_rr(r)
		local g_hh = calc_g_hh(r)
		local A_r = dr_calc_alpha(r)
		local D_rrr = dr_calc_g_rr(r)/2
		local D_rhh = dr_calc_g_hh(r)/2
		local K_rr = -d2_calc_h(r) / sqrt(g_rr)
		local K_hh = -r * d_calc_h(r) / sqrt(g_rr)
		local V_r = D_rhh / g_hh
		return {alpha, g_rr, g_hh, A_r, D_rrr, D_rhh, K_rr, K_hh, V_r}
	end
end

function ADM2DSim:calcInterfaceEigenBasis(i)
	local avgQ = {}
	for j=1,self.numStates do
		avgQ[j] = (self.qs[i-1][j] + self.qs[i][j]) / 2
	end
	
	local alpha, g_rr, g_hh, A_r, D_rrr, D_rhh, K_rr, K_hh, V_r = unpack(avgQ)
	local x = self.ixs[i]
	local f = self.calc_f(alpha)
	
	local l1 = alpha / sqrt(g_rr)
	local l2 = l1 * (f)
	
	local lambdaUr_rr = A_r + 2 * V_r - 2 * D_rhh / g_hh
	local lambdaUr_hh = D_rhh / g_rr
	
	self.eigenvalues[i] = {-l2, -l1, 0, 0, 0, 0, 0, l1, l2}
	-- row-major, math-indexed
	self.fluxMatrix[i] = {
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
	-- TODO figure these out
	self.eigenvectors[i] = {
		{0,			0,					1,	0,	0,	0,	0,	0,					0,			},
		{0,			0,					0,	1,	0,	0,	0,	0,					0,			},
		{0,			0,					0,	0,	1,	0,	0,	0,					0,			},
		{0,			-1,					0,	0,	0,	0,	0,	1,					0,			},
		{0,			0,					0,	0,	0,	1,	0,	0,					0,			},
		{-1,		0,					0,	0,	0,	0,	0,	0,					1,			},
		{0,			sqrt(f/g_rr),		0,	0,	0,	0,	0,	sqrt(f/g_rr),		0,			},
		{sqrt(g_rr),sqrt(f*g_rr)/g_hh,	0,	0,	0,	0,	0,	sqrt(f*g_rr)/g_hh,	sqrt(g_rr),	},
		{0,			-2,					0,	0,	0,	0,	0,	2,					0,			},
	}
	self.eigenvectorsInverse[i] = {
		{(g * A / f - K * sqrt(g / f)) / (2 * alpha), 0, g / (2 * f), 0, -.5 * sqrt(g / f)}, 
		{1 / alpha, 0, 0, 0, 0}, 
		{-(g * A) / (alpha * f), 0, -g / f, 1, 0}, 
		{0, 1, 0, 0, 0}, 
		{(g * A / f + K * sqrt(g / f)) / (2 * alpha), 0, g / (2 * f), 0, .5 * sqrt(g / f)}, 
	}
	-- note that because we have zero eigenvalues that the eigendecomposition cannot reconstruct the flux matrix
end

function ADM2DSim:addSourceToDerivCell(i)
	local alpha, g, A, D, K = unpack(self.qs[i])
	local f = self.calc_f(alpha)
	self.dq_dts[i][1] = self.dq_dts[i][1] - alpha * alpha * f * K / g
	self.dq_dts[i][2] = self.dq_dts[i][2] - 2 * alpha * K
	self.dq_dts[i][5] = self.dq_dts[i][5] + alpha * (A * D - K * K) / g
end

return ADM2DSim

