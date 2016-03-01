--[[
pseudo-cartesian
r = sqrt(x^2 + y^2 + z^2)
rs = schwarzschild radius
rr = rs/r

g_tt = -1 + rr
g_xx = (x^2/(1 - rr) + y^2 + z^2) / r^2
g_yy = (x^2 + y^2/(1 - rr) + z^2) / r^2
g_zz = (x^2 + y^2 + z^2/(1 - rr)) / r^2
g_xy = x y rs / (r^2 (r - rs))
g_xz = x z rs / (r^2 (r - rs))
g_yz = y z rs / (r^2 (r - rs))

...

	flux-less terms:

partial_t alpha = -alpha^2 f tr K
partial_t gamma_ij = -2 alpha K_ij

	flux terms:

partial_t A_k + alpha f gamma^ij partial_k K_ij = -alpha tr K (f + alpha f') A_k + 2 alpha f K^ij D_kij
partial_t D_kij + alpha partial_k K_ij = -alpha A_k K_ij
partial_t K_ij + alpha (gamma^km partial_k D_mij + 1/2 (delta^k_i partial_k A_j + delta^k_j partial_k A_i) + delta^k_i partial_k V_j + delta^k_j partial_k V_i) = alpha S_ij - alpha lambda^k_ij A_k + 2 alpha D_mij D_k^km
partial_t V_k = alpha P_k

	eigenfields:

lambda = 0:
	alpha
	gamma_ij
	V_i
	A_x'
	D_x'ij
	A_x - f D_x^m_m

lambda = +- alpha sqrt(gamma^xx):
	K_ix' +- sqrt(gamma_xx) (D_xix' + delta^x_i V_x' / gamma^xx) for x' != x

lambda = +- alpha sqrt(f gamma^xx):
	sqrt(f) tr(K) +- sqrt(gamma_xx) (A_x + 2 V^x / gamma^xx)

--]]
local class = require 'ext.class'
local Equation = require 'equation'
local mat33 = require 'mat33'
local symmath = require 'symmath'

local ADM3D = class(Equation)
ADM3D.name = 'ADM3D'

--[[
alpha,
gamma_xx, gamma_xy, gamma_xz, gamma_yy, gamma_yz, gamma_zz,
A_x, A_y, A_z,
D_xxx, D_xxy, D_xxz, D_xyy, D_xyz, D_xzz,
D_yxx, D_yxy, D_yxz, D_yyy, D_yyz, D_yzz,
D_zxx, D_zxy, D_zxz, D_zyy, D_zyz, D_zzz,
K_xx, K_xy, K_xz, K_yy, K_yz, K_zz,
V_x, V_y, V_z,
--]]
ADM3D.numStates = 37

function ADM3D:init(args, ...)
	local function makesym(field)
		local v = args[field]
		if not v then return end
		return symmath.clone(v):simplify() 
	end

	-- parameters that are variables of symbolic functions
	local x = assert(args.x)
	local y = args.y or symmath.var'y'
	local z = args.z or symmath.var'z'
	local vars = table{x,y,z}
	local symvars = {'xx', 'xy', 'xz', 'yy', 'yz', 'zz'}

	-- parameters that are symbolic functions -- based on coordinates 
	local exprs = table()
	exprs.alpha = assert(makesym'alpha')
	exprs.gamma = table{
		assert(makesym'gamma_xx'),	--xx
		makesym'gamma_xy' or symmath.Constant(0),	--xy
		makesym'gamma_xz' or symmath.Constant(0),	--xz
		makesym'gamma_yy' or symmath.Constant(1),	--yy
		makesym'gamma_yz' or symmath.Constant(0),	--yz
		makesym'gamma_zz' or symmath.Constant(1),	--zz
	}

	-- derived functions
	exprs.A = vars:map(function(var)
		return (exprs.alpha:diff(var) / exprs.alpha):simplify()
	end)
	exprs.K = table{
		assert(makesym'K_xx'),	--xx
		makesym'K_xy' or symmath.Constant(0),	--xy
		makesym'K_xz' or symmath.Constant(0),	--xz
		makesym'K_yy' or symmath.Constant(0),	--yy
		makesym'K_yz' or symmath.Constant(0),	--yz
		makesym'K_zz' or symmath.Constant(0),	--zz
	}
	
	local gammaUxx, gammaUxy, gammaUxz, gammaUyy, gammaUyz, gammaUzz = mat33.inv(unpack(exprs.gamma))
	exprs.gammaU = table{gammaUxx, gammaUxy, gammaUxz, gammaUyy, gammaUyz, gammaUzz}
	exprs.D = vars:map(function(x_k)
		return exprs.gamma:map(function(gamma_ij)
			return (gamma_ij:diff(x_k)/2):simplify()
		end)
	end)

	-- convert from symbolic functions to Lua functions
	local function buildCalc(expr, name)
		assert(type(expr) == 'table')
		if expr.isa and expr:isa(symmath.Expression) then
			return expr:compile{x,y,z}, name
		else
			return table.map(expr, buildCalc), name
		end
	end
	self.calc = exprs:map(buildCalc)
	
	local f_param = assert(args.f_param)

	local f = symmath.clone(assert(args.f)):simplify()
	self.calc.f = f:compile{f_param}

	local dalpha_f = f:diff(f_param):simplify()
	self.calc.dalpha_f = dalpha_f:compile{f_param}

--[=[ me trying to branch out on my own ...
	--[[
	from "Catalogue of Spacetimes"
	Schwarzschild Cartesian isotropic coordinates
	ds^2 = -((1 - rho_s/rho)/(1 + rho_s/rho))^2 dt^2 + (1 + rho_s/rho)^4 (dx^2 + dy^2 + dz^2)
	rho^2 = x^2 + y^2 + z^2
	r = rho (1 + rho_s / rho)^2

	what is r anyways?
	and is rho_s the Schwarzschild radius?
	I'm betting (as above in Schwarzschild isotropic spherical coordinates) rho_s is the Schwarzschild radius in the new isotropic radial coordinate
	 as a function of r_s in the non-isotropic radial coordinate (which is what describes the traditional Schwarzschild spherical dr)
	...which is given as ...
		r = rho (1 + rho_s / rho)^2
		rho = 1/4 (2r - rs +/- 2 sqrt(r (r - rs)) )
		rho_s = rs / 4
	... which means our xyz are no longer radially related to the 'r' coordinate of the original schwarzschild coordinates ...
	
	now for calculating extrinsic curvature...
	--]]	
	local rs = assert(tonumber(args.rs))	-- schwarzschild radius
	local rho_s = rs / 4	
	local rho = x^2 + y^2 + z^2
	local g_tt = -((1 - rho_s / rho)/(1 + rho_s/rho))^2
	local g_ii = (1 + rho_s/rho)^4
	local g_xx = g_ii
	local g_yy = g_ii
	local g_zz = g_ii
	local g_xy = symmath.Constant(0)
	local g_xz = symmath.Constant(0)
	local g_yz = symmath.Constant(0)
	local alpha = ((1 - rho_s/rho)/(1 + rho_s/rho))^2
--]=]	

	
	-- and for the graphs ...

	--[[ match 1D-3Var ...
	local get_state = index:bind(self.qs)
	local get_alpha = get_state(1)
	local get_gamma_xx = get_state(2)
	local get_A_x = get_state(8)
	local get_D_xxx = get_state(11)
	local get_K_xx = get_state(29)
	self.graphInfos = {
		{viewport={0/3, 0/3, 1/3, 1/3}, getter=get_alpha, name='alpha', color={1,0,1}},
		{viewport={0/3, 1/3, 1/3, 1/3}, getter=get_A_x, name='A_x', color={0,1,0}},
		{viewport={1/3, 0/3, 1/3, 1/3}, getter=get_gamma_xx, name='gamma_xx', color={.5,.5,1}},
		{viewport={1/3, 1/3, 1/3, 1/3}, getter=get_D_xxx, name='D_xxx', color={1,1,0}},
		{viewport={2/3, 0/3, 1/3, 1/3}, getter=get_K_xx, name='K_xx', color={0,1,1}},
		{viewport={2/3, 1/3, 1/3, 1/3}, getter=get_alpha * sqrt:compose(get_gamma_xx), name='volume', color={0,1,1}},
		{viewport={0/3, 2/3, 1/3, 1/3}, getter=log:compose(index:bind(self.eigenbasisErrors)), name='log eigenbasis error', color={1,0,0}, range={-30, 30}},
		{viewport={1/3, 2/3, 1/3, 1/3}, getter=log:compose(index:bind(self.fluxMatrixErrors)), name='log reconstuction error', color={1,0,0}, range={-30, 30}},
	}
	do return end
	--]]

	local i=0
	local j=0
	local function col() i=i+1 j=0 end
	self.graphInfos = table()

	local get_state = function(j) return function(self, i) return self.qs[i][j] end end
	local function add(args)
		local index = args.index
		for _,suffix in ipairs(args.suffix or {''}) do
			self.graphInfos:insert{
				viewport = {i,j},
				getter = get_state(index),
				name = args.name..suffix,
				color = {0,1,1},
			}
			index=index+1
			j=j+1
		end
	end

	add{name='alpha', index=1}
	do
		local state = function(self,i) return self.qs[i] end
		local alpha = state:index(1)
		local gamma_xx = state:index(2)
		local gamma_xy = state:index(3)
		local gamma_xz = state:index(4)
		local gamma_yy = state:index(5)
		local gamma_yz = state:index(6)
		local gamma_zz = state:index(7)
		self.graphInfos:insert{
			viewport = {i,j},
			-- TODO needs to subtract beta outer beta for trace -- or do det of upper and invert 
			getter = alpha * mat33.det(gamma_xx, gamma_xy, gamma_xz, gamma_yy, gamma_yz, gamma_zz),
			name = 'volume', 
			color = {0,1,1},
		}
	end
	j=j+1
	self.graphInfos:insert{
		viewport = {i,j},
		getter = function(self, i) return math.log(self.eigenbasisErrors[i]) end,
		name = 'log eigenbasis error', 
		color = {1,0,0}, 
		range = {-30, 30},
	}
	j=j+1
	self.graphInfos:insert{
		viewport = {i,j},
		getter = function(self, i) return math.log(self.fluxMatrixErrors[i]) end,
		name = 'log reconstruction error',
		color = {1,0,0},
		range = {-30, 30},
	}
	local suffix3 = {'x', 'y', 'z'}
	local suffix3x3sym = {'xx', 'xy', 'xz', 'yy', 'yz', 'zz'}
	col()
	add{name='A_', index=8, suffix=suffix3}
	add{name='V_', index=35, suffix=suffix3}
	col()
	add{name='gamma_', index=2, suffix=suffix3x3sym}
	col()
	add{name='D_x', index=11, suffix=suffix3x3sym}
	col()
	add{name='D_y', index=17, suffix=suffix3x3sym}
	col()
	add{name='D_z', index=23, suffix=suffix3x3sym}
	col()
	add{name='K_', index=29, suffix=suffix3x3sym}
	col()

	local xmax = 0
	local ymax = 0
	for _,info in ipairs(self.graphInfos) do
		xmax = max(xmax, info.viewport[1])
		ymax = max(ymax, info.viewport[2])
	end
	xmax = xmax + 1
	ymax = ymax + 1
	for _,info in ipairs(self.graphInfos) do
		local x,y = unpack(info.viewport)
		info.viewport = {x/xmax, y/ymax, 1/xmax, 1/ymax}
	end

-- [[ match the 1D 3-var layout:
	self.graphInfos = table{
		{viewport={0/3, 0/3, 1/3, 1/3}, getter=function(self,i) return self.qs[i][1] end, name='alpha', color={1,0,1}},
		{viewport={0/3, 1/3, 1/3, 1/3}, getter=function(self,i) return self.qs[i][8] end, name='A_x', color={0,1,0}},
		{viewport={1/3, 0/3, 1/3, 1/3}, getter=function(self,i) return self.qs[i][2] end, name='gamma_xx', color={.5,.5,1}},
		{viewport={1/3, 1/3, 1/3, 1/3}, getter=function(self,i) return self.qs[i][11] end, name='D_xxx', color={1,1,0}},
		{viewport={2/3, 0/3, 1/3, 1/3}, getter=function(self,i) return self.qs[i][29] end, name='K_xx', color={0,1,1}},
		{viewport={2/3, 1/3, 1/3, 1/3}, getter=function(self,i) return self.qs[i][1] * math.sqrt(self.qs[i][2]) end, name='volume', color={0,1,1}},
		{viewport={0/3, 2/3, 1/3, 1/3}, getter=function(self,i) return math.log(self.eigenbasisErrors[i]) end, name='log eigenbasis error', color={1,0,0}, range={-30, 30}},
		{viewport={1/3, 2/3, 1/3, 1/3}, getter=function(self,i) return math.log(self.fluxMatrixErrors[i]) end, name='log reconstuction error', color={1,0,0}, range={-30, 30}},
	}
--]]

	self.graphInfoForNames = self.graphInfos:map(function(info,i)
		return info, info.name
	end)
end

function ADM3D:fluxTransform(sim, i, v)
	local avgQ = {}
	for j=1,sim.numStates do 
		avgQ[j] = (sim.qs[i-1][j] + sim.qs[i][j]) / 2
	end

	-- ... is the incoming vector
	-- avgQ is the state used to make the eigenfield
	return {
		0,	--alpha
		0,0,0,0,0,0,	--gamma_ij
		0,0,0,		-- A_k
		0,0,0,0,0,0,	-- D_xij
		0,0,0,0,0,0,	-- D_yij
		0,0,0,0,0,0,	-- D_zij
		0,0,0,0,0,0,	-- K_ij
		0,0,0			-- V_k
	}
end

function ADM3D:eigenfields(sim, i, v)

	-- interface eigenfield varialbes
	local avgQ = {}
	for j=1,sim.numStates do 
		avgQ[j] = (sim.qs[i-1][j] + sim.qs[i][j]) / 2
		assert(type(avgQ[j])=='number')
	end
	local alpha = avgQ[1]
	local gamma_xx, gamma_xy, gamma_xz, gamma_yy, gamma_yz, gamma_zz = unpack(avgQ, 2, 7)
	local gamma = mat33.det(gamma_xx, gamma_xy, gamma_xz, gamma_yy, gamma_yz, gamma_zz)
	local A_x, A_y, A_z = unpack(avgQ, 8, 10)
	local D_xxx, D_xxy, D_xxz, D_xyy, D_xyz, D_xzz = unpack(avgQ, 11, 16)
	local D_yxx, D_yxy, D_yxz, D_yyy, D_yyz, D_yzz = unpack(avgQ, 17, 22)
	local D_zxx, D_zxy, D_zxz, D_zyy, D_zyz, D_zzz = unpack(avgQ, 23, 28)
	local K_xx, K_xy, K_xz, K_yy, K_yz, K_zz = unpack(avgQ, 29, 34)
	local V_x, V_y, V_z = unpack(avgQ, 35, 37)
	local gammaUxx, gammaUxy, gammaUxz, gammaUyy, gammaUyz, gammaUzz = mat33.inv(gamma_xx, gamma_xy, gamma_xz, gamma_yy, gamma_yz, gamma_zz)
	local f = self.calc.f(alpha)

	-- cell variables
	-- what if, for the ADM equations, there is no distinction?
	-- they're used for Roe's scheme for computing deltas in eigenbasis coordinates by which to scale coordinates coinciding with the lambdas ...
	-- what about creating them solely from 'v' rather than using the average whatsoever?
	-- this would mean ensuring the inputs to the eigenfields() functions were always the state variables themselves (not differences or averages)
	-- 	and deferring differences or averages til after eigenfields() is called (assuming it is a linear function)
	-- this also has an issue with eigenfieldsInverse(), which is called on a flux vector, i.e. at cell interface, which would probably need the average of cells for that input

	return {
		((((-(2 * gammaUxz * v[37])) - (gammaUxx * v[8])) + (math.sqrt(f) * (gammaUxx ^ (3 / 2)) * v[29]) + (2 * math.sqrt(f) * gammaUxy * v[30] * math.sqrt(gammaUxx)) + (2 * math.sqrt(f) * gammaUxz * v[31] * math.sqrt(gammaUxx)) + (math.sqrt(f) * gammaUyy * v[32] * math.sqrt(gammaUxx)) + (2 * math.sqrt(f) * gammaUyz * v[33] * math.sqrt(gammaUxx)) + (((math.sqrt(f) * gammaUzz * v[34] * math.sqrt(gammaUxx)) - (2 * gammaUxx * v[35])) - (2 * gammaUxy * v[36]))) / math.sqrt(gammaUxx)),
		(((-(gammaUxx * v[12])) + ((v[30] * math.sqrt(gammaUxx)) - v[36])) / math.sqrt(gammaUxx)),
		(((-(gammaUxx * v[13])) + ((v[31] * math.sqrt(gammaUxx)) - v[37])) / math.sqrt(gammaUxx)),
		((-(math.sqrt(gammaUxx) * v[14])) + v[32]),
		((-(math.sqrt(gammaUxx) * v[15])) + v[33]),
		((-(math.sqrt(gammaUxx) * v[16])) + v[34]),
		v[1],
		v[2],
		v[3],
		v[4],
		v[5],
		v[6],
		v[7],
		v[9],
		v[10],
		v[17],
		v[18],
		v[19],
		v[20],
		v[21],
		v[22],
		v[23],
		v[24],
		v[25],
		v[26],
		v[27],
		v[28],
		v[35],
		v[36],
		v[37],
		((((((v[8] - (f * gammaUxx * v[11])) - (2 * f * gammaUxy * v[12])) - (2 * f * gammaUxz * v[13])) - (f * gammaUyy * v[14])) - (2 * f * gammaUyz * v[15])) - (f * gammaUzz * v[16])),
		(((gammaUxx * v[12]) + (v[30] * math.sqrt(gammaUxx)) + v[36]) / math.sqrt(gammaUxx)),
		(((gammaUxx * v[13]) + (v[31] * math.sqrt(gammaUxx)) + v[37]) / math.sqrt(gammaUxx)),
		((math.sqrt(gammaUxx) * v[14]) + v[32]),
		((math.sqrt(gammaUxx) * v[15]) + v[33]),
		((math.sqrt(gammaUxx) * v[16]) + v[34]),
		(((2 * gammaUxz * v[37]) + (gammaUxx * v[8]) + (math.sqrt(f) * (gammaUxx ^ (3 / 2)) * v[29]) + (2 * math.sqrt(f) * gammaUxy * v[30] * math.sqrt(gammaUxx)) + (2 * math.sqrt(f) * gammaUxz * v[31] * math.sqrt(gammaUxx)) + (math.sqrt(f) * gammaUyy * v[32] * math.sqrt(gammaUxx)) + (2 * math.sqrt(f) * gammaUyz * v[33] * math.sqrt(gammaUxx)) + (math.sqrt(f) * gammaUzz * v[34] * math.sqrt(gammaUxx)) + (2 * gammaUxx * v[35]) + (2 * gammaUxy * v[36])) / math.sqrt(gammaUxx))
	}
end

function ADM3D:eigenfieldsInverse(sim, i, v)
	
	-- interface eigenfield varialbes
	local avgQ = {}
	for j=1,sim.numStates do 
		avgQ[j] = (sim.qs[i-1][j] + sim.qs[i][j]) / 2
	end
	local alpha = avgQ[1]
	local gamma_xx, gamma_xy, gamma_xz, gamma_yy, gamma_yz, gamma_zz = unpack(avgQ, 2, 7)
	local gamma = mat33.det(gamma_xx, gamma_xy, gamma_xz, gamma_yy, gamma_yz, gamma_zz)
	local A_x, A_y, A_z = unpack(avgQ, 8, 10)
	local D_xxx, D_xxy, D_xxz, D_xyy, D_xyz, D_xzz = unpack(avgQ, 11, 16)
	local D_yxx, D_yxy, D_yxz, D_yyy, D_yyz, D_yzz = unpack(avgQ, 17, 22)
	local D_zxx, D_zxy, D_zxz, D_zyy, D_zyz, D_zzz = unpack(avgQ, 23, 28)
	local K_xx, K_xy, K_xz, K_yy, K_yz, K_zz = unpack(avgQ, 29, 34)
	local V_x, V_y, V_z = unpack(avgQ, 35, 37)
	local gammaUxx, gammaUxy, gammaUxz, gammaUyy, gammaUyz, gammaUzz = mat33.inv(gamma_xx, gamma_xy, gamma_xz, gamma_yy, gamma_yz, gamma_zz)
	local f = self.calc.f(alpha)

	return {
		v[7],
		v[8],
		v[9],
		v[10],
		v[11],
		v[12],
		v[13],
		(((-v[37]) + (4 * gammaUxz * v[30] * (1 / math.sqrt(gammaUxx))) + (4 * gammaUxy * v[29] * (1 / math.sqrt(gammaUxx))) + (4 * v[28] * math.sqrt(gammaUxx)) + v[1]) / (-(2 * math.sqrt(gammaUxx)))),
		v[14],
		v[15],
		(((-v[37]) + (gammaUzz * v[36] * f) + (2 * gammaUyz * v[35] * f) + (gammaUyy * v[34] * f) + (2 * gammaUxz * v[33] * f) + (2 * gammaUxy * v[32] * f) + (2 * v[31] * math.sqrt(gammaUxx)) + ((4 * gammaUxz * v[30] * (1 / math.sqrt(gammaUxx))) - (4 * gammaUxz * f * v[30] * (1 / math.sqrt(gammaUxx)))) + ((4 * gammaUxy * v[29] * (1 / math.sqrt(gammaUxx))) - (4 * gammaUxy * f * v[29] * (1 / math.sqrt(gammaUxx)))) + ((((((4 * v[28] * math.sqrt(gammaUxx)) - (gammaUzz * v[6] * f)) - (2 * gammaUyz * v[5] * f)) - (gammaUyy * v[4] * f)) - (2 * gammaUxz * v[3] * f)) - (2 * gammaUxy * v[2] * f)) + v[1]) / (-(2 * (gammaUxx ^ (3 / 2)) * f))),
		(((-v[32]) + (2 * v[29] * (1 / math.sqrt(gammaUxx))) + v[2]) / (-(2 * math.sqrt(gammaUxx)))),
		(((-v[33]) + (2 * v[30] * (1 / math.sqrt(gammaUxx))) + v[3]) / (-(2 * math.sqrt(gammaUxx)))),
		(((-v[34]) + v[4]) / (-(2 * math.sqrt(gammaUxx)))),
		(((-v[35]) + v[5]) / (-(2 * math.sqrt(gammaUxx)))),
		(((-v[36]) + v[6]) / (-(2 * math.sqrt(gammaUxx)))),
		v[16],
		v[17],
		v[18],
		v[19],
		v[20],
		v[21],
		v[22],
		v[23],
		v[24],
		v[25],
		v[26],
		v[27],
		((((((((((((v[37] - (gammaUzz * v[36] * math.sqrt(f))) - (2 * gammaUyz * v[35] * math.sqrt(f))) - (gammaUyy * v[34] * math.sqrt(f))) - (2 * gammaUxz * v[33] * math.sqrt(f))) - (2 * gammaUxy * v[32] * math.sqrt(f))) - (gammaUzz * v[6] * math.sqrt(f))) - (2 * gammaUyz * v[5] * math.sqrt(f))) - (gammaUyy * v[4] * math.sqrt(f))) - (2 * gammaUxz * v[3] * math.sqrt(f))) - (2 * gammaUxy * v[2] * math.sqrt(f))) + v[1]) / (2 * math.sqrt(f) * gammaUxx)),
		((v[32] + v[2]) / 2),
		((v[33] + v[3]) / 2),
		((v[34] + v[4]) / 2),
		((v[35] + v[5]) / 2),
		((v[36] + v[6]) / 2),
		v[28],
		v[29],
		v[30]
	}
end


function ADM3D:initCell(sim,i)
	local x = sim.xs[i]
	local y = 0
	local z = 0
	local xs = table{x,y,z}
	local alpha = self.calc.alpha(x,y,z)
	local A = self.calc.A:map(function(A_i) return A_i(x,y,z) end)
	local gamma = self.calc.gamma:map(function(gamma_ij) return gamma_ij(x,y,z) end)
	local D = self.calc.D:map(function(D_i) return D_i:map(function(D_ijk) return D_ijk(x,y,z) end) end)
	local gammaU = self.calc.gammaU:map(function(gammaUij) return gammaUij(x,y,z) end)

	local function sym3x3(m,i,j)
		local m_xx, m_xy, m_xz, m_yy, m_yz, m_zz = unpack(m)
		return ({
			{m_xx, m_xy, m_xz},
			{m_xy, m_yy, m_yz},
			{m_xz, m_yz, m_zz},
		})[i][j]
	end
	local V = range(3):map(function(i)
		local s = 0
		for j=1,3 do
			for k=1,3 do
				local D_ijk = sym3x3(D[i],j,k)
				local D_kji = sym3x3(D[k],j,i)
				local gammaUjk = sym3x3(gammaU,j,k)
				local dg = (D_ijk - D_kji) * gammaUjk
				s = s + dg
			end
		end
		return s
	end)

	local K = {}
	for i=1,6 do
		K[i] = self.calc.K[i](x,y,z)
	end

	return {
		alpha,
		gamma[1], gamma[2], gamma[3], gamma[4], gamma[5], gamma[6],
		A[1], A[2], A[3],
		D[1][1], D[1][2], D[1][3], D[1][4], D[1][5], D[1][6],
		D[2][1], D[2][2], D[2][3], D[2][4], D[2][5], D[2][6],
		D[3][1], D[3][2], D[3][3], D[3][4], D[3][5], D[3][6],
		K[1], K[2], K[3], K[4], K[5], K[6],
		V[1], V[2], V[3],
	}
end

function ADM3D:calcInterfaceEigenBasis(sim,i,qL,qR)
	local avgQ = {}
	for j=1,self.numStates do
		avgQ[j] = (qL[j] + qR[j]) / 2
	end

	local alpha = avgQ[1]
	local gamma_xx, gamma_xy, gamma_xz, gamma_yy, gamma_yz, gamma_zz = unpack(avgQ, 2, 7)
	local A_x, A_y, A_z = unpack(avgQ, 8, 10)
	local D_xxx, D_xxy, D_xxz, D_xyy, D_xyz, D_xzz = unpack(avgQ, 11, 16)
	local D_yxx, D_yxy, D_yxz, D_yyy, D_yyz, D_yzz = unpack(avgQ, 17, 22)
	local D_zxx, D_zxy, D_zxz, D_zyy, D_zyz, D_zzz = unpack(avgQ, 23, 28)
	local K_xx, K_xy, K_xz, K_yy, K_yz, K_zz = unpack(avgQ, 29, 34)
	local V_x, V_y, V_z = unpack(avgQ, 35, 37)
	local gammaUxx, gammaUxy, gammaUxz, gammaUyy, gammaUyz, gammaUzz = mat33.inv(gamma_xx, gamma_xy, gamma_xz, gamma_yy, gamma_yz, gamma_zz)
	local f = self.calc.f(alpha)

	local lambdaLight = alpha * sqrt(gammaUxx)
	local lambdaGauge = lambdaLight * sqrt(f)
	sim.eigenvalues[i] = {
		-- gauge field
		-lambdaGauge,
		-- half of 10 along light cones ...
		-lambdaLight,
		-lambdaLight,
		-lambdaLight,
		-lambdaLight,
		-lambdaLight,
		-- 25 zeroes ...
		0,0,0,0,0,
		0,0,0,0,0,
		0,0,0,0,0,
		0,0,0,0,0,
		0,0,0,0,0,
		-- half of 10 along light cones ...
		lambdaLight,
		lambdaLight,
		lambdaLight,
		lambdaLight,
		lambdaLight,
		-- gauge field
		lambdaGauge
	}
end

function ADM3D:sourceTerm(sim, qs)
	local source = sim:newState()
	for i=1,sim.gridsize do
		local alpha = qs[i][1]
		local gamma_xx, gamma_xy, gamma_xz, gamma_yy, gamma_yz, gamma_zz = unpack(qs[i], 2, 7)
		local A_x, A_y, A_z = unpack(qs[i], 8, 10)
		local D_xxx, D_xxy, D_xxz, D_xyy, D_xyz, D_xzz = unpack(qs[i], 11, 16)
		local D_yxx, D_yxy, D_yxz, D_yyy, D_yyz, D_yzz = unpack(qs[i], 17, 22)
		local D_zxx, D_zxy, D_zxz, D_zyy, D_zyz, D_zzz = unpack(qs[i], 23, 28)
		local K_xx, K_xy, K_xz, K_yy, K_yz, K_zz = unpack(qs[i], 29, 34)
		local V_x, V_y, V_z = unpack(qs[i], 35, 37)
		local gammaUxx, gammaUxy, gammaUxz, gammaUyy, gammaUyz, gammaUzz = mat33.inv(gamma_xx, gamma_xy, gamma_xz, gamma_yy, gamma_yz, gamma_zz)
		local f = self.calc.f(alpha)

		-- constraint variables for K_ij and V_k

-- source terms
local KUL = {
{gammaUxx * K_xx + gammaUxy * K_xy + gammaUxz * K_xz,
gammaUxx * K_xy + gammaUxy * K_yy + gammaUxz * K_yz,
gammaUxx * K_xz + gammaUxy * K_yz + gammaUxz * K_zz,
},{gammaUxy * K_xx + gammaUyy * K_xy + gammaUyz * K_xz,
gammaUxy * K_xy + gammaUyy * K_yy + gammaUyz * K_yz,
gammaUxy * K_xz + gammaUyy * K_yz + gammaUyz * K_zz,
},{gammaUxz * K_xx + gammaUyz * K_xy + gammaUzz * K_xz,
gammaUxz * K_xy + gammaUyz * K_yy + gammaUzz * K_yz,
gammaUxz * K_xz + gammaUyz * K_yz + gammaUzz * K_zz,
},}
local trK = KUL[1][1] + KUL[2][2] + KUL[3][3]
local KSqSymLL = {
K_xx * KUL[1][1] + K_xy * KUL[2][1] + K_xz * KUL[3][1],
K_xx * KUL[1][2] + K_xy * KUL[2][2] + K_xz * KUL[3][2],
K_xx * KUL[1][3] + K_xy * KUL[2][3] + K_xz * KUL[3][3],
K_xy * KUL[1][2] + K_yy * KUL[2][2] + K_yz * KUL[3][2],
K_xy * KUL[1][3] + K_yy * KUL[2][3] + K_yz * KUL[3][3],
K_xz * KUL[1][3] + K_yz * KUL[2][3] + K_zz * KUL[3][3],
}
local DLUL = {
{{D_xxx * gammaUxx + D_xxy * gammaUxy + D_xxz * gammaUxz,
D_xxy * gammaUxx + D_xyy * gammaUxy + D_xyz * gammaUxz,
D_xxz * gammaUxx + D_xyz * gammaUxy + D_xzz * gammaUxz,
},{D_xxx * gammaUxy + D_xxy * gammaUyy + D_xxz * gammaUyz,
D_xxy * gammaUxy + D_xyy * gammaUyy + D_xyz * gammaUyz,
D_xxz * gammaUxy + D_xyz * gammaUyy + D_xzz * gammaUyz,
},{D_xxx * gammaUxz + D_xxy * gammaUyz + D_xxz * gammaUzz,
D_xxy * gammaUxz + D_xyy * gammaUyz + D_xyz * gammaUzz,
D_xxz * gammaUxz + D_xyz * gammaUyz + D_xzz * gammaUzz,
},},{{D_yxx * gammaUxx + D_yxy * gammaUxy + D_yxz * gammaUxz,
D_yxy * gammaUxx + D_yyy * gammaUxy + D_yyz * gammaUxz,
D_yxz * gammaUxx + D_yyz * gammaUxy + D_yzz * gammaUxz,
},{D_yxx * gammaUxy + D_yxy * gammaUyy + D_yxz * gammaUyz,
D_yxy * gammaUxy + D_yyy * gammaUyy + D_yyz * gammaUyz,
D_yxz * gammaUxy + D_yyz * gammaUyy + D_yzz * gammaUyz,
},{D_yxx * gammaUxz + D_yxy * gammaUyz + D_yxz * gammaUzz,
D_yxy * gammaUxz + D_yyy * gammaUyz + D_yyz * gammaUzz,
D_yxz * gammaUxz + D_yyz * gammaUyz + D_yzz * gammaUzz,
},},{{D_zxx * gammaUxx + D_zxy * gammaUxy + D_zxz * gammaUxz,
D_zxy * gammaUxx + D_zyy * gammaUxy + D_zyz * gammaUxz,
D_zxz * gammaUxx + D_zyz * gammaUxy + D_zzz * gammaUxz,
},{D_zxx * gammaUxy + D_zxy * gammaUyy + D_zxz * gammaUyz,
D_zxy * gammaUxy + D_zyy * gammaUyy + D_zyz * gammaUyz,
D_zxz * gammaUxy + D_zyz * gammaUyy + D_zzz * gammaUyz,
},{D_zxx * gammaUxz + D_zxy * gammaUyz + D_zxz * gammaUzz,
D_zxy * gammaUxz + D_zyy * gammaUyz + D_zyz * gammaUzz,
D_zxz * gammaUxz + D_zyz * gammaUyz + D_zzz * gammaUzz,
},},}
local D1L = {
DLUL[1][1][1] + DLUL[1][2][2] + DLUL[1][3][3],
DLUL[2][1][1] + DLUL[2][2][2] + DLUL[2][3][3],
DLUL[3][1][1] + DLUL[3][2][2] + DLUL[3][3][3],
}
local D3L = {
DLUL[1][1][1] + DLUL[2][2][1] + DLUL[3][3][1],
DLUL[1][1][2] + DLUL[2][2][2] + DLUL[3][3][2],
DLUL[1][1][3] + DLUL[2][2][3] + DLUL[3][3][3],
}
local DUUL = {
{{DLUL[1][1][1] * gammaUxx + DLUL[2][1][1] * gammaUxy + DLUL[3][1][1] * gammaUxz,
DLUL[1][1][2] * gammaUxx + DLUL[2][1][2] * gammaUxy + DLUL[3][1][2] * gammaUxz,
DLUL[1][1][3] * gammaUxx + DLUL[2][1][3] * gammaUxy + DLUL[3][1][3] * gammaUxz,
},{DLUL[1][2][1] * gammaUxx + DLUL[2][2][1] * gammaUxy + DLUL[3][2][1] * gammaUxz,
DLUL[1][2][2] * gammaUxx + DLUL[2][2][2] * gammaUxy + DLUL[3][2][2] * gammaUxz,
DLUL[1][2][3] * gammaUxx + DLUL[2][2][3] * gammaUxy + DLUL[3][2][3] * gammaUxz,
},{DLUL[1][3][1] * gammaUxx + DLUL[2][3][1] * gammaUxy + DLUL[3][3][1] * gammaUxz,
DLUL[1][3][2] * gammaUxx + DLUL[2][3][2] * gammaUxy + DLUL[3][3][2] * gammaUxz,
DLUL[1][3][3] * gammaUxx + DLUL[2][3][3] * gammaUxy + DLUL[3][3][3] * gammaUxz,
},},{{DLUL[1][1][1] * gammaUxy + DLUL[2][1][1] * gammaUyy + DLUL[3][1][1] * gammaUyz,
DLUL[1][1][2] * gammaUxy + DLUL[2][1][2] * gammaUyy + DLUL[3][1][2] * gammaUyz,
DLUL[1][1][3] * gammaUxy + DLUL[2][1][3] * gammaUyy + DLUL[3][1][3] * gammaUyz,
},{DLUL[1][2][1] * gammaUxy + DLUL[2][2][1] * gammaUyy + DLUL[3][2][1] * gammaUyz,
DLUL[1][2][2] * gammaUxy + DLUL[2][2][2] * gammaUyy + DLUL[3][2][2] * gammaUyz,
DLUL[1][2][3] * gammaUxy + DLUL[2][2][3] * gammaUyy + DLUL[3][2][3] * gammaUyz,
},{DLUL[1][3][1] * gammaUxy + DLUL[2][3][1] * gammaUyy + DLUL[3][3][1] * gammaUyz,
DLUL[1][3][2] * gammaUxy + DLUL[2][3][2] * gammaUyy + DLUL[3][3][2] * gammaUyz,
DLUL[1][3][3] * gammaUxy + DLUL[2][3][3] * gammaUyy + DLUL[3][3][3] * gammaUyz,
},},{{DLUL[1][1][1] * gammaUxz + DLUL[2][1][1] * gammaUyz + DLUL[3][1][1] * gammaUzz,
DLUL[1][1][2] * gammaUxz + DLUL[2][1][2] * gammaUyz + DLUL[3][1][2] * gammaUzz,
DLUL[1][1][3] * gammaUxz + DLUL[2][1][3] * gammaUyz + DLUL[3][1][3] * gammaUzz,
},{DLUL[1][2][1] * gammaUxz + DLUL[2][2][1] * gammaUyz + DLUL[3][2][1] * gammaUzz,
DLUL[1][2][2] * gammaUxz + DLUL[2][2][2] * gammaUyz + DLUL[3][2][2] * gammaUzz,
DLUL[1][2][3] * gammaUxz + DLUL[2][2][3] * gammaUyz + DLUL[3][2][3] * gammaUzz,
},{DLUL[1][3][1] * gammaUxz + DLUL[2][3][1] * gammaUyz + DLUL[3][3][1] * gammaUzz,
DLUL[1][3][2] * gammaUxz + DLUL[2][3][2] * gammaUyz + DLUL[3][3][2] * gammaUzz,
DLUL[1][3][3] * gammaUxz + DLUL[2][3][3] * gammaUyz + DLUL[3][3][3] * gammaUzz,
},},}
local D12SymLL = {
D_xxx * DUUL[1][1][1] + D_xxy * DUUL[1][2][1] + D_xxz * DUUL[1][3][1] + D_yxx * DUUL[2][1][1] + D_yxy * DUUL[2][2][1] + D_yxz * DUUL[2][3][1] + D_zxx * DUUL[3][1][1] + D_zxy * DUUL[3][2][1] + D_zxz * DUUL[3][3][1],
D_xxy * DUUL[1][1][1] + D_xyy * DUUL[1][2][1] + D_xyz * DUUL[1][3][1] + D_yxy * DUUL[2][1][1] + D_yyy * DUUL[2][2][1] + D_yyz * DUUL[2][3][1] + D_zxy * DUUL[3][1][1] + D_zyy * DUUL[3][2][1] + D_zyz * DUUL[3][3][1],
D_xxz * DUUL[1][1][1] + D_xyz * DUUL[1][2][1] + D_xzz * DUUL[1][3][1] + D_yxz * DUUL[2][1][1] + D_yyz * DUUL[2][2][1] + D_yzz * DUUL[2][3][1] + D_zxz * DUUL[3][1][1] + D_zyz * DUUL[3][2][1] + D_zzz * DUUL[3][3][1],
D_xxy * DUUL[1][1][2] + D_xyy * DUUL[1][2][2] + D_xyz * DUUL[1][3][2] + D_yxy * DUUL[2][1][2] + D_yyy * DUUL[2][2][2] + D_yyz * DUUL[2][3][2] + D_zxy * DUUL[3][1][2] + D_zyy * DUUL[3][2][2] + D_zyz * DUUL[3][3][2],
D_xxz * DUUL[1][1][2] + D_xyz * DUUL[1][2][2] + D_xzz * DUUL[1][3][2] + D_yxz * DUUL[2][1][2] + D_yyz * DUUL[2][2][2] + D_yzz * DUUL[2][3][2] + D_zxz * DUUL[3][1][2] + D_zyz * DUUL[3][2][2] + D_zzz * DUUL[3][3][2],
D_xxz * DUUL[1][1][3] + D_xyz * DUUL[1][2][3] + D_xzz * DUUL[1][3][3] + D_yxz * DUUL[2][1][3] + D_yyz * DUUL[2][2][3] + D_yzz * DUUL[2][3][3] + D_zxz * DUUL[3][1][3] + D_zyz * DUUL[3][2][3] + D_zzz * DUUL[3][3][3],
}
local GammaLSymLL = {
{D_xxx,
D_yxx,
D_zxx,
((2 * D_yxy) - D_xyy),
(D_zxy + (D_yxz - D_xyz)),
((2 * D_zxz) - D_xzz),
},{((2 * D_xxy) - D_yxx),
D_xyy,
(D_zxy + (D_xyz - D_yxz)),
D_yyy,
D_zyy,
((2 * D_zyz) - D_yzz),
},{((2 * D_xxz) - D_zxx),
(D_yxz + (D_xyz - D_zxy)),
D_xzz,
((2 * D_yyz) - D_zyy),
D_yzz,
D_zzz,
},}
local GammaUSymLL = {
{gammaUxx * GammaLSymLL[1][1] + gammaUxy * GammaLSymLL[2][1] + gammaUxz * GammaLSymLL[3][1],
gammaUxx * GammaLSymLL[1][2] + gammaUxy * GammaLSymLL[2][2] + gammaUxz * GammaLSymLL[3][2],
gammaUxx * GammaLSymLL[1][3] + gammaUxy * GammaLSymLL[2][3] + gammaUxz * GammaLSymLL[3][3],
gammaUxx * GammaLSymLL[1][4] + gammaUxy * GammaLSymLL[2][4] + gammaUxz * GammaLSymLL[3][4],
gammaUxx * GammaLSymLL[1][5] + gammaUxy * GammaLSymLL[2][5] + gammaUxz * GammaLSymLL[3][5],
gammaUxx * GammaLSymLL[1][6] + gammaUxy * GammaLSymLL[2][6] + gammaUxz * GammaLSymLL[3][6],
},{gammaUxy * GammaLSymLL[1][1] + gammaUyy * GammaLSymLL[2][1] + gammaUyz * GammaLSymLL[3][1],
gammaUxy * GammaLSymLL[1][2] + gammaUyy * GammaLSymLL[2][2] + gammaUyz * GammaLSymLL[3][2],
gammaUxy * GammaLSymLL[1][3] + gammaUyy * GammaLSymLL[2][3] + gammaUyz * GammaLSymLL[3][3],
gammaUxy * GammaLSymLL[1][4] + gammaUyy * GammaLSymLL[2][4] + gammaUyz * GammaLSymLL[3][4],
gammaUxy * GammaLSymLL[1][5] + gammaUyy * GammaLSymLL[2][5] + gammaUyz * GammaLSymLL[3][5],
gammaUxy * GammaLSymLL[1][6] + gammaUyy * GammaLSymLL[2][6] + gammaUyz * GammaLSymLL[3][6],
},{gammaUxz * GammaLSymLL[1][1] + gammaUyz * GammaLSymLL[2][1] + gammaUzz * GammaLSymLL[3][1],
gammaUxz * GammaLSymLL[1][2] + gammaUyz * GammaLSymLL[2][2] + gammaUzz * GammaLSymLL[3][2],
gammaUxz * GammaLSymLL[1][3] + gammaUyz * GammaLSymLL[2][3] + gammaUzz * GammaLSymLL[3][3],
gammaUxz * GammaLSymLL[1][4] + gammaUyz * GammaLSymLL[2][4] + gammaUzz * GammaLSymLL[3][4],
gammaUxz * GammaLSymLL[1][5] + gammaUyz * GammaLSymLL[2][5] + gammaUzz * GammaLSymLL[3][5],
gammaUxz * GammaLSymLL[1][6] + gammaUyz * GammaLSymLL[2][6] + gammaUzz * GammaLSymLL[3][6],
},}
local Gamma3L = {
GammaUSymLL[1][1] + GammaUSymLL[2][2] + GammaUSymLL[3][3],
GammaUSymLL[1][2] + GammaUSymLL[2][4] + GammaUSymLL[3][5],
GammaUSymLL[1][3] + GammaUSymLL[2][5] + GammaUSymLL[3][6],
}
local Gamma31SymLL = {
Gamma3L[1] * GammaUSymLL[1][1] + Gamma3L[2] * GammaUSymLL[2][1] + Gamma3L[3] * GammaUSymLL[3][1],
Gamma3L[1] * GammaUSymLL[1][2] + Gamma3L[2] * GammaUSymLL[2][2] + Gamma3L[3] * GammaUSymLL[3][2],
Gamma3L[1] * GammaUSymLL[1][3] + Gamma3L[2] * GammaUSymLL[2][3] + Gamma3L[3] * GammaUSymLL[3][3],
Gamma3L[1] * GammaUSymLL[1][4] + Gamma3L[2] * GammaUSymLL[2][4] + Gamma3L[3] * GammaUSymLL[3][4],
Gamma3L[1] * GammaUSymLL[1][5] + Gamma3L[2] * GammaUSymLL[2][5] + Gamma3L[3] * GammaUSymLL[3][5],
Gamma3L[1] * GammaUSymLL[1][6] + Gamma3L[2] * GammaUSymLL[2][6] + Gamma3L[3] * GammaUSymLL[3][6],
}
local GammaLUL = {
{{gammaUxx * GammaLSymLL[1][1] + gammaUxy * GammaLSymLL[1][2] + gammaUxz * GammaLSymLL[1][3],
gammaUxx * GammaLSymLL[1][2] + gammaUxy * GammaLSymLL[1][4] + gammaUxz * GammaLSymLL[1][5],
gammaUxx * GammaLSymLL[1][3] + gammaUxy * GammaLSymLL[1][5] + gammaUxz * GammaLSymLL[1][6],
},{gammaUxy * GammaLSymLL[1][1] + gammaUyy * GammaLSymLL[1][2] + gammaUyz * GammaLSymLL[1][3],
gammaUxy * GammaLSymLL[1][2] + gammaUyy * GammaLSymLL[1][4] + gammaUyz * GammaLSymLL[1][5],
gammaUxy * GammaLSymLL[1][3] + gammaUyy * GammaLSymLL[1][5] + gammaUyz * GammaLSymLL[1][6],
},{gammaUxz * GammaLSymLL[1][1] + gammaUyz * GammaLSymLL[1][2] + gammaUzz * GammaLSymLL[1][3],
gammaUxz * GammaLSymLL[1][2] + gammaUyz * GammaLSymLL[1][4] + gammaUzz * GammaLSymLL[1][5],
gammaUxz * GammaLSymLL[1][3] + gammaUyz * GammaLSymLL[1][5] + gammaUzz * GammaLSymLL[1][6],
},},{{gammaUxx * GammaLSymLL[2][1] + gammaUxy * GammaLSymLL[2][2] + gammaUxz * GammaLSymLL[2][3],
gammaUxx * GammaLSymLL[2][2] + gammaUxy * GammaLSymLL[2][4] + gammaUxz * GammaLSymLL[2][5],
gammaUxx * GammaLSymLL[2][3] + gammaUxy * GammaLSymLL[2][5] + gammaUxz * GammaLSymLL[2][6],
},{gammaUxy * GammaLSymLL[2][1] + gammaUyy * GammaLSymLL[2][2] + gammaUyz * GammaLSymLL[2][3],
gammaUxy * GammaLSymLL[2][2] + gammaUyy * GammaLSymLL[2][4] + gammaUyz * GammaLSymLL[2][5],
gammaUxy * GammaLSymLL[2][3] + gammaUyy * GammaLSymLL[2][5] + gammaUyz * GammaLSymLL[2][6],
},{gammaUxz * GammaLSymLL[2][1] + gammaUyz * GammaLSymLL[2][2] + gammaUzz * GammaLSymLL[2][3],
gammaUxz * GammaLSymLL[2][2] + gammaUyz * GammaLSymLL[2][4] + gammaUzz * GammaLSymLL[2][5],
gammaUxz * GammaLSymLL[2][3] + gammaUyz * GammaLSymLL[2][5] + gammaUzz * GammaLSymLL[2][6],
},},{{gammaUxx * GammaLSymLL[3][1] + gammaUxy * GammaLSymLL[3][2] + gammaUxz * GammaLSymLL[3][3],
gammaUxx * GammaLSymLL[3][2] + gammaUxy * GammaLSymLL[3][4] + gammaUxz * GammaLSymLL[3][5],
gammaUxx * GammaLSymLL[3][3] + gammaUxy * GammaLSymLL[3][5] + gammaUxz * GammaLSymLL[3][6],
},{gammaUxy * GammaLSymLL[3][1] + gammaUyy * GammaLSymLL[3][2] + gammaUyz * GammaLSymLL[3][3],
gammaUxy * GammaLSymLL[3][2] + gammaUyy * GammaLSymLL[3][4] + gammaUyz * GammaLSymLL[3][5],
gammaUxy * GammaLSymLL[3][3] + gammaUyy * GammaLSymLL[3][5] + gammaUyz * GammaLSymLL[3][6],
},{gammaUxz * GammaLSymLL[3][1] + gammaUyz * GammaLSymLL[3][2] + gammaUzz * GammaLSymLL[3][3],
gammaUxz * GammaLSymLL[3][2] + gammaUyz * GammaLSymLL[3][4] + gammaUzz * GammaLSymLL[3][5],
gammaUxz * GammaLSymLL[3][3] + gammaUyz * GammaLSymLL[3][5] + gammaUzz * GammaLSymLL[3][6],
},},}
local GammaLSymUU = {
{gammaUxx * GammaLUL[1][1][1] + gammaUxy * GammaLUL[1][1][2] + gammaUxz * GammaLUL[1][1][3],
gammaUxy * GammaLUL[1][1][1] + gammaUyy * GammaLUL[1][1][2] + gammaUyz * GammaLUL[1][1][3],
gammaUxz * GammaLUL[1][1][1] + gammaUyz * GammaLUL[1][1][2] + gammaUzz * GammaLUL[1][1][3],
gammaUxy * GammaLUL[1][2][1] + gammaUyy * GammaLUL[1][2][2] + gammaUyz * GammaLUL[1][2][3],
gammaUxz * GammaLUL[1][2][1] + gammaUyz * GammaLUL[1][2][2] + gammaUzz * GammaLUL[1][2][3],
gammaUxz * GammaLUL[1][3][1] + gammaUyz * GammaLUL[1][3][2] + gammaUzz * GammaLUL[1][3][3],
},{gammaUxx * GammaLUL[2][1][1] + gammaUxy * GammaLUL[2][1][2] + gammaUxz * GammaLUL[2][1][3],
gammaUxy * GammaLUL[2][1][1] + gammaUyy * GammaLUL[2][1][2] + gammaUyz * GammaLUL[2][1][3],
gammaUxz * GammaLUL[2][1][1] + gammaUyz * GammaLUL[2][1][2] + gammaUzz * GammaLUL[2][1][3],
gammaUxy * GammaLUL[2][2][1] + gammaUyy * GammaLUL[2][2][2] + gammaUyz * GammaLUL[2][2][3],
gammaUxz * GammaLUL[2][2][1] + gammaUyz * GammaLUL[2][2][2] + gammaUzz * GammaLUL[2][2][3],
gammaUxz * GammaLUL[2][3][1] + gammaUyz * GammaLUL[2][3][2] + gammaUzz * GammaLUL[2][3][3],
},{gammaUxx * GammaLUL[3][1][1] + gammaUxy * GammaLUL[3][1][2] + gammaUxz * GammaLUL[3][1][3],
gammaUxy * GammaLUL[3][1][1] + gammaUyy * GammaLUL[3][1][2] + gammaUyz * GammaLUL[3][1][3],
gammaUxz * GammaLUL[3][1][1] + gammaUyz * GammaLUL[3][1][2] + gammaUzz * GammaLUL[3][1][3],
gammaUxy * GammaLUL[3][2][1] + gammaUyy * GammaLUL[3][2][2] + gammaUyz * GammaLUL[3][2][3],
gammaUxz * GammaLUL[3][2][1] + gammaUyz * GammaLUL[3][2][2] + gammaUzz * GammaLUL[3][2][3],
gammaUxz * GammaLUL[3][3][1] + gammaUyz * GammaLUL[3][3][2] + gammaUzz * GammaLUL[3][3][3],
},}
local Gamma11SymLL = {
GammaLSymLL[1][1] * GammaLSymUU[1][1] + GammaLSymLL[1][2] * GammaLSymUU[1][2] + GammaLSymLL[1][3] * GammaLSymUU[1][3] + GammaLSymLL[1][2] * GammaLSymUU[1][2] + GammaLSymLL[1][4] * GammaLSymUU[1][4] + GammaLSymLL[1][5] * GammaLSymUU[1][5] + GammaLSymLL[1][3] * GammaLSymUU[1][3] + GammaLSymLL[1][5] * GammaLSymUU[1][5] + GammaLSymLL[1][6] * GammaLSymUU[1][6],
GammaLSymLL[1][1] * GammaLSymUU[2][1] + GammaLSymLL[1][2] * GammaLSymUU[2][2] + GammaLSymLL[1][3] * GammaLSymUU[2][3] + GammaLSymLL[1][2] * GammaLSymUU[2][2] + GammaLSymLL[1][4] * GammaLSymUU[2][4] + GammaLSymLL[1][5] * GammaLSymUU[2][5] + GammaLSymLL[1][3] * GammaLSymUU[2][3] + GammaLSymLL[1][5] * GammaLSymUU[2][5] + GammaLSymLL[1][6] * GammaLSymUU[2][6],
GammaLSymLL[1][1] * GammaLSymUU[3][1] + GammaLSymLL[1][2] * GammaLSymUU[3][2] + GammaLSymLL[1][3] * GammaLSymUU[3][3] + GammaLSymLL[1][2] * GammaLSymUU[3][2] + GammaLSymLL[1][4] * GammaLSymUU[3][4] + GammaLSymLL[1][5] * GammaLSymUU[3][5] + GammaLSymLL[1][3] * GammaLSymUU[3][3] + GammaLSymLL[1][5] * GammaLSymUU[3][5] + GammaLSymLL[1][6] * GammaLSymUU[3][6],
GammaLSymLL[2][1] * GammaLSymUU[2][1] + GammaLSymLL[2][2] * GammaLSymUU[2][2] + GammaLSymLL[2][3] * GammaLSymUU[2][3] + GammaLSymLL[2][2] * GammaLSymUU[2][2] + GammaLSymLL[2][4] * GammaLSymUU[2][4] + GammaLSymLL[2][5] * GammaLSymUU[2][5] + GammaLSymLL[2][3] * GammaLSymUU[2][3] + GammaLSymLL[2][5] * GammaLSymUU[2][5] + GammaLSymLL[2][6] * GammaLSymUU[2][6],
GammaLSymLL[2][1] * GammaLSymUU[3][1] + GammaLSymLL[2][2] * GammaLSymUU[3][2] + GammaLSymLL[2][3] * GammaLSymUU[3][3] + GammaLSymLL[2][2] * GammaLSymUU[3][2] + GammaLSymLL[2][4] * GammaLSymUU[3][4] + GammaLSymLL[2][5] * GammaLSymUU[3][5] + GammaLSymLL[2][3] * GammaLSymUU[3][3] + GammaLSymLL[2][5] * GammaLSymUU[3][5] + GammaLSymLL[2][6] * GammaLSymUU[3][6],
GammaLSymLL[3][1] * GammaLSymUU[3][1] + GammaLSymLL[3][2] * GammaLSymUU[3][2] + GammaLSymLL[3][3] * GammaLSymUU[3][3] + GammaLSymLL[3][2] * GammaLSymUU[3][2] + GammaLSymLL[3][4] * GammaLSymUU[3][4] + GammaLSymLL[3][5] * GammaLSymUU[3][5] + GammaLSymLL[3][3] * GammaLSymUU[3][3] + GammaLSymLL[3][5] * GammaLSymUU[3][5] + GammaLSymLL[3][6] * GammaLSymUU[3][6],
}
local ADL = {
A_x - 2 * D3L[1],
A_y - 2 * D3L[2],
A_z - 2 * D3L[3],
}
local ADU = {
gammaUxx * ADL[1] + gammaUxy * ADL[2] + gammaUxz * ADL[3],
gammaUxy * ADL[1] + gammaUyy * ADL[2] + gammaUyz * ADL[3],
gammaUxz * ADL[1] + gammaUyz * ADL[2] + gammaUzz * ADL[3],
}
local ADDSymLL = {
ADU[1] * (2 * D_xxx) + ADU[2] * (2 * D_xxy) + ADU[3] * (2 * D_xxz),
ADU[1] * (D_xxy + D_yxx) + ADU[2] * (D_xyy + D_yxy) + ADU[3] * (D_xyz + D_yxz),
ADU[1] * (D_xxz + D_zxx) + ADU[2] * (D_xyz + D_zxy) + ADU[3] * (D_xzz + D_zxz),
ADU[1] * (2 * D_yxy) + ADU[2] * (2 * D_yyy) + ADU[3] * (2 * D_yyz),
ADU[1] * (D_yxz + D_zxy) + ADU[2] * (D_yyz + D_zyy) + ADU[3] * (D_yzz + D_zyz),
ADU[1] * (2 * D_zxz) + ADU[2] * (2 * D_zyz) + ADU[3] * (2 * D_zzz),
}
local R4SymLL = {
0,
0,
0,
0,
0,
0,
}
local SL = {
-R4SymLL[1] + trK * K_xx - 2 * KSqSymLL[1] + 4 * D12SymLL[1] + Gamma31SymLL[1] - Gamma11SymLL[1] + ADDSymLL[1] + (A_x * ((2 * V_x) - D1L[1])),
-R4SymLL[2] + trK * K_xy - 2 * KSqSymLL[2] + 4 * D12SymLL[2] + Gamma31SymLL[2] - Gamma11SymLL[2] + ADDSymLL[2] + ((((2 * A_y * V_x) - (A_y * D1L[1])) + ((2 * A_x * V_y) - (A_x * D1L[2]))) / 2),
-R4SymLL[3] + trK * K_xz - 2 * KSqSymLL[3] + 4 * D12SymLL[3] + Gamma31SymLL[3] - Gamma11SymLL[3] + ADDSymLL[3] + ((((2 * A_z * V_x) - (A_z * D1L[1])) + ((2 * A_x * V_z) - (A_x * D1L[3]))) / 2),
-R4SymLL[4] + trK * K_yy - 2 * KSqSymLL[4] + 4 * D12SymLL[4] + Gamma31SymLL[4] - Gamma11SymLL[4] + ADDSymLL[4] + (A_y * ((2 * V_y) - D1L[2])),
-R4SymLL[5] + trK * K_yz - 2 * KSqSymLL[5] + 4 * D12SymLL[5] + Gamma31SymLL[5] - Gamma11SymLL[5] + ADDSymLL[5] + ((((2 * A_z * V_y) - (A_z * D1L[2])) + ((2 * A_y * V_z) - (A_y * D1L[3]))) / 2),
-R4SymLL[6] + trK * K_zz - 2 * KSqSymLL[6] + 4 * D12SymLL[6] + Gamma31SymLL[6] - Gamma11SymLL[6] + ADDSymLL[6] + (A_z * ((2 * V_z) - D1L[3])),
}
local GU0L = {
0,
0,
0,
}
local AKL = {
A_x * KUL[1][1] + A_y * KUL[2][1] + A_z * KUL[3][1],
A_x * KUL[1][2] + A_y * KUL[2][2] + A_z * KUL[3][2],
A_x * KUL[1][3] + A_y * KUL[2][3] + A_z * KUL[3][3],
}
local K12D23L = {
KUL[1][1] * DLUL[1][1][1] +KUL[1][2] * DLUL[1][2][1] +KUL[1][3] * DLUL[1][3][1] + KUL[2][1] * DLUL[1][1][2] +KUL[2][2] * DLUL[1][2][2] +KUL[2][3] * DLUL[1][3][2] + KUL[3][1] * DLUL[1][1][3] +KUL[3][2] * DLUL[1][2][3] +KUL[3][3] * DLUL[1][3][3],
KUL[1][1] * DLUL[2][1][1] +KUL[1][2] * DLUL[2][2][1] +KUL[1][3] * DLUL[2][3][1] + KUL[2][1] * DLUL[2][1][2] +KUL[2][2] * DLUL[2][2][2] +KUL[2][3] * DLUL[2][3][2] + KUL[3][1] * DLUL[2][1][3] +KUL[3][2] * DLUL[2][2][3] +KUL[3][3] * DLUL[2][3][3],
KUL[1][1] * DLUL[3][1][1] +KUL[1][2] * DLUL[3][2][1] +KUL[1][3] * DLUL[3][3][1] + KUL[2][1] * DLUL[3][1][2] +KUL[2][2] * DLUL[3][2][2] +KUL[2][3] * DLUL[3][3][2] + KUL[3][1] * DLUL[3][1][3] +KUL[3][2] * DLUL[3][2][3] +KUL[3][3] * DLUL[3][3][3],
}
local KD23L = {
KUL[1][1] * D1L[1] + KUL[2][1] * D1L[2] + KUL[3][1] * D1L[3],
KUL[1][2] * D1L[1] + KUL[2][2] * D1L[2] + KUL[3][2] * D1L[3],
KUL[1][3] * D1L[1] + KUL[2][3] * D1L[2] + KUL[3][3] * D1L[3],
}
local K12D12L = {
KUL[1][1] * DLUL[1][1][1] + KUL[1][2] * DLUL[1][2][1] + KUL[1][3] * DLUL[1][3][1] + KUL[2][1] * DLUL[2][1][1] + KUL[2][2] * DLUL[2][2][1] + KUL[2][3] * DLUL[2][3][1] + KUL[3][1] * DLUL[3][1][1] + KUL[3][2] * DLUL[3][2][1] + KUL[3][3] * DLUL[3][3][1],
KUL[1][1] * DLUL[1][1][2] + KUL[1][2] * DLUL[1][2][2] + KUL[1][3] * DLUL[1][3][2] + KUL[2][1] * DLUL[2][1][2] + KUL[2][2] * DLUL[2][2][2] + KUL[2][3] * DLUL[2][3][2] + KUL[3][1] * DLUL[3][1][2] + KUL[3][2] * DLUL[3][2][2] + KUL[3][3] * DLUL[3][3][2],
KUL[1][1] * DLUL[1][1][3] + KUL[1][2] * DLUL[1][2][3] + KUL[1][3] * DLUL[1][3][3] + KUL[2][1] * DLUL[2][1][3] + KUL[2][2] * DLUL[2][2][3] + KUL[2][3] * DLUL[2][3][3] + KUL[3][1] * DLUL[3][1][3] + KUL[3][2] * DLUL[3][2][3] + KUL[3][3] * DLUL[3][3][3],
}
local KD12L = {
KUL[1][1] * D3L[1] + KUL[2][1] * D3L[2] + KUL[3][1] * D3L[3],
KUL[1][2] * D3L[1] + KUL[2][2] * D3L[2] + KUL[3][2] * D3L[3],
KUL[1][3] * D3L[1] + KUL[2][3] * D3L[2] + KUL[3][3] * D3L[3],
}
local PL = {
GU0L[1] + AKL[1] - A_x * trK + K12D23L[1] + KD23L[1] - 2 * K12D12L[1] + 2 * KD12L[1],
GU0L[2] + AKL[2] - A_y * trK + K12D23L[2] + KD23L[2] - 2 * K12D12L[2] + 2 * KD12L[2],
GU0L[3] + AKL[3] - A_z * trK + K12D23L[3] + KD23L[3] - 2 * K12D12L[3] + 2 * KD12L[3],
}

		source[i][1] = -alpha * alpha * f * trK
		source[i][2] = -2 * alpha * K_xx
		source[i][3] = -2 * alpha * K_xy
		source[i][4] = -2 * alpha * K_xz
		source[i][5] = -2 * alpha * K_yy
		source[i][6] = -2 * alpha * K_yz
		source[i][7] = -2 * alpha * K_zz
		source[i][29] = alpha * SL[1]
		source[i][30] = alpha * SL[2]
		source[i][31] = alpha * SL[3]
		source[i][32] = alpha * SL[4]
		source[i][33] = alpha * SL[5]
		source[i][34] = alpha * SL[6]
		source[i][35] = alpha * PL[1]
		source[i][36] = alpha * PL[2]
		source[i][37] = alpha * PL[3]
	end
	return source
end

-- enforce constraint V_k = (D_kmn - D_mnk) gamma^mn
function ADM3D:postIterate(sim, qs)
	for i=1,sim.gridsize do
		-- [[ direct assign (seems like this would be constantly overwriting the V_k source term contribution
		local gamma_xx, gamma_xy, gamma_xz, gamma_yy, gamma_yz, gamma_zz = unpack(qs[i], 2, 7)
		local D_xxx, D_xxy, D_xxz, D_xyy, D_xyz, D_xzz = unpack(qs[i], 11, 16)
		local D_yxx, D_yxy, D_yxz, D_yyy, D_yyz, D_yzz = unpack(qs[i], 17, 22)
		local D_zxx, D_zxy, D_zxz, D_zyy, D_zyz, D_zzz = unpack(qs[i], 23, 28)
		local gammaUxx, gammaUxy, gammaUxz, gammaUyy, gammaUyz, gammaUzz = mat33.inv(gamma_xx, gamma_xy, gamma_xz, gamma_yy, gamma_yz, gamma_zz)
		qs[i][35] = 
			(D_xxy - D_yxx) * gammaUxy
			+ (D_xxz - D_zxx) * gammaUxz
			+ (D_xyy - D_yxy) * gammaUyy
			+ (D_xyz - D_yxz) * gammaUyz
			+ (D_xyz - D_zxy) * gammaUyz
			+ (D_xzz - D_zxz) * gammaUzz
		qs[i][36] = 
			(D_yxx - D_xxy) * gammaUxx
			+ (D_yxy - D_xyy) * gammaUxy
			+ (D_yxz - D_xyz) * gammaUxz
			+ (D_yxz - D_zxy) * gammaUxz
			+ (D_yyz - D_zyy) * gammaUyz
			+ (D_yzz - D_zyz) * gammaUzz
		qs[i][37] = 
			(D_zxx - D_xxz) * gammaUxx
			+ (D_zxy - D_xyz) * gammaUxy
			+ (D_zxy - D_yxz) * gammaUxy
			+ (D_zxz - D_xzz) * gammaUxz
			+ (D_zyy - D_yyz) * gammaUyy
			+ (D_zyz - D_yzz) * gammaUyz
		--]]
		--[[ TODO projection 
		--]]
	end
end

return ADM3D
