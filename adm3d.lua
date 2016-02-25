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
partial_t g_ij = -2 alpha K_ij

	flux terms:

partial_t A_k + alpha f g^ij partial_k K_ij = -alpha tr K (f + alpha f') A_k + 2 alpha f K^ij D_kij
partial_t D_kij + alpha partial_k K_ij = -alpha A_k K_ij
partial_t K_ij + alpha (g^km partial_k D_mij + 1/2 (delta^k_i partial_k A_j + delta^k_j partial_k A_i) + delta^k_i partial_k V_j + delta^k_j partial_k V_i) = alpha S_ij - alpha lambda^k_ij A_k + 2 alpha D_mij D_k^km
partial_t V_k = alpha P_k

	eigenfields:

lambda = 0:
	alpha
	g_ij
	V_i
	A_x'
	D_x'ij
	A_x - f D_x^m_m

lambda = +- alpha sqrt(g^xx):
	K_ix' +- sqrt(g_xx) (D_xix' + delta^x_i V_x' / g^xx) for x' != x

lambda = +- alpha sqrt(f g^xx):
	sqrt(f) tr(K) +- sqrt(g_xx) (A_x + 2 V^x / g^xx)

--]]
local class = require 'ext.class'
local Equation = require 'equation'
local mat33 = require 'mat33'
local symmath = require 'symmath'

local ADM3D = class(Equation)
ADM3D.name = 'ADM3D'

--[[
alpha,
g_xx, g_xy, g_xz, g_yy, g_yz, g_zz,
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
	exprs.g = table{
		assert(makesym'g_xx'),	--xx
		makesym'g_xy' or symmath.Constant(0),	--xy
		makesym'g_xz' or symmath.Constant(0),	--xz
		makesym'g_yy' or symmath.Constant(1),	--yy
		makesym'g_yz' or symmath.Constant(0),	--yz
		makesym'g_zz' or symmath.Constant(1),	--zz
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
	
	local gUxx, gUxy, gUxz, gUyy, gUyz, gUzz = mat33.inv(unpack(exprs.g))
	exprs.gU = table{gUxx, gUxy, gUxz, gUyy, gUyz, gUzz}
	exprs.D = vars:map(function(x_k)
		return exprs.g:map(function(g_ij)
			return (g_ij:diff(x_k)/2):simplify()
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
	local get_g_xx = get_state(2)
	local get_A_x = get_state(8)
	local get_D_xxx = get_state(11)
	local get_K_xx = get_state(29)
	self.graphInfos = {
		{viewport={0/3, 0/3, 1/3, 1/3}, getter=get_alpha, name='alpha', color={1,0,1}},
		{viewport={0/3, 1/3, 1/3, 1/3}, getter=get_A_x, name='A_x', color={0,1,0}},
		{viewport={1/3, 0/3, 1/3, 1/3}, getter=get_g_xx, name='g_xx', color={.5,.5,1}},
		{viewport={1/3, 1/3, 1/3, 1/3}, getter=get_D_xxx, name='D_xxx', color={1,1,0}},
		{viewport={2/3, 0/3, 1/3, 1/3}, getter=get_K_xx, name='K_xx', color={0,1,1}},
		{viewport={2/3, 1/3, 1/3, 1/3}, getter=get_alpha * sqrt:compose(get_g_xx), name='volume', color={0,1,1}},
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
		local g_xx = state:index(2)
		local g_xy = state:index(3)
		local g_xz = state:index(4)
		local g_yy = state:index(5)
		local g_yz = state:index(6)
		local g_zz = state:index(7)
		self.graphInfos:insert{
			viewport = {i,j},
			getter = alpha * mat33.det(g_xx, g_xy, g_xz, g_yy, g_yz, g_zz),
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
	add{name='g_', index=2, suffix=suffix3x3sym}
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
		{viewport={1/3, 0/3, 1/3, 1/3}, getter=function(self,i) return self.qs[i][2] end, name='g_xx', color={.5,.5,1}},
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
		0,0,0,0,0,0,	--g_ij
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
	local g_xx, g_xy, g_xz, g_yy, g_yz, g_zz = unpack(avgQ, 2, 7)
	local g = mat33.det(g_xx, g_xy, g_xz, g_yy, g_yz, g_zz)
	local A_x, A_y, A_z = unpack(avgQ, 8, 10)
	local D_xxx, D_xxy, D_xxz, D_xyy, D_xyz, D_xzz = unpack(avgQ, 11, 16)
	local D_yxx, D_yxy, D_yxz, D_yyy, D_yyz, D_yzz = unpack(avgQ, 17, 22)
	local D_zxx, D_zxy, D_zxz, D_zyy, D_zyz, D_zzz = unpack(avgQ, 23, 28)
	local K_xx, K_xy, K_xz, K_yy, K_yz, K_zz = unpack(avgQ, 29, 34)
	local V_x, V_y, V_z = unpack(avgQ, 35, 37)
	local gUxx, gUxy, gUxz, gUyy, gUyz, gUzz = mat33.inv(g_xx, g_xy, g_xz, g_yy, g_yz, g_zz)
	local f = self.calc.f(alpha)

	-- cell variables
	-- what if, for the ADM equations, there is no distinction?
	-- they're used for Roe's scheme for computing deltas in eigenbasis coordinates by which to scale coordinates coinciding with the lambdas ...
	-- what about creating them solely from 'v' rather than using the average whatsoever?
	-- this would mean ensuring the inputs to the eigenfields() functions were always the state variables themselves (not differences or averages)
	-- 	and deferring differences or averages til after eigenfields() is called (assuming it is a linear function)
	-- this also has an issue with eigenfieldsInverse(), which is called on a flux vector, i.e. at cell interface, which would probably need the average of cells for that input

	return {
		((((-(2 * gUxz * v[37])) - (gUxx * v[8])) + (math.sqrt(f) * (gUxx ^ (3 / 2)) * v[29]) + (2 * math.sqrt(f) * gUxy * v[30] * math.sqrt(gUxx)) + (2 * math.sqrt(f) * gUxz * v[31] * math.sqrt(gUxx)) + (math.sqrt(f) * gUyy * v[32] * math.sqrt(gUxx)) + (2 * math.sqrt(f) * gUyz * v[33] * math.sqrt(gUxx)) + (((math.sqrt(f) * gUzz * v[34] * math.sqrt(gUxx)) - (2 * gUxx * v[35])) - (2 * gUxy * v[36]))) / math.sqrt(gUxx)),
		(((-(gUxx * v[12])) + ((v[30] * math.sqrt(gUxx)) - v[36])) / math.sqrt(gUxx)),
		(((-(gUxx * v[13])) + ((v[31] * math.sqrt(gUxx)) - v[37])) / math.sqrt(gUxx)),
		((-(math.sqrt(gUxx) * v[14])) + v[32]),
		((-(math.sqrt(gUxx) * v[15])) + v[33]),
		((-(math.sqrt(gUxx) * v[16])) + v[34]),
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
		((((((v[8] - (f * gUxx * v[11])) - (2 * f * gUxy * v[12])) - (2 * f * gUxz * v[13])) - (f * gUyy * v[14])) - (2 * f * gUyz * v[15])) - (f * gUzz * v[16])),
		(((gUxx * v[12]) + (v[30] * math.sqrt(gUxx)) + v[36]) / math.sqrt(gUxx)),
		(((gUxx * v[13]) + (v[31] * math.sqrt(gUxx)) + v[37]) / math.sqrt(gUxx)),
		((math.sqrt(gUxx) * v[14]) + v[32]),
		((math.sqrt(gUxx) * v[15]) + v[33]),
		((math.sqrt(gUxx) * v[16]) + v[34]),
		(((2 * gUxz * v[37]) + (gUxx * v[8]) + (math.sqrt(f) * (gUxx ^ (3 / 2)) * v[29]) + (2 * math.sqrt(f) * gUxy * v[30] * math.sqrt(gUxx)) + (2 * math.sqrt(f) * gUxz * v[31] * math.sqrt(gUxx)) + (math.sqrt(f) * gUyy * v[32] * math.sqrt(gUxx)) + (2 * math.sqrt(f) * gUyz * v[33] * math.sqrt(gUxx)) + (math.sqrt(f) * gUzz * v[34] * math.sqrt(gUxx)) + (2 * gUxx * v[35]) + (2 * gUxy * v[36])) / math.sqrt(gUxx))
	}

end

function ADM3D:eigenfieldsInverse(sim, i, v)
	
	-- interface eigenfield varialbes
	local avgQ = {}
	for j=1,sim.numStates do 
		avgQ[j] = (sim.qs[i-1][j] + sim.qs[i][j]) / 2
	end
	local alpha = avgQ[1]
	local g_xx, g_xy, g_xz, g_yy, g_yz, g_zz = unpack(avgQ, 2, 7)
	local g = mat33.det(g_xx, g_xy, g_xz, g_yy, g_yz, g_zz)
	local A_x, A_y, A_z = unpack(avgQ, 8, 10)
	local D_xxx, D_xxy, D_xxz, D_xyy, D_xyz, D_xzz = unpack(avgQ, 11, 16)
	local D_yxx, D_yxy, D_yxz, D_yyy, D_yyz, D_yzz = unpack(avgQ, 17, 22)
	local D_zxx, D_zxy, D_zxz, D_zyy, D_zyz, D_zzz = unpack(avgQ, 23, 28)
	local K_xx, K_xy, K_xz, K_yy, K_yz, K_zz = unpack(avgQ, 29, 34)
	local V_x, V_y, V_z = unpack(avgQ, 35, 37)
	local gUxx, gUxy, gUxz, gUyy, gUyz, gUzz = mat33.inv(g_xx, g_xy, g_xz, g_yy, g_yz, g_zz)
	local f = self.calc.f(alpha)

	return {
		v[7],
		v[8],
		v[9],
		v[10],
		v[11],
		v[12],
		v[13],
		(((-v[37]) + (4 * gUxz * v[30] * (1 / math.sqrt(gUxx))) + (4 * gUxy * v[29] * (1 / math.sqrt(gUxx))) + (4 * v[28] * math.sqrt(gUxx)) + v[1]) / (-(2 * math.sqrt(gUxx)))),
		v[14],
		v[15],
		(((-v[37]) + (gUzz * v[36] * f) + (2 * gUyz * v[35] * f) + (gUyy * v[34] * f) + (2 * gUxz * v[33] * f) + (2 * gUxy * v[32] * f) + (2 * v[31] * math.sqrt(gUxx)) + ((4 * gUxz * v[30] * (1 / math.sqrt(gUxx))) - (4 * gUxz * f * v[30] * (1 / math.sqrt(gUxx)))) + ((4 * gUxy * v[29] * (1 / math.sqrt(gUxx))) - (4 * gUxy * f * v[29] * (1 / math.sqrt(gUxx)))) + ((((((4 * v[28] * math.sqrt(gUxx)) - (gUzz * v[6] * f)) - (2 * gUyz * v[5] * f)) - (gUyy * v[4] * f)) - (2 * gUxz * v[3] * f)) - (2 * gUxy * v[2] * f)) + v[1]) / (-(2 * (gUxx ^ (3 / 2)) * f))),
		(((-v[32]) + (2 * v[29] * (1 / math.sqrt(gUxx))) + v[2]) / (-(2 * math.sqrt(gUxx)))),
		(((-v[33]) + (2 * v[30] * (1 / math.sqrt(gUxx))) + v[3]) / (-(2 * math.sqrt(gUxx)))),
		(((-v[34]) + v[4]) / (-(2 * math.sqrt(gUxx)))),
		(((-v[35]) + v[5]) / (-(2 * math.sqrt(gUxx)))),
		(((-v[36]) + v[6]) / (-(2 * math.sqrt(gUxx)))),
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
		((((((((((((v[37] - (gUzz * v[36] * math.sqrt(f))) - (2 * gUyz * v[35] * math.sqrt(f))) - (gUyy * v[34] * math.sqrt(f))) - (2 * gUxz * v[33] * math.sqrt(f))) - (2 * gUxy * v[32] * math.sqrt(f))) - (gUzz * v[6] * math.sqrt(f))) - (2 * gUyz * v[5] * math.sqrt(f))) - (gUyy * v[4] * math.sqrt(f))) - (2 * gUxz * v[3] * math.sqrt(f))) - (2 * gUxy * v[2] * math.sqrt(f))) + v[1]) / (2 * math.sqrt(f) * gUxx)),
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
	local g = self.calc.g:map(function(g_ij) return g_ij(x,y,z) end)
	local D = self.calc.D:map(function(D_i) return D_i:map(function(D_ijk) return D_ijk(x,y,z) end) end)
	local gU = self.calc.gU:map(function(gUij) return gUij(x,y,z) end)

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
				local gUjk = sym3x3(gU,j,k)
				local dg = (D_ijk - D_kji) * gUjk
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
		g[1], g[2], g[3], g[4], g[5], g[6],
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
	local g_xx, g_xy, g_xz, g_yy, g_yz, g_zz = unpack(avgQ, 2, 7)
	local A_x, A_y, A_z = unpack(avgQ, 8, 10)
	local D_xxx, D_xxy, D_xxz, D_xyy, D_xyz, D_xzz = unpack(avgQ, 11, 16)
	local D_yxx, D_yxy, D_yxz, D_yyy, D_yyz, D_yzz = unpack(avgQ, 17, 22)
	local D_zxx, D_zxy, D_zxz, D_zyy, D_zyz, D_zzz = unpack(avgQ, 23, 28)
	local K_xx, K_xy, K_xz, K_yy, K_yz, K_zz = unpack(avgQ, 29, 34)
	local V_x, V_y, V_z = unpack(avgQ, 35, 37)
	local gUxx, gUxy, gUxz, gUyy, gUyz, gUzz = mat33.inv(g_xx, g_xy, g_xz, g_yy, g_yz, g_zz)
	local f = self.calc.f(alpha)

	local lambdaLight = alpha * sqrt(gUxx)
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
		local g_xx, g_xy, g_xz, g_yy, g_yz, g_zz = unpack(qs[i], 2, 7)
		local A_x, A_y, A_z = unpack(qs[i], 8, 10)
		local D_xxx, D_xxy, D_xxz, D_xyy, D_xyz, D_xzz = unpack(qs[i], 11, 16)
		local D_yxx, D_yxy, D_yxz, D_yyy, D_yyz, D_yzz = unpack(qs[i], 17, 22)
		local D_zxx, D_zxy, D_zxz, D_zyy, D_zyz, D_zzz = unpack(qs[i], 23, 28)
		local K_xx, K_xy, K_xz, K_yy, K_yz, K_zz = unpack(qs[i], 29, 34)
		local V_x, V_y, V_z = unpack(qs[i], 35, 37)
		local gUxx, gUxy, gUxz, gUyy, gUyz, gUzz = mat33.inv(g_xx, g_xy, g_xz, g_yy, g_yz, g_zz)
		local f = self.calc.f(alpha)

		local tr_K = 
			K_xx * gUxx
			+ K_yy * gUyy
			+ K_zz * gUzz
			+ 2 * K_xy * gUxy
			+ 2 * K_yz * gUyz
			+ 2 * K_xz * gUxz

		source[i][1] = -alpha * alpha * f * tr_K
		source[i][2] = -2 * alpha * K_xx
		source[i][3] = -2 * alpha * K_xy
		source[i][4] = -2 * alpha * K_xz
		source[i][5] = -2 * alpha * K_yy
		source[i][6] = -2 * alpha * K_yz
		source[i][7] = -2 * alpha * K_zz
		-- TODO dq_dts[i][29..34] == partial_t K_ij += alpha * S_ij
		-- partial_t V_k = alpha * P_k
	end
	return source
end

-- enforce constraint V_k = (D_kmn - D_mnk) g^mn
function ADM3D:postIterate(sim, qs)
	for i=1,sim.gridsize do
		-- [[ direct assign (seems like this would be constantly overwriting the V_k source term contribution
		local g_xx, g_xy, g_xz, g_yy, g_yz, g_zz = unpack(qs[i], 2, 7)
		local D_xxx, D_xxy, D_xxz, D_xyy, D_xyz, D_xzz = unpack(qs[i], 11, 16)
		local D_yxx, D_yxy, D_yxz, D_yyy, D_yyz, D_yzz = unpack(qs[i], 17, 22)
		local D_zxx, D_zxy, D_zxz, D_zyy, D_zyz, D_zzz = unpack(qs[i], 23, 28)
		local gUxx, gUxy, gUxz, gUyy, gUyz, gUzz = mat33.inv(g_xx, g_xy, g_xz, g_yy, g_yz, g_zz)
		qs[i][35] = 
			(D_xxy - D_yxx) * gUxy
			+ (D_xxz - D_zxx) * gUxz
			+ (D_xyy - D_yxy) * gUyy
			+ (D_xyz - D_yxz) * gUyz
			+ (D_xyz - D_zxy) * gUyz
			+ (D_xzz - D_zxz) * gUzz
		qs[i][36] = 
			(D_yxx - D_xxy) * gUxx
			+ (D_yxy - D_xyy) * gUxy
			+ (D_yxz - D_xyz) * gUxz
			+ (D_yxz - D_zxy) * gUxz
			+ (D_yyz - D_zyy) * gUyz
			+ (D_yzz - D_zyz) * gUzz
		qs[i][37] = 
			(D_zxx - D_xxz) * gUxx
			+ (D_zxy - D_xyz) * gUxy
			+ (D_zxy - D_yxz) * gUxy
			+ (D_zxz - D_xzz) * gUxz
			+ (D_zyy - D_yyz) * gUyy
			+ (D_zyz - D_yzz) * gUyz
		--]]
		--[[ TODO make a giant linear system out of it and project it to the nullspace 
		--]]
	end
end

return ADM3D
