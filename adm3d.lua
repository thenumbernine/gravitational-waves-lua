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

eigenfields:

lambda = 0:
	alpha
	g_ij
	A_x'	<- what is the x <-> x' transformation?
	D_x'ij
	V_i
	A_x - f D_x^m_m

lambda = +- alpha sqrt(g^xx):
	K_ix' +- sqrt(g_xx) (D_xix' + delta^x_i V_x' / g^xx) for x' != x (what does that mean? what is x'?)

lambda = +- alpha sqrt(f g^xx):
	sqrt(f) tr(K) +- sqrt(g_xx) (A_x + 2 V^x / g^xx)

eigenvector decomposition
--]]
require 'ext'
local Simulation = require 'simulation'
local ADM3DSimulation = class(Simulation)

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
ADM3DSimulation.numStates = 37

local function det3x3sym(xx, xy, xz, yy, yz, zz)
	return xx * yy * zz
		+ xy * yz * xz
		+ xz * xy * yz
		- xz * yy * xz
		- yz * yz * xx
		- zz * xy * xy
end

function ADM3DSimulation:init(args, ...)
	ADM3DSimulation.super.init(self, args, ...)

	local x = assert(args.x)
	local y = assert(args.y)
	local z = assert(args.z)
	local vars = table{x,y,z}

	local symmath = require 'symmath'

	local alpha = symmath.clone(assert(args.alpha)):simplify()
	self.calc_alpha = alpha:compile{x,y,z}

	local diff_alpha = vars:map(function(var) return alpha:diff(var):simplify() end)
	self.calc_diff_alpha = diff_alpha:map(function(eqn) return eqn:compile(vars) end)

	-- pseudo-Cartesian Schwarzschild coordinates from "Catalogue of Spacetimes"
	-- I need to avoid coordinate singularities or something
	local r = (x^2 + y^2 + z^2)^.5

	local rs = assert(args.rs)	-- schwarzschild radius

	local g_xx = (x^2 / (1 - rs/r) + y^2 + z^2) / r^2
	local g_yy = (x^2 + y^2 / (1 - rs/r) + z^2) / r^2
	local g_zz = (x^2 + y^2 + z^2 / (1 - rs/r)) / r^2
	local g_xy = x * y * rs / (r^2 * (r - rs))
	local g_xz = x * z * rs / (r^2 * (r - rs))
	local g_yz = y * z * rs / (r^2 * (r - rs))
	local g = table{g_xx, g_xy, g_xz, g_yy, g_yz, g_zz}
	self.calc_g = g:map(function(g_ij) return g_ij:compile(vars) end)
	self.calc_diff_g = vars:map(function(x_i)
		return g:map(function(g_ij)
			return g_ij:diff(x_i):simplify():compile(vars)
		end)
	end)

	-- local gInvMat = symmat.Matrix.inverse( matrix of g ):simplify()
	-- or just do it manually ...
	local det = det3x3sym(g_xx, g_xy, g_xz, g_yy, g_yz, g_zz)
	local gInv = {
		(g_yy * g_zz - g_yz^2) / det,	-- xx
		(g_xz * g_yz - g_xy * g_zz)  / det,	-- xy
		(g_xy * g_yz - g_xz * g_yy) / det,	-- xz
		(g_xx * g_zz - g_xz^2) / det,	-- yy
		(g_xz * g_xy - g_xx * g_yz) / det,	-- yz
		(g_xx * g_yy - g_xy^2) / det,	-- zz
	}
	self.calc_gInv = gInv:map(function(gInv_ij) return gInv_ij:compile(vars) end)

	-- and for the graphs ...

	local i=0
	local j=0
	local function col() i=i+1 j=0 end
	self.graphInfos = table()

	local get_state = index:bind(self.qs)
	local function add(args)
		local index = args.index
		for k=1,(args.count or 1) do
			self.graphInfos:insert{
				viewport = {i,j},
				getter = get_state:index(index),
				name = args.name,
				color = {0,1,1},
			}
			index=index+1
			j=j+1
		end
	end

	add{name='alpha', index=1}
	do
		local state = index:bind(self.qs)
		local alpha = state:index(1)
		local g_xx = state:index(2)
		local g_xy = state:index(3)
		local g_xz = state:index(4)
		local g_yy = state:index(5)
		local g_yz = state:index(6)
		local g_zz = state:index(7)
		self.graphInfos:insert{
			viewport = {i,j},
			getter = alpha * det3x3sym(g_xx, g_xy, g_xz, g_yy, g_yz, g_zz),
			name = 'volume', 
			color = {0,1,1},
		}
	end
	j=j+1
	self.graphInfos:insert{
		viewport = {i,j},
		getter = log:compose(index:bind(self.eigenbasisErrors)),
		name = 'log eigenbasis error', 
		color = {1,0,0}, 
		range = {-30, 30},
	}
	j=j+1
	self.graphInfos:insert{
		viewport = {i,j},
		getter = log:compose(index:bind(self.fluxMatrixErrors)),
		name = 'log reconstruction error',
		color = {1,0,0},
		range = {-30, 30},
	}
	col()
	add{name='A', index=8, count=3}
	add{name='V', index=35, count=3}
	col()
	add{name='g', index=2, count=6}
	col()
	add{name='D_x', index=11, count=6}
	col()
	add{name='D_y', index=17, count=6}
	col()
	add{name='D_z', index=23, count=6}
	col()
	add{name='K', index=29, count=6}
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
end

function ADM3DSimulation:initCell(i)
	local x = self.xs[i]
	local y = 0
	local z = 0
	local xs = table{x,y,z}
	local alpha = self.calc_alpha(x)
	local diff_alpha = self.calc_diff_alpha:map(function(f) return f(x,y,z) end)
	local A = diff_alpha:map(function(diff_alpha_i) return diff_alpha_i / alpha end)
	local g = self.calc_g:map(function(f) return f(x,y,z) end)
	local D = self.calc_diff_g:map(function(fs)
		return fs:map(function(f) return .5 * f(x,y,z) end)
	end)

	local gInv = self.calc_gInv:map(function(f) return f(x,y,z) end)

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
				s = s + (sym3x3(D[i],j,k) - sym3x3(D[k],j,i)) * sym3x3(gInv,j,k)
			end
		end
		return s
	end)

	-- hmm ... 
	local K = {}
	for i=1,6 do
		K[i] = 0
	end

	return {
		alpha,
		g[1], g[2], g[3], g[4], g[5], g[6],
		A[1], A[2], A[3],
		D[1][1], D[1][2], D[1][3], D[1][4], D[1][5], D[1][6],
		D[2][1], D[2][2], D[2][3], D[2][4], D[2][5], D[2][6],
		D[3][1], D[3][2], D[3][3], D[3][4], D[3][5], D[3][6],
		K_xx, K_xy, K_xz, K_yy, K_yz, K_zz,
		V[1], V[2], V[3],
	}
end

function ADM3DSimulation:calcInterfaceEigenBasis(i)
	local avgQ = {}
	for j=1,self.numStates do
		avgQ[j] = (self.qs[i-1][j] + self.qs[i][j]) / 2
	end

	local alpha = avgQ[1]
	local g_xx, g_xy, g_xz, g_yy, g_yz, g_zz = unpack(avgQ, 5, 10)
	local A_x, A_y, A_z = unpack(avgQ, 2, 4)
	local D_xxx, D_xxy, D_xxz, D_xyy, D_xyz, D_xzz = unpack(avgQ, 11, 16)
	local D_yxx, D_yxy, D_yxz, D_yyy, D_yyz, D_yzz = unpack(avgQ, 17, 22)
	local D_zxx, D_zxy, D_zxz, D_zyy, D_zyz, D_zzz = unpack(avgQ, 23, 28)
	local K_xx, K_xy, K_xz, K_yy, K_yz, K_zz = unpack(avgQ, 29, 34)
	local V_x, V_y, V_z = unpack(avgQ, 35, 37)

end

