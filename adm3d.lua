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
	
local function inv3x3sym(xx, xy, xz, yy, yz, zz, det)
	if not det then det = det3x3sym(xx, xy, xz, yy, yz, zz) end
	return {
		(yy * zz - yz^2) / det,	-- xx
		(xz * yz - xy * zz)  / det,	-- xy
		(xy * yz - xz * yy) / det,	-- xz
		(xx * zz - xz^2) / det,	-- yy
		(xz * xy - xx * yz) / det,	-- yz
		(xx * yy - xy^2) / det,	-- zz
	}
end

function ADM3DSimulation:init(args, ...)
	ADM3DSimulation.super.init(self, args, ...)

	local symmath = require 'symmath'

	local x = assert(args.x)
	local y = assert(args.y)
	local z = assert(args.z)
	local vars = table{x,y,z}

	local alpha = symmath.clone(assert(args.alpha)):simplify()
	self.calc_alpha = alpha:compile(vars)

	local diff_alpha = vars:map(function(var) return alpha:diff(var):simplify() end)
	self.calc_diff_alpha = diff_alpha:map(function(eqn) return (eqn:compile(vars)) end)

	local f_param = assert(args.f_param)

	local f = symmath.clone(assert(args.f)):simplify()
	self.calc_f = f:compile{f_param}

	local dalpha_f = f:diff(f_param):simplify()
	self.calc_dalpha_f = dalpha_f:compile{f_param}

	-- pseudo-Cartesian Schwarzschild coordinates from "Catalogue of Spacetimes"
	-- I need to avoid coordinate singularities or something
	local r = (x^2 + y^2 + z^2)^.5

	local rs = assert(tonumber(args.rs))	-- schwarzschild radius

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
	--]]	
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
	local g = table{g_xx, g_xy, g_xz, g_yy, g_yz, g_zz}
	self.calc_g = g:map(function(g_ij) return (g_ij:compile(vars)) end)
	self.calc_diff_g = vars:map(function(x_i, i)
		return g:map(function(g_ij,j)
			return (g_ij:diff(x_i):simplify():compile(vars))
		end)
	end)
	
	local gU = inv3x3sym(g_xx, g_xy, g_xz, g_yy, g_yz, g_zz)
	self.calc_gU = table.map(gU, function(gUij) return (gUij:compile(vars)) end)

	-- and for the graphs ...

	local i=0
	local j=0
	local function col() i=i+1 j=0 end
	self.graphInfos = table()

	local get_state = index:bind(self.qs)
	local function add(args)
		local index = args.index
		for _,suffix in ipairs(args.suffix or {''}) do
			self.graphInfos:insert{
				viewport = {i,j},
				getter = get_state:index(index),
				name = args.name..suffix,
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

	self:buildFields{
		-- not going to worry about flux transform for now, it's only used for error calculations
		fluxTransform = function(i, ...)
			local avgQ = {}
			for j=1,self.numStates do 
				avgQ[j] = (self.qs[i-1][j] + self.qs[i][j]) / 2
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
				0,0,0,			-- V_k
			}
		end,
		eigenfields = function(i, ...)
	
			-- interface eigenfield varialbes
			local avgQ = {}
			for j=1,self.numStates do 
				avgQ[j] = (self.qs[i-1][j] + self.qs[i][j]) / 2
			end
			
			local q_alpha = avgQ[1]
			local q_g_xx, q_g_xy, q_g_xz, q_g_yy, q_g_yz, q_g_zz = unpack(avgQ, 5, 10)
			print('q_g_ij',q_g_xx, q_g_xy, q_g_xz, q_g_yy, q_g_yz, q_g_zz)
			local q_g = det3x3sym(q_g_xx, q_g_xy, q_g_xz, q_g_yy, q_g_yz, q_g_zz)
			print('q_g',q_g)
			local q_A_x, q_A_y, q_A_z = unpack(avgQ, 2, 4)
			local q_D_xxx, q_D_xxy, q_D_xxz, q_D_xyy, q_D_xyz, q_D_xzz = unpack(avgQ, 11, 16)
			local q_D_yxx, q_D_yxy, q_D_yxz, q_D_yyy, q_D_yyz, q_D_yzz = unpack(avgQ, 17, 22)
			local q_D_zxx, q_D_zxy, q_D_zxz, q_D_zyy, q_D_zyz, q_D_zzz = unpack(avgQ, 23, 28)
			local q_K_xx, q_K_xy, q_K_xz, q_K_yy, q_K_yz, q_K_zz = unpack(avgQ, 29, 34)
			local q_V_x, q_V_y, q_V_z = unpack(avgQ, 35, 37)
			local q_gUxx, q_gUxy, q_gUxz, q_gUyy, q_gUyz, q_gUzz = unpack(inv3x3sym(q_g_xx, q_g_xy, q_g_xz, q_g_yy, q_g_yz, q_g_zz))
			local q_f = self.calc_f(q_alpha)
			assert(q_gUxx == q_gUxx)

			-- cell variables
			-- what if, for the ADM equations, there is no distinction?
			-- they're used for Roe's scheme for computing deltas in eigenbasis coordinates by which to scale coordinates coinciding with the lambdas ...
			-- what about creating them solely from 'v' rather than using the average whatsoever?
			-- this would mean ensuring the inputs to the eigenfields() functions were always the state variables themselves (not differences or averages)
			-- 	and deferring differences or averages til after eigenfields() is called (assuming it is a linear function)
			-- this also has an issue with eigenfieldsInverse(), which is called on a flux vector, i.e. at cell interface, which would probably need the average of cells for that input
			local v = {...}
			local v_alpha = avgQ[1]
			local v_g_xx, v_g_xy, v_g_xz, v_g_yy, v_g_yz, v_g_zz = unpack(avgQ, 5, 10)
			local v_A_x, v_A_y, v_A_z = unpack(avgQ, 2, 4)
			local v_D_xxx, v_D_xxy, v_D_xxz, v_D_xyy, v_D_xyz, v_D_xzz = unpack(avgQ, 11, 16)
			local v_D_yxx, v_D_yxy, v_D_yxz, v_D_yyy, v_D_yyz, v_D_yzz = unpack(avgQ, 17, 22)
			local v_D_zxx, v_D_zxy, v_D_zxz, v_D_zyy, v_D_zyz, v_D_zzz = unpack(avgQ, 23, 28)
			local v_K_xx, v_K_xy, v_K_xz, v_K_yy, v_K_yz, v_K_zz = unpack(avgQ, 29, 34)
			local v_V_x, v_V_y, v_V_z = unpack(avgQ, 35, 37)
			local v_gUxx, v_gUxy, v_gUxz, v_gUyy, v_gUyz, v_gUzz = unpack(inv3x3sym(v_g_xx, v_g_xy, v_g_xz, v_g_yy, v_g_yz, v_g_zz))

			return {
				-- negative gauge
				sqrt(q_f) * (q_gUxx * v_K_xx + q_gUxy * v_K_xy + q_gUxz * v_K_xz + q_gUyy * v_K_yy + q_gUyz * v_K_yz + q_gUzz * v_K_zz) - sqrt(q_gUxx) * (v_K_Ax + 2 * (v[35] + q_gUxy/q_gUxx*v[36] + q_gUxz/q_gUxx*v[37])),
				-- negative light cone
				-- zero
				-- positive light cone
				-- positive gauge
			}
		end,
	}
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

	local gU = self.calc_gU:map(function(f) return f(x,y,z) end)

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
		K[1], K[2], K[3], K[4], K[5], K[6],
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
	local gUxx, gUxy, gUxz, gUyy, gUyz, gUzz = unpack(inv3x3sym(g_xx, g_xy, g_xz, g_yy, g_yz, g_zz))
	local f = self.calc_f(alpha)

	local lambdaLight = alpha * sqrt(gUxx)
	local lambdaGauge = lambdaLight * sqrt(f)

	self.eigenvalues[i] = {
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
		lambdaGauge,
	}
end

return ADM3DSimulation
