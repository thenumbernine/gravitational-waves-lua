--[[
this is taking the eigenvectors and doing an inverse,
using those in a Roe solver
just like in my ADM3D

variables are 
alpha		1	1
gamma_ij	6	7
a_i			3	10
d_ijk		18	28
k_ij		6	34
theta		1	35
z_i			3	38
--]]

local class = require 'ext.class'
local Equation = require 'equation'
local mat33sym = require 'mat33sym'
local symmath = require 'symmath'

local Z43D = class(Equation)
Z43D.name = 'Z4-3D'
Z43D.numStates = 38
Z43D.numWaves = 31

function Z43D:init(args, ...)
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
	
	local gammaUxx, gammaUxy, gammaUxz, gammaUyy, gammaUyz, gammaUzz = mat33sym.inv(unpack(exprs.gamma))
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


	-- graphs

	local getters = table()
	local q = function(self,i) return self.qs[i] end
	
	local function add(args)
		local index = args.index
		for _,suffix in ipairs(args.suffix or {''}) do
			getters:insert{[args.name..suffix] = q:_(index)}
			index=index+1
		end
	end

	add{name='alpha', index=1}
	do
		local alpha = q:_(1)
		local gamma_xx = q:_(2)
		local gamma_xy = q:_(3)
		local gamma_xz = q:_(4)
		local gamma_yy = q:_(5)
		local gamma_yz = q:_(6)
		local gamma_zz = q:_(7)
		getters:insert{volume = alpha * math.sqrt:o(mat33sym.det(gamma_xx, gamma_xy, gamma_xz, gamma_yy, gamma_yz, gamma_zz))}
	end
	local suffix3 = {'x', 'y', 'z'}
	local suffix3x3sym = {'xx', 'xy', 'xz', 'yy', 'yz', 'zz'}
	add{name='a_', index=8, suffix=suffix3}
	
	if self.useMomentumConstraints then	
		add{name='V_', index=35, suffix=suffix3}
	else
		if self.useContractedGammaLower then 
			add{name='Gamma_', index=35, suffix=suffix3}
		else
			add{name='Gamma^', index=35, suffix=suffix3}
		end
	end

	add{name='gamma_', index=2, suffix=suffix3x3sym}
	add{name='d_x', index=11, suffix=suffix3x3sym}
	add{name='d_y', index=17, suffix=suffix3x3sym}
	add{name='d_z', index=23, suffix=suffix3x3sym}
	add{name='K_', index=29, suffix=suffix3x3sym}


	self:buildGraphInfos(getters)
end

local function sym3x3(m,i,j)
	local m_xx, m_xy, m_xz, m_yy, m_yz, m_zz = unpack(m)
	return ({
		{m_xx, m_xy, m_xz},
		{m_xy, m_yy, m_yz},
		{m_xz, m_yz, m_zz},
	})[i][j]
end

function Z43D:initCell(solver,i)
	local x = solver.xs[i]
	local y = 0
	local z = 0
	local xs = table{x,y,z}
	local alpha = self.calc.alpha(x,y,z)
	local A = self.calc.A:map(function(a_i) return a_i(x,y,z) end)
	local gamma = self.calc.gamma:map(function(gamma_ij) return gamma_ij(x,y,z) end)
	local D = self.calc.D:map(function(d_i) return d_i:map(function(d_ijk) return d_ijk(x,y,z) end) end)
	local gammaU = self.calc.gammaU:map(function(gammaUij) return gammaUij(x,y,z) end)

	local V = self.useMomentumConstraints and range(3):map(function(i)
		local s = 0
		for j=1,3 do
			for k=1,3 do
				local d_ijk = sym3x3(D[i],j,k)
				local d_kji = sym3x3(D[k],j,i)
				local gammaUjk = sym3x3(gammaU,j,k)
				local dg = (d_ijk - d_kji) * gammaUjk
				s = s + dg
			end
		end
		return s
	end) or nil
	-- Gamma_i = gamma^jk (2 d_jki - d_ijk)
	local GammaL = not self.useMomentumConstraints and range(3):map(function(i)
		local s = 0
		for j=1,3 do
			for k=1,3 do
				local d_jki = sym3x3(D[j], k, i)
				local d_ijk = sym3x3(D[i], j, k)
				local gammaUjk = sym3x3(gammaU, j, k)
				s = s + gammaUjk * (2 * d_jki - d_ijk)
			end
		end
		return s
	end) or nil
	local GammaU = not self.useMomentumConstraints and not self.useContractedGammaLower and range(3):map(function(i)
		local s = 0
		for j=1,3 do
			local gammaUij = sym3x3(gammaU, i, j)
			s = s + gammaUij * GammaL[j]
		end
		return s
	end) or nil

	local K = {}
	for i=1,6 do
		K[i] = self.calc.K[i](x,y,z)
	end

	local Theta = 0
	local Z = {0,0,0}

	-- what do we initialize Theta or Z to?

	return {
		alpha,
		gamma[1], gamma[2], gamma[3], gamma[4], gamma[5], gamma[6],
		A[1], A[2], A[3],
		D[1][1], D[1][2], D[1][3], D[1][4], D[1][5], D[1][6],
		D[2][1], D[2][2], D[2][3], D[2][4], D[2][5], D[2][6],
		D[3][1], D[3][2], D[3][3], D[3][4], D[3][5], D[3][6],
		K[1], K[2], K[3], K[4], K[5], K[6],
		Theta,
		Z[1], Z[2], Z[3],
	}
end

function Z43D:fluxMatrixTransform(solver, avgQ, v)
	-- avgQ is the state used to make the flux 
	-- v is the incoming vector
	return {
		0,	--alpha
		0,0,0,0,0,0,	--gamma_ij
		0,0,0,		-- a_k
		0,0,0,0,0,0,	-- d_xij
		0,0,0,0,0,0,	-- d_yij
		0,0,0,0,0,0,	-- d_zij
		0,0,0,0,0,0,	-- K_ij
		0,				-- Theta
		0,0,0			-- Z_k
	}
end

local lambda = 2

function Z43D:eigenLeftTransform(solver, avgQ, v)
	-- interface eigenvector variables
	local alpha = avgQ[1]
	local gamma_xx, gamma_xy, gamma_xz, gamma_yy, gamma_yz, gamma_zz = unpack(avgQ, 2, 7)
	local gamma = mat33sym.det(gamma_xx, gamma_xy, gamma_xz, gamma_yy, gamma_yz, gamma_zz)
	local a_x, a_y, a_z = unpack(avgQ, 8, 10)
	local d_xxx, d_xxy, d_xxz, d_xyy, d_xyz, d_xzz = unpack(avgQ, 11, 16)
	local d_yxx, d_yxy, d_yxz, d_yyy, d_yyz, d_yzz = unpack(avgQ, 17, 22)
	local d_zxx, d_zxy, d_zxz, d_zyy, d_zyz, d_zzz = unpack(avgQ, 23, 28)
	local K_xx, K_xy, K_xz, K_yy, K_yz, K_zz = unpack(avgQ, 29, 34)
	local Theta = avgQ[35]
	local Z_x, Z_y, Z_z = unpack(avgQ, 36, 38)
	local gammaUxx, gammaUxy, gammaUxz, gammaUyy, gammaUyz, gammaUzz = mat33sym.inv(gamma_xx, gamma_xy, gamma_xz, gamma_yy, gamma_yz, gamma_zz)
	local f = self.calc.f(alpha)

	local aUx = a_x * gammaUxx + a_y * gammaUxy + a_z * gammaUxz

	local K = K_xx * gammaUxx + K_yy * gammaUyy + K_zz * gammaUzz + 2 * (K_xy * gammaUxy + K_xz * gammaUxz + K_yz * gammaUyz)

	local sqrt_f = math.sqrt(f)
	local param1 = (2 - lambda) / (f - 1)
	local param2 = (2 * f - lambda) / (f - 1)

	local d_xmUm = d_xxx * gammaUxx + d_xyy * gammaUyy + d_xzz * gammaUzz + 2 * (d_xxy * gammaUxy + d_xxz * gammaUxz + d_xyz * gammaUyz)
	local d_ymUm = d_yxx * gammaUxx + d_yyy * gammaUyy + d_yzz * gammaUzz + 2 * (d_yxy * gammaUxy + d_yxz * gammaUxz + d_yyz * gammaUyz)
	local d_zmUm = d_zxx * gammaUxx + d_zyy * gammaUyy + d_zzz * gammaUzz + 2 * (d_zxy * gammaUxy + d_zxz * gammaUxz + d_zyz * gammaUyz)

	local dUm_mx = d_xxx * gammaUxx + d_yxy * gammaUyy + d_zxz * gammaUzz + 2 * (d_xxy * gammaUxy + d_xxz * gammaUxz + d_yxz * gammaUyz)
	local dUm_my = d_xxy * gammaUxx + d_yyy * gammaUyy + d_zyz * gammaUzz + 2 * (d_xyy * gammaUxy + d_xyz * gammaUxz + d_yyz * gammaUyz)
	local dUm_mz = d_xxz * gammaUxx + d_yyz * gammaUyy + d_zzz * gammaUzz + 2 * (d_xyz * gammaUxy + d_xzz * gammaUxz + d_yzz * gammaUyz)

	local V_x = d_xmUm - dUm_mx - Z_x
	local V_y = d_ymUm - dUm_my - Z_y
	local V_z = d_zmUm - dUm_mz - Z_z

	local VUx = gammaUxx * V_x + gammaUxy * V_y + gammaUxz * V_z
	local VUy = gammaUxy * V_x + gammaUyy * V_y + gammaUyz * V_z
	local VUz = gammaUxz * V_x + gammaUyz * V_y + gammaUzz * V_z

	--[[
	Lambda^x_ij = d^x_ij + 2 (d_(ij)^x - d^m_m(i delta_j)^x) + delta^x_(i (a_j) - d_j)m^m + V_j))
	Lambda^x_ij = d^x_ij + d_ij^x + d_ji^x - d^m_mi delta_j^x - d^m_mj delta_i^x 
			+ 1/2 ( delta^x_i (a_j - d_jm^m + V_j) + delta^x_j (a_i - d_im^m + V_i) )
	--]]
	local dUx_xy = gammaUxx * d_xxy + gammaUxy * d_yxy + gammaUxz * d_zxy
	local dUx_xz = gammaUxx * d_xxz + gammaUxy * d_yxz + gammaUxz * d_zxz
	local dUx_yy = gammaUxx * d_xyy + gammaUxy * d_yyy + gammaUxz * d_zyy
	local dUx_yz = gammaUxx * d_xyz + gammaUxy * d_yyz + gammaUxz * d_zyz
	local dUx_zz = gammaUxx * d_xzz + gammaUxy * d_yzz + gammaUxz * d_zzz

	local d_xyUx = d_xxy * gammaUxx + d_xyy * gammaUxy + d_xyz * gammaUxz
	local d_xzUx = d_xxz * gammaUxx + d_xyz * gammaUxy + d_xzz * gammaUxz
	local d_yxUx = d_yxx * gammaUxx + d_yxy * gammaUxy + d_yxz * gammaUxz
	local d_yyUx = d_yxy * gammaUxx + d_yyy * gammaUxy + d_yyz * gammaUxz
	local d_yzUx = d_yxz * gammaUxx + d_yyz * gammaUxy + d_yzz * gammaUxz
	local d_zxUx = d_zxx * gammaUxx + d_zxy * gammaUxy + d_zxz * gammaUxz
	local d_zyUx = d_zxy * gammaUxx + d_zyy * gammaUxy + d_zyz * gammaUxz
	local d_zzUx = d_zxz * gammaUxx + d_zyz * gammaUxy + d_zzz * gammaUxz

	local LambdaUx_xy = dUx_xy + d_xyUx + d_yxUx - dUm_my + .5 * (a_y - d_ymUm + V_y)
	local LambdaUx_xz = dUx_xz + d_xzUx + d_zxUx - dUm_mz + .5 * (a_z - d_zmUm + V_z)
	local LambdaUx_yy = dUx_yy + 2 * d_yyUx
	local LambdaUx_yz = dUx_yz + d_yzUx + d_zyUx
	local LambdaUx_zz = dUx_zz + 2 * d_zzUx

	return {
		-- gauge -
		sqrt_f * (K + param1 * Theta) - (aUx + param2 * VUx),
		-- energy -
		Theta - VUx,
		-- light - x5
		K_xy - LambdaUx_xy,
		K_xz - LambdaUx_xz,
		K_yy - LambdaUx_yy,
		K_yz - LambdaUx_yz,
		K_zz - LambdaUx_zz,
		-- time x17
		a_y,
		a_z,
		d_yxx,
		d_yxy,
		d_yxz,
		d_yyy,
		d_yyz,
		d_yzz,
		d_zxx,
		d_zxy,
		d_zxz,
		d_zyy,
		d_zyz,
		d_zzz,
		a_x - f * d_xmUm + lambda * V_x,
		a_y - f * d_ymUm + lambda * V_y,
		a_z - f * d_zmUm + lambda * V_z,
		-- light + x5
		K_xy + LambdaUx_xy,
		K_xz + LambdaUx_xz,
		K_yy + LambdaUx_yy,
		K_yz + LambdaUx_yz,
		K_zz + LambdaUx_zz,
		-- energy +
		Theta + VUx,
		-- gauge +
		sqrt_f * (K + param1 * Theta) + (aUx + param2 * VUx),
	}
end

function Z43D:eigenRightTransform(solver, avgQ, v)
	-- interface eigenvector varialbes
	local alpha = avgQ[1]
	local gamma_xx, gamma_xy, gamma_xz, gamma_yy, gamma_yz, gamma_zz = unpack(avgQ, 2, 7)
	local gamma = mat33sym.det(gamma_xx, gamma_xy, gamma_xz, gamma_yy, gamma_yz, gamma_zz)
	local a_x, a_y, a_z = unpack(avgQ, 8, 10)
	local d_xxx, d_xxy, d_xxz, d_xyy, d_xyz, d_xzz = unpack(avgQ, 11, 16)
	local d_yxx, d_yxy, d_yxz, d_yyy, d_yyz, d_yzz = unpack(avgQ, 17, 22)
	local d_zxx, d_zxy, d_zxz, d_zyy, d_zyz, d_zzz = unpack(avgQ, 23, 28)
	local K_xx, K_xy, K_xz, K_yy, K_yz, K_zz = unpack(avgQ, 29, 34)
	local Theta = avgQ[35]
	local Z_x, Z_y, Z_z = unpack(avgQ, 36, 38)
	local gammaUxx, gammaUxy, gammaUxz, gammaUyy, gammaUyz, gammaUzz = mat33sym.inv(gamma_xx, gamma_xy, gamma_xz, gamma_yy, gamma_yz, gamma_zz)
	local f = self.calc.f(alpha)
	
	local sqrt_f = math.sqrt(f)
	local param1 = (2 - lambda) / (f - 1)
	local param2 = (2 * f - lambda) / (f - 1)
	
	return {
		0,
		0,0,0,0,0,0,
		(-.5 * v[1] 
			+ .5 * param2 * v[2]
			- gammaUxy * v[8]
			- gammaUxz * v[9]
			- .5 * param2 * v[30]
			+ .5 * v[31]) / gammaUxx,
		v[8],
		v[9],
		(-.5 / f * v[1]
			+ .5 * (1 + param2) / f * v[2]
			+ gammaUxy * v[3]
			+ gammaUxz * v[4]
			+ .5 * gammaUyy * v[5]
			+ gammaUyz * v[6]
			+ .5 * gammaUzz * v[7]
			+ 3 * gammaUxy * v[8]
			+ 3 * gammaUxz * v[9]
			- 2 * gammaUxx * gammaUxy * (1 + f) * v[10]
			- 2 * gammaUxy * gammaUxy * (1 + 2 * f) * v[11]
			- 2 * gammaUxy * gammaUxz * (1 + 2 * f) * v[12]
			- gammaUxy * gammaUyy * (1 + 2 * f) * v[13]
			- 2 * gammaUxy * gammaUyz * (1 + 2 * f) * v[14]
			- gammaUxy * gammaUzz * (1 + 2 * f) * v[15]
			- 2 * gammaUxx * gammaUxz * (1 + f) * v[16]
			- 2 * gammaUxy * gammaUxz * (1 + 2 * f) * v[17]
			- 2 * gammaUxz * gammaUxz * (1 + 2 * f) * v[18]
			- gammaUyy * gammaUxz * (1 + 2 * f) * v[19]
			- 2 * gammaUxz * gammaUyz * (1 + 2 * f) * v[20]
			- gammaUxz * gammaUzz * (1 + 2 * f) * v[21]
			- gammaUxx / f * v[22]
			- gammaUxy * (1 / f + 2) * v[23]
			- gammaUxz * (1 / f + 2) * v[24]
			- gammaUxy * v[25]
			- gammaUxz * v[26]
			- .5 * gammaUyy * v[27]
			- gammaUyz * v[28]
			- .5 * gammaUzz * v[29]
			- .5 * (1 + param2) / f * v[30]
			+ .5 / f * v[31]) / (gammaUxx * gammaUxx),
		(-.5 * v[3]
			- 1.5 * v[8]
			+ .5 * (1 + 2 * f) * gammaUxx * v[10]
			+ 2 * gammaUxy * f * v[11]
			+ gammaUxz * (1 + 2 * f) * v[12]
			+ .5 * gammaUyy * (1 + 2 * f) * v[13]
			+ gammaUyz * (1 + 2 * f) * v[14]
			+ .5 * gammaUzz * (1 + 2 * f) * v[15]
			- gammaUxz * v[17]
			+ v[23]
			+ .5 * v[25]) / gammaUxx,
		(-.5 * v[4]
			- 1.5 * v[9]
			- gammaUxy * v[12]
			+ .5 * gammaUxx * (1 + 2 * f) * v[16]
			+ gammaUxy * (1 + 2 * f) * v[17]
			+ 2 * gammaUxz * f * v[18]
			+ .5 * gammaUyy * (1 + 2 * f) * v[19]
			+ gammaUyz * (1 + 2 * f) * v[20]
			+ .5 * gammaUzz * (1 + 2 * f) * v[21]
			+ v[24]
			+ .5 * v[26]) / gammaUxx,
		(-.5 * v[5]
			- gammaUxy * v[13]
			- gammaUxz * v[19]
			+ .5 * v[27]) / gammaUxx,
		(-.5 * v[6]
			- gammaUxy * v[14]
			- gammaUxz * v[10]
			+ .5 * v[28]) / gammaUxx,
		(-.5 * v[7]
			- gammaUxy * v[15]
			- gammaUxz * v[21]
			+ .5 * v[29]) / gammaUxx,
		v[10],
		v[11],
		v[12],
		v[13],
		v[14],
		v[15],
		v[16],
		v[17],
		v[18],
		v[19],
		v[20],
		v[21],
		(.5 / sqrt_f * v[1]
			- .5 * param1 * v[2]
			- gammaUxy * v[3]
			- gammaUxz * v[4]
			- .5 * gammaUyy * v[5]
			- gammaUyz * v[6]
			- gammaUzz * v[7]
			- gammaUxy * v[25]
			- gammaUxz * v[26]
			- .5 * gammaUyy * v[27]
			- gammaUyz * v[28]
			- .5 * gammaUzz * v[29]
			- .5 * param1 * v[30]
			- .5 / sqrt_f * v[31]) / gammaUxx,
		(v[3] + v[25]) * .5,
		(v[4] + v[26]) * .5,
		(v[5] + v[27]) * .5,
		(v[6] + v[28]) * .5,
		(v[7] + v[29]) * .5,
		(v[2] + v[30]) * .5,
		(.5 * v[2]
			- .5 * gammaUxy * v[3]
			- .5 * gammaUxz * v[4]
			- .5 * gammaUyy * v[5]
			- gammaUyz * v[6]
			- .5 * gammaUzz * v[7]
			- .5 * gammaUxy * v[8]
			- .5 * gammaUxz * v[9]
			- .5 * gammaUxx * gammaUxy * v[10]
			- gammaUxx * gammaUyy * v[11]
			- gammaUxx * gammaUyz * v[12]
			- .5 * gammaUxy * gammaUyy * v[13]
			- gammaUxy * gammaUyz * v[14]
			- .5 * gammaUxy * gammaUzz * v[15]
			- .5 * gammaUxx * gammaUxz * v[16]
			- gammaUxx * gammaUyz * v[17]
			- gammaUxx * gammaUzz * v[18]
			- .5 * gammaUxz * gammaUyy * v[19]
			- gammaUxz * gammaUyz * v[20]
			- .5 * gammaUxz * gammaUzz * v[21]
			+ .5 * gammaUxy * v[25]
			+ .5 * gammaUxz * v[26]
			+ .5 * gammaUyy * v[27]
			+ gammaUyz * v[28]
			+ .5 * gammaUzz * v[29]
			- .5 * v[30]) / gammaUxx,
		(.5 * gammaUxx * v[3]
			+ .5 * gammaUxy * v[5]
			+ .5 * gammaUxz * v[6]
			+ .5 * gammaUxx * v[8]
			+ .5 * gammaUxx * gammaUxx * v[10]
			+ gammaUxx * gammaUxy * v[11]
			+ gammaUxx * gammaUxz * v[12]
			- .5 * (gammaUxx * gammaUyy - 2 * gammaUxy * gammaUxy) * v[13]
			+ gammaUxy * gammaUxz * v[14]
			+ .5 * gammaUxy * gammaUzz * v[15]
			- (gammaUxx * gammaUyz - gammaUxz * gammaUxy) * v[19]
			+ (gammaUxz * gammaUxz * gammaUzz * gammaUxx) * v[20]
			- .5 * gammaUxx * v[25]
			- .5 * gammaUxy * v[27]
			- .5 * gammaUxz * v[28]) / gammaUxx,
		(.5 * gammaUxx * v[4]
			+ .5 * gammaUxy * v[6]
			+ .5 * gammaUxz * v[7]
			+ .5 * gammaUxx * v[9]
			- (gammaUxx * gammaUyy - gammaUxy * gammaUxy) * v[14]
			- (gammaUxx * gammaUyz - gammaUxz * gammaUxy) * v[15]
			+ .5 * gammaUxx * gammaUxx * v[16]
			+ gammaUxx * gammaUxy * v[17]
			+ gammaUxx * gammaUxz * v[18]
			+ .5 * gammaUxx * gammaUyy * v[19]
			+ gammaUxy * gammaUxz * v[20]
			- .5 * (gammaUxx * gammaUzz - 2 * gammaUxz * gammaUxz) * v[21]
			- .5 * gammaUxx * v[26]
			- .5 * gammaUxy * v[28]
			- .5 * gammaUxz * v[29]) / gammaUxx,
	}
end

function Z43D:calcEigenBasis(eigenvalues, rightEigenvectors, leftEigenvectors, fluxMatrix, ...)
	fill(eigenvalues, self:calcEigenvaluesFromCons(...))
	fill(leftEigenvectors, ...)
	fill(rightEigenvectors, ...)
	if fluxMatrix then fill(fluxMatrix, ...) end
end

function Z43D:calcEigenvaluesFromCons(
		alpha,
		gamma_xx, gamma_xy, gamma_xz, gamma_yy, gamma_yz, gamma_zz,
		...)
	local gammaUxx, gammaUxy, gammaUxz, gammaUyy, gammaUyz, gammaUzz = mat33sym.inv(gamma_xx, gamma_xy, gamma_xz, gamma_yy, gamma_yz, gamma_zz)
	local f = self.calc.f(alpha)
	local lambdaLight = alpha * math.sqrt(gammaUxx)
	local lambdaGauge = lambdaLight * math.sqrt(f)

	return
		-- gauge field
		-lambdaGauge,
		-- half of 12 along light cones ...
		-lambdaLight,
		-lambdaLight,
		-lambdaLight,
		-lambdaLight,
		-lambdaLight,
		-lambdaLight,
		-- 17 zeroes ...
		0,0,0,0,0,
		0,0,0,0,0,
		0,0,0,0,0,
		0,0,
		-- half of 12 along light cones ...
		lambdaLight,
		lambdaLight,
		lambdaLight,
		lambdaLight,
		lambdaLight,
		lambdaLight,
		-- gauge field
		lambdaGauge
end

function Z43D:sourceTerm(solver, qs)
	local source = solver:newState()
	for i=1,solver.gridsize do
		local alpha = qs[i][1]
		local gamma_xx, gamma_xy, gamma_xz, gamma_yy, gamma_yz, gamma_zz = unpack(qs[i], 2, 7)
		local a_x, a_y, a_z = unpack(qs[i], 8, 10)
		local d_xxx, d_xxy, d_xxz, d_xyy, d_xyz, d_xzz = unpack(qs[i], 11, 16)
		local d_yxx, d_yxy, d_yxz, d_yyy, d_yyz, d_yzz = unpack(qs[i], 17, 22)
		local d_zxx, d_zxy, d_zxz, d_zyy, d_zyz, d_zzz = unpack(qs[i], 23, 28)
		local K_xx, K_xy, K_xz, K_yy, K_yz, K_zz = unpack(qs[i], 29, 34)
		local Theta = qs[i][35]
		local Z_x, Z_y, Z_z = unpack(qs[i], 36, 38)
		local gammaUxx, gammaUxy, gammaUxz, gammaUyy, gammaUyz, gammaUzz = mat33sym.inv(gamma_xx, gamma_xy, gamma_xz, gamma_yy, gamma_yz, gamma_zz)
		local f = self.calc.f(alpha)
		local K = K_xx * gammaUxx + K_yy * gammaUyy + K_zz * gammaUzz + 2 * (K_xy * gammaUxy + K_xz * gammaUxz + K_yz * gammaUyz)
		
		source[i][1] = -alpha * alpha * f * K
		source[i][2] = -2 * alpha * K_xx
		source[i][3] = -2 * alpha * K_xy
		source[i][4] = -2 * alpha * K_xz
		source[i][5] = -2 * alpha * K_yy
		source[i][6] = -2 * alpha * K_yz
		source[i][7] = -2 * alpha * K_zz
	
		-- how about the other source terms?
	end
	return source
end

return Z43D 
