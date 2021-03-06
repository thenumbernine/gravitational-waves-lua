--[[
see the symmath/tests/numerical_relativity_codegen.lua for more info
--]]
local class = require 'ext.class'
local Equation = require 'equation'
local mat33sym = require 'mat33sym'
local symmath = require 'symmath'

local ADM3D = class(Equation)
ADM3D.name = 'ADM3D'

--[[
alpha,
gamma_xx, gamma_xy, gamma_xz, gamma_yy, gamma_yz, gamma_zz,
a_x, a_y, a_z,
d_xxx, d_xxy, d_xxz, d_xyy, d_xyz, d_xzz,
d_yxx, d_yxy, d_yxz, d_yyy, d_yyz, d_yzz,
d_zxx, d_zxy, d_zxz, d_zyy, d_zyz, d_zzz,
K_xx, K_xy, K_xz, K_yy, K_yz, K_zz,
V_x, V_y, V_z,
--]]
ADM3D.numStates = 37

--[[
true means use V_i (requires constraints as well)
false means use Gamma^i (doesn't require constraints)
	Gamma^i = gamma^ij Gamma^i_jk = gamma^jk gamma^il Gamma_ljk
		= 1/2 gamma^jk gamma^il (gamma_lj,k + gamma_lk,j - gamma_jk,l)
		= 1/2 (gamma^il gamma_lj,k gamma^jk + gamma^il gamma_lk,j gamma^kj - gamma^jk gamma_jk,l gamma^li)
		= -(gamma^ij_,j - 1/2 gamma_jk gamma^jk_,l gamma^li)
		= gamma^jk (2 d_jk^i - d^i_jk)
		= gamma^il gamma^jk (2 d_jkl - d_ljk)

0 = 3,k = delta^i_i,j = (gamma^ij gamma_ij),k = gamma^ij_,k gamma_ij + gamma^ij gamma_ij,k 
so gamma^ij_,k gamma_ij = -gamma^ij gamma_ij,k
--]]

ADM3D.useMomentumConstraints = true	-- advect V_i
--ADM3D.useMomentumConstraints = false	-- advect Gamma^i

-- if useMomentumConstraints == true
ADM3D.momentumConstraintMethod = 'directAssign'	--directly assign V_i = d_im^m + d^m_mi.  technically this is ignoring the V_i's altogether ... but this matches the ADM1D5Var solver
--ADM3D.momentumConstraintMethod = 'linearProject' 	-- linear project V_i and d_ijk

-- if useMomentumConstraints == false
-- true means use Gamma_i, false means use Gamma^i
--ADM3D.useContractedGammaLower = true	-- technically not hyperbolic ... and maybe that's why it's not working?
ADM3D.useContractedGammaLower = false

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
		if symmath.Expression:isa(expr) then
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

--[[ match the 1D 3-var layout:
	getters = table{
		{alpha = function(self,i) return self.qs[i][1] end},
		{a_x = function(self,i) return self.qs[i][8] end},
		{gamma_xx = function(self,i) return self.qs[i][2] end},
		{d_xxx = function(self,i) return self.qs[i][11] end},
		{K_xx = function(self,i) return self.qs[i][29] end},
		{volume = function(self,i) return self.qs[i][1] * math.sqrt(self.qs[i][2]) end},
	}
--]]

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

function ADM3D:initCell(solver,i)
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

	local last = self.useMomentumConstraints and V or (self.useContractedGammaLower and GammaL or GammaU)

	return {
		alpha,
		gamma[1], gamma[2], gamma[3], gamma[4], gamma[5], gamma[6],
		A[1], A[2], A[3],
		D[1][1], D[1][2], D[1][3], D[1][4], D[1][5], D[1][6],
		D[2][1], D[2][2], D[2][3], D[2][4], D[2][5], D[2][6],
		D[3][1], D[3][2], D[3][3], D[3][4], D[3][5], D[3][6],
		K[1], K[2], K[3], K[4], K[5], K[6],
		last[1], last[2], last[3],
	}
end

function ADM3D:fluxMatrixTransform(solver, avgQ, v)
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
		0,0,0			-- V_k
	}
end

function ADM3D:eigenLeftTransform(solver, avgQ, v)
	-- interface eigenvector variables
	local alpha = avgQ[1]
	local gamma_xx, gamma_xy, gamma_xz, gamma_yy, gamma_yz, gamma_zz = unpack(avgQ, 2, 7)
	local gamma = mat33sym.det(gamma_xx, gamma_xy, gamma_xz, gamma_yy, gamma_yz, gamma_zz)
	local a_x, a_y, a_z = unpack(avgQ, 8, 10)
	local d_xxx, d_xxy, d_xxz, d_xyy, d_xyz, d_xzz = unpack(avgQ, 11, 16)
	local d_yxx, d_yxy, d_yxz, d_yyy, d_yyz, d_yzz = unpack(avgQ, 17, 22)
	local d_zxx, d_zxy, d_zxz, d_zyy, d_zyz, d_zzz = unpack(avgQ, 23, 28)
	local K_xx, K_xy, K_xz, K_yy, K_yz, K_zz = unpack(avgQ, 29, 34)
	--local V_x, V_y, V_z = unpack(avgQ, 35, 37)
	local gammaUxx, gammaUxy, gammaUxz, gammaUyy, gammaUyz, gammaUzz = mat33sym.inv(gamma_xx, gamma_xy, gamma_xz, gamma_yy, gamma_yz, gamma_zz)
	local f = self.calc.f(alpha)

	-- cell variables
	-- what if, for the ADM equations, there is no distinction?
	-- they're used for Roe's scheme for computing deltas in eigenbasis coordinates by which to scale coordinates coinciding with the lambdas ...
	-- what about creating them solely from 'v' rather than using the average whatsoever?
	-- this would mean ensuring the inputs to the applyLeftEigenvectors() functions were always the state variables themselves (not differences or averages)
	-- 	and deferring differences or averages til after applyLeftEigenvectors() is called (assuming it is a linear function)
	-- this also has an issue with applyRightEigenvectors(), which is called on a flux vector, i.e. at cell interface, which would probably need the average of cells for that input

	if self.useMomentumConstraints then
		-- left eigenvectors in x:
		return {
			((math.sqrt(f) * (gammaUxx ^ (3 / 2)) * v[29]) + (2 * math.sqrt(f) * math.sqrt(gammaUxx) * gammaUxy * v[30]) + (2 * math.sqrt(f) * math.sqrt(gammaUxx) * gammaUxz * v[31]) + (math.sqrt(f) * math.sqrt(gammaUxx) * gammaUyy * v[32]) + (2 * math.sqrt(f) * math.sqrt(gammaUxx) * gammaUyz * v[33]) + (((((((math.sqrt(f) * math.sqrt(gammaUxx) * gammaUzz * v[34]) - (gammaUxx * v[8])) - (2 * gammaUxx * v[35])) - (gammaUxy * v[9])) - (2 * gammaUxy * v[36])) - (gammaUxz * v[10])) - (2 * gammaUxz * v[37]))),
			((-(v[9] + (2 * v[36]) + ((((2 * gammaUxx * v[12]) - (gammaUxx * v[17])) - (2 * math.sqrt(gammaUxx) * v[30])) - (2 * gammaUxz * v[19])) + ((((2 * gammaUxz * v[24]) - (gammaUyy * v[20])) - (2 * gammaUyz * v[21])) - (gammaUzz * v[22])))) / 2),
			((-(v[10] + (2 * v[37]) + (((2 * gammaUxx * v[13]) - (gammaUxx * v[23])) - (2 * math.sqrt(gammaUxx) * v[31])) + (((((2 * gammaUxy * v[19]) - (2 * gammaUxy * v[24])) - (gammaUyy * v[26])) - (2 * gammaUyz * v[27])) - (gammaUzz * v[28])))) / 2),
			(-(((gammaUxx * v[14]) - (math.sqrt(gammaUxx) * v[32])) + (gammaUxy * v[20]) + (gammaUxz * v[26]))),
			(-(((gammaUxx * v[15]) - (math.sqrt(gammaUxx) * v[33])) + (gammaUxy * v[21]) + (gammaUxz * v[27]))),
			(-(((gammaUxx * v[16]) - (math.sqrt(gammaUxx) * v[34])) + (gammaUxy * v[22]) + (gammaUxz * v[28]))),
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
			((v[9] + (2 * v[36]) + ((2 * gammaUxx * v[12]) - (gammaUxx * v[17])) + ((2 * math.sqrt(gammaUxx) * v[30]) - (2 * gammaUxz * v[19])) + ((((2 * gammaUxz * v[24]) - (gammaUyy * v[20])) - (2 * gammaUyz * v[21])) - (gammaUzz * v[22]))) / 2),
			((v[10] + (2 * v[37]) + ((2 * gammaUxx * v[13]) - (gammaUxx * v[23])) + (2 * math.sqrt(gammaUxx) * v[31]) + (((((2 * gammaUxy * v[19]) - (2 * gammaUxy * v[24])) - (gammaUyy * v[26])) - (2 * gammaUyz * v[27])) - (gammaUzz * v[28]))) / 2),
			((gammaUxx * v[14]) + (math.sqrt(gammaUxx) * v[32]) + (gammaUxy * v[20]) + (gammaUxz * v[26])),
			((gammaUxx * v[15]) + (math.sqrt(gammaUxx) * v[33]) + (gammaUxy * v[21]) + (gammaUxz * v[27])),
			((gammaUxx * v[16]) + (math.sqrt(gammaUxx) * v[34]) + (gammaUxy * v[22]) + (gammaUxz * v[28])),
			((math.sqrt(f) * (gammaUxx ^ (3 / 2)) * v[29]) + (2 * math.sqrt(f) * math.sqrt(gammaUxx) * gammaUxy * v[30]) + (2 * math.sqrt(f) * math.sqrt(gammaUxx) * gammaUxz * v[31]) + (math.sqrt(f) * math.sqrt(gammaUxx) * gammaUyy * v[32]) + (2 * math.sqrt(f) * math.sqrt(gammaUxx) * gammaUyz * v[33]) + (math.sqrt(f) * math.sqrt(gammaUxx) * gammaUzz * v[34]) + (gammaUxx * v[8]) + (2 * gammaUxx * v[35]) + (gammaUxy * v[9]) + (2 * gammaUxy * v[36]) + (gammaUxz * v[10]) + (2 * gammaUxz * v[37]))
		}
	else
		if self.useContractedGammaLower then
			return {
				((math.sqrt(f) * (gammaUxx ^ (3 / 2)) * v[29]) + (2 * math.sqrt(f) * math.sqrt(gammaUxx) * gammaUxy * v[30]) + (2 * math.sqrt(f) * math.sqrt(gammaUxx) * gammaUxz * v[31]) + (math.sqrt(f) * math.sqrt(gammaUxx) * gammaUyy * v[32]) + (2 * math.sqrt(f) * math.sqrt(gammaUxx) * gammaUyz * v[33]) + (((math.sqrt(f) * math.sqrt(gammaUxx) * gammaUzz * v[34]) - (gammaUxx * v[8])) - ((gammaUxx ^ 2) * v[11])) + (((((((gammaUxx * v[35]) - (2 * gammaUxx * gammaUxy * v[12])) - (gammaUxx * gammaUxy * v[17])) - (2 * gammaUxx * gammaUxz * v[13])) - (gammaUxx * gammaUxz * v[23])) - (gammaUxy * v[9])) - (2 * (gammaUxy ^ 2) * v[18])) + ((((gammaUxy * v[36]) - (2 * gammaUxy * gammaUxz * v[19])) - (gammaUxz * v[10])) - (2 * (gammaUxz ^ 2) * v[25])) + (((((((((((gammaUxz * v[37]) - (2 * gammaUxz * gammaUxy * v[24])) - (gammaUyy * gammaUxx * v[14])) - (gammaUyy * gammaUxy * v[20])) - (gammaUyy * gammaUxz * v[26])) - (2 * gammaUyz * gammaUxx * v[15])) - (2 * gammaUyz * gammaUxy * v[21])) - (2 * gammaUyz * gammaUxz * v[27])) - (gammaUzz * gammaUxx * v[16])) - (gammaUzz * gammaUxy * v[22])) - (gammaUzz * gammaUxz * v[28]))),
				((-((v[9] - v[36]) + ((2 * gammaUxx * v[12]) - (2 * math.sqrt(gammaUxx) * v[30])) + (2 * gammaUxy * v[18]) + (2 * gammaUxz * v[24]))) / 2),
				((-((v[10] - v[37]) + ((2 * gammaUxx * v[13]) - (2 * math.sqrt(gammaUxx) * v[31])) + (2 * gammaUxy * v[19]) + (2 * gammaUxz * v[25]))) / 2),
				(-(((gammaUxx * v[14]) - (math.sqrt(gammaUxx) * v[32])) + (gammaUxy * v[20]) + (gammaUxz * v[26]))),
				(-(((gammaUxx * v[15]) - (math.sqrt(gammaUxx) * v[33])) + (gammaUxy * v[21]) + (gammaUxz * v[27]))),
				(-(((gammaUxx * v[16]) - (math.sqrt(gammaUxx) * v[34])) + (gammaUxy * v[22]) + (gammaUxz * v[28]))),
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
				((-((((((v[35] - (gammaUxx * v[11])) - (2 * gammaUxy * v[12])) - (2 * gammaUxz * v[13])) - (gammaUyy * v[14])) - (2 * gammaUyz * v[15])) - (gammaUzz * v[16]))) / 2),
				((-((((((v[36] - (gammaUxx * v[17])) - (2 * gammaUxy * v[18])) - (2 * gammaUxz * v[19])) - (gammaUyy * v[20])) - (2 * gammaUyz * v[21])) - (gammaUzz * v[22]))) / 2),
				((-((((((v[37] - (gammaUxx * v[23])) - (2 * gammaUxy * v[24])) - (2 * gammaUxz * v[25])) - (gammaUyy * v[26])) - (2 * gammaUyz * v[27])) - (gammaUzz * v[28]))) / 2),
				((((((v[8] - (f * gammaUxx * v[11])) - (2 * f * gammaUxy * v[12])) - (2 * f * gammaUxz * v[13])) - (f * gammaUyy * v[14])) - (2 * f * gammaUyz * v[15])) - (f * gammaUzz * v[16])),
				(((v[9] - v[36]) + (2 * gammaUxx * v[12]) + (2 * math.sqrt(gammaUxx) * v[30]) + (2 * gammaUxy * v[18]) + (2 * gammaUxz * v[24])) / 2),
				(((v[10] - v[37]) + (2 * gammaUxx * v[13]) + (2 * math.sqrt(gammaUxx) * v[31]) + (2 * gammaUxy * v[19]) + (2 * gammaUxz * v[25])) / 2),
				((gammaUxx * v[14]) + (math.sqrt(gammaUxx) * v[32]) + (gammaUxy * v[20]) + (gammaUxz * v[26])),
				((gammaUxx * v[15]) + (math.sqrt(gammaUxx) * v[33]) + (gammaUxy * v[21]) + (gammaUxz * v[27])),
				((gammaUxx * v[16]) + (math.sqrt(gammaUxx) * v[34]) + (gammaUxy * v[22]) + (gammaUxz * v[28])),
				((math.sqrt(f) * (gammaUxx ^ (3 / 2)) * v[29]) + (2 * math.sqrt(f) * math.sqrt(gammaUxx) * gammaUxy * v[30]) + (2 * math.sqrt(f) * math.sqrt(gammaUxx) * gammaUxz * v[31]) + (math.sqrt(f) * math.sqrt(gammaUxx) * gammaUyy * v[32]) + (2 * math.sqrt(f) * math.sqrt(gammaUxx) * gammaUyz * v[33]) + (math.sqrt(f) * math.sqrt(gammaUxx) * gammaUzz * v[34]) + (gammaUxx * v[8]) + (((gammaUxx ^ 2) * v[11]) - (gammaUxx * v[35])) + (2 * gammaUxx * gammaUxy * v[12]) + (gammaUxx * gammaUxy * v[17]) + (2 * gammaUxx * gammaUxz * v[13]) + (gammaUxx * gammaUxz * v[23]) + (gammaUxy * v[9]) + ((2 * (gammaUxy ^ 2) * v[18]) - (gammaUxy * v[36])) + (2 * gammaUxy * gammaUxz * v[19]) + (gammaUxz * v[10]) + ((2 * (gammaUxz ^ 2) * v[25]) - (gammaUxz * v[37])) + (2 * gammaUxz * gammaUxy * v[24]) + (gammaUyy * gammaUxx * v[14]) + (gammaUyy * gammaUxy * v[20]) + (gammaUyy * gammaUxz * v[26]) + (2 * gammaUyz * gammaUxx * v[15]) + (2 * gammaUyz * gammaUxy * v[21]) + (2 * gammaUyz * gammaUxz * v[27]) + (gammaUzz * gammaUxx * v[16]) + (gammaUzz * gammaUxy * v[22]) + (gammaUzz * gammaUxz * v[28]))
			}
		else
			return {
				(v[35] + (math.sqrt(f) * (gammaUxx ^ (3 / 2)) * v[29]) + (2 * math.sqrt(f) * math.sqrt(gammaUxx) * gammaUxy * v[30]) + (2 * math.sqrt(f) * math.sqrt(gammaUxx) * gammaUxz * v[31]) + (math.sqrt(f) * math.sqrt(gammaUxx) * gammaUyy * v[32]) + (2 * math.sqrt(f) * math.sqrt(gammaUxx) * gammaUyz * v[33]) + ((((((((((((((((((((((math.sqrt(f) * math.sqrt(gammaUxx) * gammaUzz * v[34]) - (gammaUxx * v[8])) - ((gammaUxx ^ 2) * v[11])) - (2 * gammaUxx * gammaUxy * v[12])) - (gammaUxx * gammaUxy * v[17])) - (2 * gammaUxx * gammaUxz * v[13])) - (gammaUxx * gammaUxz * v[23])) - (gammaUxy * v[9])) - (2 * (gammaUxy ^ 2) * v[18])) - (2 * gammaUxy * gammaUxz * v[19])) - (gammaUxz * v[10])) - (2 * (gammaUxz ^ 2) * v[25])) - (2 * gammaUxz * gammaUxy * v[24])) - (gammaUyy * gammaUxx * v[14])) - (gammaUyy * gammaUxy * v[20])) - (gammaUyy * gammaUxz * v[26])) - (2 * gammaUyz * gammaUxx * v[15])) - (2 * gammaUyz * gammaUxy * v[21])) - (2 * gammaUyz * gammaUxz * v[27])) - (gammaUzz * gammaUxx * v[16])) - (gammaUzz * gammaUxy * v[22])) - (gammaUzz * gammaUxz * v[28]))),
				((-(v[9] + ((2 * gammaUxx * v[12]) - (2 * math.sqrt(gammaUxx) * v[30])) + (2 * gammaUxy * v[18]) + ((((2 * gammaUxz * v[24]) - (v[3] * v[35])) - (v[5] * v[36])) - (v[6] * v[37])))) / 2),
				((-(v[10] + ((2 * gammaUxx * v[13]) - (2 * math.sqrt(gammaUxx) * v[31])) + (2 * gammaUxy * v[19]) + ((((2 * gammaUxz * v[25]) - (v[4] * v[35])) - (v[6] * v[36])) - (v[7] * v[37])))) / 2),
				(-(((gammaUxx * v[14]) - (math.sqrt(gammaUxx) * v[32])) + (gammaUxy * v[20]) + (gammaUxz * v[26]))),
				(-(((gammaUxx * v[15]) - (math.sqrt(gammaUxx) * v[33])) + (gammaUxy * v[21]) + (gammaUxz * v[27]))),
				(-(((gammaUxx * v[16]) - (math.sqrt(gammaUxx) * v[34])) + (gammaUxy * v[22]) + (gammaUxz * v[28]))),
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
				(((gammaUxx * v[11]) + (2 * gammaUxy * v[12]) + (2 * gammaUxz * v[13]) + (gammaUyy * v[14]) + (2 * gammaUyz * v[15]) + ((((gammaUzz * v[16]) - (v[2] * v[35])) - (v[3] * v[36])) - (v[4] * v[37]))) / 2),
				(((gammaUxx * v[17]) + (2 * gammaUxy * v[18]) + (2 * gammaUxz * v[19]) + (gammaUyy * v[20]) + (2 * gammaUyz * v[21]) + ((((gammaUzz * v[22]) - (v[3] * v[35])) - (v[5] * v[36])) - (v[6] * v[37]))) / 2),
				(((gammaUxx * v[23]) + (2 * gammaUxy * v[24]) + (2 * gammaUxz * v[25]) + (gammaUyy * v[26]) + (2 * gammaUyz * v[27]) + ((((gammaUzz * v[28]) - (v[4] * v[35])) - (v[6] * v[36])) - (v[7] * v[37]))) / 2),
				((((((v[8] - (f * gammaUxx * v[11])) - (2 * f * gammaUxy * v[12])) - (2 * f * gammaUxz * v[13])) - (f * gammaUyy * v[14])) - (2 * f * gammaUyz * v[15])) - (f * gammaUzz * v[16])),
				((v[9] + (2 * gammaUxx * v[12]) + (2 * math.sqrt(gammaUxx) * v[30]) + (2 * gammaUxy * v[18]) + ((((2 * gammaUxz * v[24]) - (v[3] * v[35])) - (v[5] * v[36])) - (v[6] * v[37]))) / 2),
				((v[10] + (2 * gammaUxx * v[13]) + (2 * math.sqrt(gammaUxx) * v[31]) + (2 * gammaUxy * v[19]) + ((((2 * gammaUxz * v[25]) - (v[4] * v[35])) - (v[6] * v[36])) - (v[7] * v[37]))) / 2),
				((gammaUxx * v[14]) + (math.sqrt(gammaUxx) * v[32]) + (gammaUxy * v[20]) + (gammaUxz * v[26])),
				((gammaUxx * v[15]) + (math.sqrt(gammaUxx) * v[33]) + (gammaUxy * v[21]) + (gammaUxz * v[27])),
				((gammaUxx * v[16]) + (math.sqrt(gammaUxx) * v[34]) + (gammaUxy * v[22]) + (gammaUxz * v[28])),
				(-(((((((((((((((((((((((((((v[35] - (math.sqrt(f) * (gammaUxx ^ (3 / 2)) * v[29])) - (2 * math.sqrt(f) * math.sqrt(gammaUxx) * gammaUxy * v[30])) - (2 * math.sqrt(f) * math.sqrt(gammaUxx) * gammaUxz * v[31])) - (math.sqrt(f) * math.sqrt(gammaUxx) * gammaUyy * v[32])) - (2 * math.sqrt(f) * math.sqrt(gammaUxx) * gammaUyz * v[33])) - (math.sqrt(f) * math.sqrt(gammaUxx) * gammaUzz * v[34])) - (gammaUxx * v[8])) - ((gammaUxx ^ 2) * v[11])) - (2 * gammaUxx * gammaUxy * v[12])) - (gammaUxx * gammaUxy * v[17])) - (2 * gammaUxx * gammaUxz * v[13])) - (gammaUxx * gammaUxz * v[23])) - (gammaUxy * v[9])) - (2 * (gammaUxy ^ 2) * v[18])) - (2 * gammaUxy * gammaUxz * v[19])) - (gammaUxz * v[10])) - (2 * (gammaUxz ^ 2) * v[25])) - (2 * gammaUxz * gammaUxy * v[24])) - (gammaUyy * gammaUxx * v[14])) - (gammaUyy * gammaUxy * v[20])) - (gammaUyy * gammaUxz * v[26])) - (2 * gammaUyz * gammaUxx * v[15])) - (2 * gammaUyz * gammaUxy * v[21])) - (2 * gammaUyz * gammaUxz * v[27])) - (gammaUzz * gammaUxx * v[16])) - (gammaUzz * gammaUxy * v[22])) - (gammaUzz * gammaUxz * v[28])))
			}
		end
	end
end

function ADM3D:eigenRightTransform(solver, avgQ, v)
	-- interface eigenvector varialbes
	local alpha = avgQ[1]
	local gamma_xx, gamma_xy, gamma_xz, gamma_yy, gamma_yz, gamma_zz = unpack(avgQ, 2, 7)
	local gamma = mat33sym.det(gamma_xx, gamma_xy, gamma_xz, gamma_yy, gamma_yz, gamma_zz)
	local a_x, a_y, a_z = unpack(avgQ, 8, 10)
	local d_xxx, d_xxy, d_xxz, d_xyy, d_xyz, d_xzz = unpack(avgQ, 11, 16)
	local d_yxx, d_yxy, d_yxz, d_yyy, d_yyz, d_yzz = unpack(avgQ, 17, 22)
	local d_zxx, d_zxy, d_zxz, d_zyy, d_zyz, d_zzz = unpack(avgQ, 23, 28)
	local K_xx, K_xy, K_xz, K_yy, K_yz, K_zz = unpack(avgQ, 29, 34)
	--local V_x, V_y, V_z = unpack(avgQ, 35, 37)
	local gammaUxx, gammaUxy, gammaUxz, gammaUyy, gammaUyz, gammaUzz = mat33sym.inv(gamma_xx, gamma_xy, gamma_xz, gamma_yy, gamma_yz, gamma_zz)
	local f = self.calc.f(alpha)

	-- right eigenvectors in x:
	if self.useMomentumConstraints then
		return {
			v[7],
			v[8],
			v[9],
			v[10],
			v[11],
			v[12],
			v[13],
			((((4 * v[28] * gammaUxx) - v[37]) + v[1] + (2 * gammaUxy * v[14]) + (4 * gammaUxy * v[29]) + (2 * gammaUxz * v[15]) + (4 * gammaUxz * v[30])) / (-(2 * gammaUxx))),
			v[14],
			v[15],
			((-((4 * v[28] * gammaUxx) + ((2 * v[31] * gammaUxx) - v[37]) + v[1] + (2 * gammaUxy * v[14]) + (2 * gammaUxy * v[16] * gammaUxx * f) + (4 * gammaUxy * v[29]) + ((((2 * gammaUxy * v[32] * f) - (2 * gammaUxy * f * v[14])) - (4 * gammaUxy * f * v[29])) - (2 * gammaUxy * v[2] * f)) + (2 * gammaUxz * v[15]) + (2 * gammaUxz * v[22] * gammaUxx * f) + (4 * gammaUxz * v[30]) + ((((2 * gammaUxz * v[33] * f) - (2 * gammaUxz * f * v[15])) - (4 * gammaUxz * f * v[30])) - (2 * gammaUxz * v[3] * f)) + ((gammaUyy * v[34] * f) - (gammaUyy * v[4] * f)) + ((2 * gammaUyz * v[35] * f) - (2 * gammaUyz * v[5] * f)) + ((gammaUzz * v[36] * f) - (gammaUzz * v[6] * f)))) / (2 * (gammaUxx ^ 2) * f)),
			(((v[14] - (v[16] * gammaUxx)) + (((2 * v[29]) - v[32]) - (2 * gammaUxz * v[18])) + ((((2 * gammaUxz * v[23]) - (gammaUyy * v[19])) - (2 * gammaUyz * v[20])) - (gammaUzz * v[21])) + v[2]) / (-(2 * gammaUxx))),
			(((v[15] - (v[22] * gammaUxx)) + ((2 * v[30]) - v[33]) + (((((2 * gammaUxy * v[18]) - (2 * gammaUxy * v[23])) - (gammaUyy * v[25])) - (2 * gammaUyz * v[26])) - (gammaUzz * v[27])) + v[3]) / (-(2 * gammaUxx))),
			((((v[34] - (2 * gammaUxy * v[19])) - (2 * gammaUxz * v[25])) - v[4]) / (2 * gammaUxx)),
			((((v[35] - (2 * gammaUxy * v[20])) - (2 * gammaUxz * v[26])) - v[5]) / (2 * gammaUxx)),
			((((v[36] - (2 * gammaUxy * v[21])) - (2 * gammaUxz * v[27])) - v[6]) / (2 * gammaUxx)),
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
			((v[37] + ((((((((((v[1] - (2 * gammaUxy * v[32] * math.sqrt(f))) - (2 * gammaUxy * v[2] * math.sqrt(f))) - (2 * gammaUxz * v[33] * math.sqrt(f))) - (2 * gammaUxz * v[3] * math.sqrt(f))) - (gammaUyy * v[34] * math.sqrt(f))) - (gammaUyy * v[4] * math.sqrt(f))) - (2 * gammaUyz * v[35] * math.sqrt(f))) - (2 * gammaUyz * v[5] * math.sqrt(f))) - (gammaUzz * v[36] * math.sqrt(f))) - (gammaUzz * v[6] * math.sqrt(f)))) / (2 * math.sqrt(f) * (gammaUxx ^ (3 / 2)))),
			((v[32] + v[2]) / (2 * math.sqrt(gammaUxx))),
			((v[33] + v[3]) / (2 * math.sqrt(gammaUxx))),
			((v[34] + v[4]) / (2 * math.sqrt(gammaUxx))),
			((v[35] + v[5]) / (2 * math.sqrt(gammaUxx))),
			((v[36] + v[6]) / (2 * math.sqrt(gammaUxx))),
			v[28],
			v[29],
			v[30]
		}
	else
		if self.useContractedGammaLower then
			return {
				v[7],
				v[8],
				v[9],
				v[10],
				v[11],
				v[12],
				v[13],
				v[14],
				v[15],
				(((((((v[37] - (4 * v[30] * gammaUzz)) - v[1]) - (2 * gammaUxz * v[14])) - (4 * gammaUxz * v[28])) - (2 * gammaUyz * v[15])) - (4 * gammaUyz * v[29])) / (2 * gammaUzz)),
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
				((((v[32] - (2 * gammaUxz * v[16])) - (2 * gammaUyz * v[22])) - v[2]) / (2 * gammaUzz)),
				((((v[33] - (2 * gammaUxz * v[17])) - (2 * gammaUyz * v[23])) - v[3]) / (2 * gammaUzz)),
				(((v[14] - (v[21] * gammaUzz)) + ((((((2 * v[28]) - v[34]) - (gammaUxx * v[16])) - (2 * gammaUxy * v[17])) - (gammaUyy * v[19])) - (2 * gammaUyz * v[20])) + (2 * gammaUyz * v[24]) + v[4]) / (-(2 * gammaUzz))),
				((((v[35] - (2 * gammaUxz * v[19])) - (2 * gammaUyz * v[25])) - v[5]) / (2 * gammaUzz)),
				((((v[15] - (v[27] * gammaUzz)) - v[36]) + (((2 * v[29]) - (gammaUxx * v[22])) - (2 * gammaUxy * v[23])) + (((2 * gammaUxz * v[20]) - (2 * gammaUxz * v[24])) - (gammaUyy * v[25])) + v[6]) / (-(2 * gammaUzz))),
				((((((v[37] - (4 * v[30] * gammaUzz)) - (2 * v[31] * gammaUzz)) - v[1]) - (gammaUxx * v[32] * f)) + ((gammaUxx * v[2] * f) - (2 * gammaUxy * v[33] * f)) + (((((2 * gammaUxy * v[3] * f) - (2 * gammaUxz * v[14])) - (2 * gammaUxz * v[21] * gammaUzz * f)) - (4 * gammaUxz * v[28])) - (2 * gammaUxz * v[34] * f)) + (2 * gammaUxz * f * v[14]) + (4 * gammaUxz * f * v[28]) + ((2 * gammaUxz * v[4] * f) - (gammaUyy * v[35] * f)) + (((((gammaUyy * v[5] * f) - (2 * gammaUyz * v[15])) - (2 * gammaUyz * v[27] * gammaUzz * f)) - (2 * gammaUyz * v[36] * f)) - (4 * gammaUyz * v[29])) + (2 * gammaUyz * f * v[15]) + (4 * gammaUyz * f * v[29]) + (2 * gammaUyz * v[6] * f)) / (2 * (gammaUzz ^ 2) * f)),
				((v[32] + v[2]) / (2 * math.sqrt(gammaUzz))),
				((v[33] + v[3]) / (2 * math.sqrt(gammaUzz))),
				((v[34] + v[4]) / (2 * math.sqrt(gammaUzz))),
				((v[35] + v[5]) / (2 * math.sqrt(gammaUzz))),
				((v[36] + v[6]) / (2 * math.sqrt(gammaUzz))),
				((v[37] + ((((((((((v[1] - (gammaUxx * v[32] * math.sqrt(f))) - (gammaUxx * v[2] * math.sqrt(f))) - (2 * gammaUxy * v[33] * math.sqrt(f))) - (2 * gammaUxy * v[3] * math.sqrt(f))) - (2 * gammaUxz * v[34] * math.sqrt(f))) - (2 * gammaUxz * v[4] * math.sqrt(f))) - (gammaUyy * v[35] * math.sqrt(f))) - (gammaUyy * v[5] * math.sqrt(f))) - (2 * gammaUyz * v[36] * math.sqrt(f))) - (2 * gammaUyz * v[6] * math.sqrt(f)))) / (2 * math.sqrt(f) * (gammaUzz ^ (3 / 2)))),
				(-(((((((2 * v[28]) - (gammaUxx * v[16])) - (2 * gammaUxy * v[17])) - (2 * gammaUxz * v[18])) - (gammaUyy * v[19])) - (2 * gammaUyz * v[20])) - (gammaUzz * v[21]))),
				(-(((((((2 * v[29]) - (gammaUxx * v[22])) - (2 * gammaUxy * v[23])) - (2 * gammaUxz * v[24])) - (gammaUyy * v[25])) - (2 * gammaUyz * v[26])) - (gammaUzz * v[27]))),
				(((((((((v[37] - (4 * v[30] * gammaUzz)) - (2 * v[31] * gammaUzz)) - v[1]) - (4 * f * v[30] * gammaUzz)) - (2 * gammaUxz * v[14])) - (4 * gammaUxz * v[28])) - (2 * gammaUyz * v[15])) - (4 * gammaUyz * v[29])) / (2 * f * gammaUzz))
			}
		else
			return {
				(v[35] + (math.sqrt(f) * (gammaUxx ^ (3 / 2)) * v[29]) + (2 * math.sqrt(f) * math.sqrt(gammaUxx) * gammaUxy * v[30]) + (2 * math.sqrt(f) * math.sqrt(gammaUxx) * gammaUxz * v[31]) + (math.sqrt(f) * math.sqrt(gammaUxx) * gammaUyy * v[32]) + (2 * math.sqrt(f) * math.sqrt(gammaUxx) * gammaUyz * v[33]) + ((((((((((((((((((((((math.sqrt(f) * math.sqrt(gammaUxx) * gammaUzz * v[34]) - (gammaUxx * v[8])) - ((gammaUxx ^ 2) * v[11])) - (2 * gammaUxx * gammaUxy * v[12])) - (gammaUxx * gammaUxy * v[17])) - (2 * gammaUxx * gammaUxz * v[13])) - (gammaUxx * gammaUxz * v[23])) - (gammaUxy * v[9])) - (2 * (gammaUxy ^ 2) * v[18])) - (2 * gammaUxy * gammaUxz * v[19])) - (gammaUxz * v[10])) - (2 * (gammaUxz ^ 2) * v[25])) - (2 * gammaUxz * gammaUxy * v[24])) - (gammaUyy * gammaUxx * v[14])) - (gammaUyy * gammaUxy * v[20])) - (gammaUyy * gammaUxz * v[26])) - (2 * gammaUyz * gammaUxx * v[15])) - (2 * gammaUyz * gammaUxy * v[21])) - (2 * gammaUyz * gammaUxz * v[27])) - (gammaUzz * gammaUxx * v[16])) - (gammaUzz * gammaUxy * v[22])) - (gammaUzz * gammaUxz * v[28]))),
				((-(v[9] + ((2 * gammaUxx * v[12]) - (2 * math.sqrt(gammaUxx) * v[30])) + (2 * gammaUxy * v[18]) + ((((2 * gammaUxz * v[24]) - (v[3] * v[35])) - (v[5] * v[36])) - (v[6] * v[37])))) / 2),
				((-(v[10] + ((2 * gammaUxx * v[13]) - (2 * math.sqrt(gammaUxx) * v[31])) + (2 * gammaUxy * v[19]) + ((((2 * gammaUxz * v[25]) - (v[4] * v[35])) - (v[6] * v[36])) - (v[7] * v[37])))) / 2),
				(-(((gammaUxx * v[14]) - (math.sqrt(gammaUxx) * v[32])) + (gammaUxy * v[20]) + (gammaUxz * v[26]))),
				(-(((gammaUxx * v[15]) - (math.sqrt(gammaUxx) * v[33])) + (gammaUxy * v[21]) + (gammaUxz * v[27]))),
				(-(((gammaUxx * v[16]) - (math.sqrt(gammaUxx) * v[34])) + (gammaUxy * v[22]) + (gammaUxz * v[28]))),
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
				(((gammaUxx * v[11]) + (2 * gammaUxy * v[12]) + (2 * gammaUxz * v[13]) + (gammaUyy * v[14]) + (2 * gammaUyz * v[15]) + ((((gammaUzz * v[16]) - (v[2] * v[35])) - (v[3] * v[36])) - (v[4] * v[37]))) / 2),
				(((gammaUxx * v[17]) + (2 * gammaUxy * v[18]) + (2 * gammaUxz * v[19]) + (gammaUyy * v[20]) + (2 * gammaUyz * v[21]) + ((((gammaUzz * v[22]) - (v[3] * v[35])) - (v[5] * v[36])) - (v[6] * v[37]))) / 2),
				(((gammaUxx * v[23]) + (2 * gammaUxy * v[24]) + (2 * gammaUxz * v[25]) + (gammaUyy * v[26]) + (2 * gammaUyz * v[27]) + ((((gammaUzz * v[28]) - (v[4] * v[35])) - (v[6] * v[36])) - (v[7] * v[37]))) / 2),
				((((((v[8] - (f * gammaUxx * v[11])) - (2 * f * gammaUxy * v[12])) - (2 * f * gammaUxz * v[13])) - (f * gammaUyy * v[14])) - (2 * f * gammaUyz * v[15])) - (f * gammaUzz * v[16])),
				((v[9] + (2 * gammaUxx * v[12]) + (2 * math.sqrt(gammaUxx) * v[30]) + (2 * gammaUxy * v[18]) + ((((2 * gammaUxz * v[24]) - (v[3] * v[35])) - (v[5] * v[36])) - (v[6] * v[37]))) / 2),
				((v[10] + (2 * gammaUxx * v[13]) + (2 * math.sqrt(gammaUxx) * v[31]) + (2 * gammaUxy * v[19]) + ((((2 * gammaUxz * v[25]) - (v[4] * v[35])) - (v[6] * v[36])) - (v[7] * v[37]))) / 2),
				((gammaUxx * v[14]) + (math.sqrt(gammaUxx) * v[32]) + (gammaUxy * v[20]) + (gammaUxz * v[26])),
				((gammaUxx * v[15]) + (math.sqrt(gammaUxx) * v[33]) + (gammaUxy * v[21]) + (gammaUxz * v[27])),
				((gammaUxx * v[16]) + (math.sqrt(gammaUxx) * v[34]) + (gammaUxy * v[22]) + (gammaUxz * v[28])),
				(-(((((((((((((((((((((((((((v[35] - (math.sqrt(f) * (gammaUxx ^ (3 / 2)) * v[29])) - (2 * math.sqrt(f) * math.sqrt(gammaUxx) * gammaUxy * v[30])) - (2 * math.sqrt(f) * math.sqrt(gammaUxx) * gammaUxz * v[31])) - (math.sqrt(f) * math.sqrt(gammaUxx) * gammaUyy * v[32])) - (2 * math.sqrt(f) * math.sqrt(gammaUxx) * gammaUyz * v[33])) - (math.sqrt(f) * math.sqrt(gammaUxx) * gammaUzz * v[34])) - (gammaUxx * v[8])) - ((gammaUxx ^ 2) * v[11])) - (2 * gammaUxx * gammaUxy * v[12])) - (gammaUxx * gammaUxy * v[17])) - (2 * gammaUxx * gammaUxz * v[13])) - (gammaUxx * gammaUxz * v[23])) - (gammaUxy * v[9])) - (2 * (gammaUxy ^ 2) * v[18])) - (2 * gammaUxy * gammaUxz * v[19])) - (gammaUxz * v[10])) - (2 * (gammaUxz ^ 2) * v[25])) - (2 * gammaUxz * gammaUxy * v[24])) - (gammaUyy * gammaUxx * v[14])) - (gammaUyy * gammaUxy * v[20])) - (gammaUyy * gammaUxz * v[26])) - (2 * gammaUyz * gammaUxx * v[15])) - (2 * gammaUyz * gammaUxy * v[21])) - (2 * gammaUyz * gammaUxz * v[27])) - (gammaUzz * gammaUxx * v[16])) - (gammaUzz * gammaUxy * v[22])) - (gammaUzz * gammaUxz * v[28])))
			}
		end
	end
end

function ADM3D:calcEigenBasis(eigenvalues, rightEigenvectors, leftEigenvectors, fluxMatrix, ...)
	fill(eigenvalues, self:calcEigenvaluesFromCons(...))
	fill(leftEigenvectors, ...)
	fill(rightEigenvectors, ...)
	if fluxMatrix then fill(fluxMatrix, ...) end
end

function ADM3D:calcEigenvaluesFromCons(
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
end

function ADM3D:get_V_from_state(q, gammaUxx, gammaUxy, gammaUxz, gammaUyy, gammaUyz, gammaUzz)
	if self.useMomentumConstraints then
		return unpack(q, 35, 37)
	end
		
	local Gamma_x, Gamma_y, Gamma_z
	if self.useContractedGammaLower then
		Gamma_x, Gamma_y, Gamma_z = unpack(q, 35, 37)
	else
		local gamma_xx, gamma_xy, gamma_xz, gamma_yy, gamma_yz, gamma_zz = unpack(q, 2, 7)
		local GammaUx, GammaUy, GammaUz = unpack(q, 35, 37)
		Gamma_x, Gamma_y, Gamma_z = mat33sym.mul(
			gamma_xx, gamma_xy, gamma_xz, gamma_yy, gamma_yz, gamma_zz,
			GammaUx, GammaUy, GammaUz)
	end
	
	local D = {{unpack(q, 11, 16)}, {unpack(q, 17, 22)}, {unpack(q, 23, 28)}}
	local gammaU = {gammaUxx, gammaUxy, gammaUxz, gammaUyy, gammaUyz, gammaUzz}
	local V_x, V_y, V_z = 0, 0, 0
	for i=1,3 do
		for j=1,3 do
			local gammaUij = sym3x3(gammaU, i, j)
			V_x = V_x + (sym3x3(D[1], i, j) * gammaUij - Gamma_x) / 2
			V_y = V_y + (sym3x3(D[2], i, j) * gammaUij - Gamma_y) / 2
			V_z = V_z + (sym3x3(D[3], i, j) * gammaUij - Gamma_z) / 2
		end
	end
	return V_x, V_y, V_z
end

function ADM3D:sourceTerm(solver, qs, dt)
	local source = solver:newState()
	for i=2,solver.gridsize-1 do
		local alpha = qs[i][1]
		local gamma_xx, gamma_xy, gamma_xz, gamma_yy, gamma_yz, gamma_zz = unpack(qs[i], 2, 7)
		local a_x, a_y, a_z = unpack(qs[i], 8, 10)
		local d_xxx, d_xxy, d_xxz, d_xyy, d_xyz, d_xzz = unpack(qs[i], 11, 16)
		local d_yxx, d_yxy, d_yxz, d_yyy, d_yyz, d_yzz = unpack(qs[i], 17, 22)
		local d_zxx, d_zxy, d_zxz, d_zyy, d_zyz, d_zzz = unpack(qs[i], 23, 28)
		local K_xx, K_xy, K_xz, K_yy, K_yz, K_zz = unpack(qs[i], 29, 34)
		local gammaUxx, gammaUxy, gammaUxz, gammaUyy, gammaUyz, gammaUzz = mat33sym.inv(gamma_xx, gamma_xy, gamma_xz, gamma_yy, gamma_yz, gamma_zz)
		local f = self.calc.f(alpha)
		
		local V_x, V_y, V_z = self:get_V_from_state(qs[i], gammaUxx, gammaUxy, gammaUxz, gammaUyy, gammaUyz, gammaUzz)

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
},};
local trK = KUL[1][1] + KUL[2][2] + KUL[3][3];
local KSqSymLL = {
K_xx * KUL[1][1] + K_xy * KUL[2][1] + K_xz * KUL[3][1],
K_xx * KUL[1][2] + K_xy * KUL[2][2] + K_xz * KUL[3][2],
K_xx * KUL[1][3] + K_xy * KUL[2][3] + K_xz * KUL[3][3],
K_xy * KUL[1][2] + K_yy * KUL[2][2] + K_yz * KUL[3][2],
K_xy * KUL[1][3] + K_yy * KUL[2][3] + K_yz * KUL[3][3],
K_xz * KUL[1][3] + K_yz * KUL[2][3] + K_zz * KUL[3][3],
};
local DLUL = {
{{d_xxx * gammaUxx + d_xxy * gammaUxy + d_xxz * gammaUxz,
d_xxy * gammaUxx + d_xyy * gammaUxy + d_xyz * gammaUxz,
d_xxz * gammaUxx + d_xyz * gammaUxy + d_xzz * gammaUxz,
},{d_xxx * gammaUxy + d_xxy * gammaUyy + d_xxz * gammaUyz,
d_xxy * gammaUxy + d_xyy * gammaUyy + d_xyz * gammaUyz,
d_xxz * gammaUxy + d_xyz * gammaUyy + d_xzz * gammaUyz,
},{d_xxx * gammaUxz + d_xxy * gammaUyz + d_xxz * gammaUzz,
d_xxy * gammaUxz + d_xyy * gammaUyz + d_xyz * gammaUzz,
d_xxz * gammaUxz + d_xyz * gammaUyz + d_xzz * gammaUzz,
},},{{d_yxx * gammaUxx + d_yxy * gammaUxy + d_yxz * gammaUxz,
d_yxy * gammaUxx + d_yyy * gammaUxy + d_yyz * gammaUxz,
d_yxz * gammaUxx + d_yyz * gammaUxy + d_yzz * gammaUxz,
},{d_yxx * gammaUxy + d_yxy * gammaUyy + d_yxz * gammaUyz,
d_yxy * gammaUxy + d_yyy * gammaUyy + d_yyz * gammaUyz,
d_yxz * gammaUxy + d_yyz * gammaUyy + d_yzz * gammaUyz,
},{d_yxx * gammaUxz + d_yxy * gammaUyz + d_yxz * gammaUzz,
d_yxy * gammaUxz + d_yyy * gammaUyz + d_yyz * gammaUzz,
d_yxz * gammaUxz + d_yyz * gammaUyz + d_yzz * gammaUzz,
},},{{d_zxx * gammaUxx + d_zxy * gammaUxy + d_zxz * gammaUxz,
d_zxy * gammaUxx + d_zyy * gammaUxy + d_zyz * gammaUxz,
d_zxz * gammaUxx + d_zyz * gammaUxy + d_zzz * gammaUxz,
},{d_zxx * gammaUxy + d_zxy * gammaUyy + d_zxz * gammaUyz,
d_zxy * gammaUxy + d_zyy * gammaUyy + d_zyz * gammaUyz,
d_zxz * gammaUxy + d_zyz * gammaUyy + d_zzz * gammaUyz,
},{d_zxx * gammaUxz + d_zxy * gammaUyz + d_zxz * gammaUzz,
d_zxy * gammaUxz + d_zyy * gammaUyz + d_zyz * gammaUzz,
d_zxz * gammaUxz + d_zyz * gammaUyz + d_zzz * gammaUzz,
},},};
local D1L = {
DLUL[1][1][1] + DLUL[1][2][2] + DLUL[1][3][3],
DLUL[2][1][1] + DLUL[2][2][2] + DLUL[2][3][3],
DLUL[3][1][1] + DLUL[3][2][2] + DLUL[3][3][3],
};
local D3L = {
DLUL[1][1][1] + DLUL[2][2][1] + DLUL[3][3][1],
DLUL[1][1][2] + DLUL[2][2][2] + DLUL[3][3][2],
DLUL[1][1][3] + DLUL[2][2][3] + DLUL[3][3][3],
};
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
},},};
local D12SymLL = {
d_xxx * DUUL[1][1][1] + d_xxy * DUUL[1][2][1] + d_xxz * DUUL[1][3][1] + d_yxx * DUUL[2][1][1] + d_yxy * DUUL[2][2][1] + d_yxz * DUUL[2][3][1] + d_zxx * DUUL[3][1][1] + d_zxy * DUUL[3][2][1] + d_zxz * DUUL[3][3][1],
d_xxy * DUUL[1][1][1] + d_xyy * DUUL[1][2][1] + d_xyz * DUUL[1][3][1] + d_yxy * DUUL[2][1][1] + d_yyy * DUUL[2][2][1] + d_yyz * DUUL[2][3][1] + d_zxy * DUUL[3][1][1] + d_zyy * DUUL[3][2][1] + d_zyz * DUUL[3][3][1],
d_xxz * DUUL[1][1][1] + d_xyz * DUUL[1][2][1] + d_xzz * DUUL[1][3][1] + d_yxz * DUUL[2][1][1] + d_yyz * DUUL[2][2][1] + d_yzz * DUUL[2][3][1] + d_zxz * DUUL[3][1][1] + d_zyz * DUUL[3][2][1] + d_zzz * DUUL[3][3][1],
d_xxy * DUUL[1][1][2] + d_xyy * DUUL[1][2][2] + d_xyz * DUUL[1][3][2] + d_yxy * DUUL[2][1][2] + d_yyy * DUUL[2][2][2] + d_yyz * DUUL[2][3][2] + d_zxy * DUUL[3][1][2] + d_zyy * DUUL[3][2][2] + d_zyz * DUUL[3][3][2],
d_xxz * DUUL[1][1][2] + d_xyz * DUUL[1][2][2] + d_xzz * DUUL[1][3][2] + d_yxz * DUUL[2][1][2] + d_yyz * DUUL[2][2][2] + d_yzz * DUUL[2][3][2] + d_zxz * DUUL[3][1][2] + d_zyz * DUUL[3][2][2] + d_zzz * DUUL[3][3][2],
d_xxz * DUUL[1][1][3] + d_xyz * DUUL[1][2][3] + d_xzz * DUUL[1][3][3] + d_yxz * DUUL[2][1][3] + d_yyz * DUUL[2][2][3] + d_yzz * DUUL[2][3][3] + d_zxz * DUUL[3][1][3] + d_zyz * DUUL[3][2][3] + d_zzz * DUUL[3][3][3],
};
local GammaLSymLL = {
{d_xxx,
d_yxx,
d_zxx,
((2 * d_yxy) - d_xyy),
(d_zxy + (d_yxz - d_xyz)),
((2 * d_zxz) - d_xzz),
},{((2 * d_xxy) - d_yxx),
d_xyy,
(d_zxy + (d_xyz - d_yxz)),
d_yyy,
d_zyy,
((2 * d_zyz) - d_yzz),
},{((2 * d_xxz) - d_zxx),
(d_yxz + (d_xyz - d_zxy)),
d_xzz,
((2 * d_yyz) - d_zyy),
d_yzz,
d_zzz,
},};
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
},};
local Gamma3L = {
GammaUSymLL[1][1] + GammaUSymLL[2][2] + GammaUSymLL[3][3],
GammaUSymLL[1][2] + GammaUSymLL[2][4] + GammaUSymLL[3][5],
GammaUSymLL[1][3] + GammaUSymLL[2][5] + GammaUSymLL[3][6],
};
local Gamma31SymLL = {
Gamma3L[1] * GammaUSymLL[1][1] + Gamma3L[2] * GammaUSymLL[2][1] + Gamma3L[3] * GammaUSymLL[3][1],
Gamma3L[1] * GammaUSymLL[1][2] + Gamma3L[2] * GammaUSymLL[2][2] + Gamma3L[3] * GammaUSymLL[3][2],
Gamma3L[1] * GammaUSymLL[1][3] + Gamma3L[2] * GammaUSymLL[2][3] + Gamma3L[3] * GammaUSymLL[3][3],
Gamma3L[1] * GammaUSymLL[1][4] + Gamma3L[2] * GammaUSymLL[2][4] + Gamma3L[3] * GammaUSymLL[3][4],
Gamma3L[1] * GammaUSymLL[1][5] + Gamma3L[2] * GammaUSymLL[2][5] + Gamma3L[3] * GammaUSymLL[3][5],
Gamma3L[1] * GammaUSymLL[1][6] + Gamma3L[2] * GammaUSymLL[2][6] + Gamma3L[3] * GammaUSymLL[3][6],
};
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
},},};
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
},};
local Gamma11SymLL = {
GammaLSymLL[1][1] * GammaLSymUU[1][1] + GammaLSymLL[1][2] * GammaLSymUU[1][2] + GammaLSymLL[1][3] * GammaLSymUU[1][3] + GammaLSymLL[1][2] * GammaLSymUU[1][2] + GammaLSymLL[1][4] * GammaLSymUU[1][4] + GammaLSymLL[1][5] * GammaLSymUU[1][5] + GammaLSymLL[1][3] * GammaLSymUU[1][3] + GammaLSymLL[1][5] * GammaLSymUU[1][5] + GammaLSymLL[1][6] * GammaLSymUU[1][6],
GammaLSymLL[1][1] * GammaLSymUU[2][1] + GammaLSymLL[1][2] * GammaLSymUU[2][2] + GammaLSymLL[1][3] * GammaLSymUU[2][3] + GammaLSymLL[1][2] * GammaLSymUU[2][2] + GammaLSymLL[1][4] * GammaLSymUU[2][4] + GammaLSymLL[1][5] * GammaLSymUU[2][5] + GammaLSymLL[1][3] * GammaLSymUU[2][3] + GammaLSymLL[1][5] * GammaLSymUU[2][5] + GammaLSymLL[1][6] * GammaLSymUU[2][6],
GammaLSymLL[1][1] * GammaLSymUU[3][1] + GammaLSymLL[1][2] * GammaLSymUU[3][2] + GammaLSymLL[1][3] * GammaLSymUU[3][3] + GammaLSymLL[1][2] * GammaLSymUU[3][2] + GammaLSymLL[1][4] * GammaLSymUU[3][4] + GammaLSymLL[1][5] * GammaLSymUU[3][5] + GammaLSymLL[1][3] * GammaLSymUU[3][3] + GammaLSymLL[1][5] * GammaLSymUU[3][5] + GammaLSymLL[1][6] * GammaLSymUU[3][6],
GammaLSymLL[2][1] * GammaLSymUU[2][1] + GammaLSymLL[2][2] * GammaLSymUU[2][2] + GammaLSymLL[2][3] * GammaLSymUU[2][3] + GammaLSymLL[2][2] * GammaLSymUU[2][2] + GammaLSymLL[2][4] * GammaLSymUU[2][4] + GammaLSymLL[2][5] * GammaLSymUU[2][5] + GammaLSymLL[2][3] * GammaLSymUU[2][3] + GammaLSymLL[2][5] * GammaLSymUU[2][5] + GammaLSymLL[2][6] * GammaLSymUU[2][6],
GammaLSymLL[2][1] * GammaLSymUU[3][1] + GammaLSymLL[2][2] * GammaLSymUU[3][2] + GammaLSymLL[2][3] * GammaLSymUU[3][3] + GammaLSymLL[2][2] * GammaLSymUU[3][2] + GammaLSymLL[2][4] * GammaLSymUU[3][4] + GammaLSymLL[2][5] * GammaLSymUU[3][5] + GammaLSymLL[2][3] * GammaLSymUU[3][3] + GammaLSymLL[2][5] * GammaLSymUU[3][5] + GammaLSymLL[2][6] * GammaLSymUU[3][6],
GammaLSymLL[3][1] * GammaLSymUU[3][1] + GammaLSymLL[3][2] * GammaLSymUU[3][2] + GammaLSymLL[3][3] * GammaLSymUU[3][3] + GammaLSymLL[3][2] * GammaLSymUU[3][2] + GammaLSymLL[3][4] * GammaLSymUU[3][4] + GammaLSymLL[3][5] * GammaLSymUU[3][5] + GammaLSymLL[3][3] * GammaLSymUU[3][3] + GammaLSymLL[3][5] * GammaLSymUU[3][5] + GammaLSymLL[3][6] * GammaLSymUU[3][6],
};
local ADL = {
a_x - 2 * D3L[1],
a_y - 2 * D3L[2],
a_z - 2 * D3L[3],
};
local ADU = {
gammaUxx * ADL[1] + gammaUxy * ADL[2] + gammaUxz * ADL[3],
gammaUxy * ADL[1] + gammaUyy * ADL[2] + gammaUyz * ADL[3],
gammaUxz * ADL[1] + gammaUyz * ADL[2] + gammaUzz * ADL[3],
};
local ADDSymLL = {
ADU[1] * (2 * d_xxx) + ADU[2] * (2 * d_xxy) + ADU[3] * (2 * d_xxz),
ADU[1] * (d_xxy + d_yxx) + ADU[2] * (d_xyy + d_yxy) + ADU[3] * (d_xyz + d_yxz),
ADU[1] * (d_xxz + d_zxx) + ADU[2] * (d_xyz + d_zxy) + ADU[3] * (d_xzz + d_zxz),
ADU[1] * (2 * d_yxy) + ADU[2] * (2 * d_yyy) + ADU[3] * (2 * d_yyz),
ADU[1] * (d_yxz + d_zxy) + ADU[2] * (d_yyz + d_zyy) + ADU[3] * (d_yzz + d_zyz),
ADU[1] * (2 * d_zxz) + ADU[2] * (2 * d_zyz) + ADU[3] * (2 * d_zzz),
};
local R4SymLL = {
0,
0,
0,
0,
0,
0,
};
local SSymLL = {
-R4SymLL[1] + trK * K_xx - 2 * KSqSymLL[1] + 4 * D12SymLL[1] + Gamma31SymLL[1] - Gamma11SymLL[1] + ADDSymLL[1] + (a_x * ((2 * V_x) - D1L[1])),
-R4SymLL[2] + trK * K_xy - 2 * KSqSymLL[2] + 4 * D12SymLL[2] + Gamma31SymLL[2] - Gamma11SymLL[2] + ADDSymLL[2] + ((((2 * a_y * V_x) - (a_y * D1L[1])) + ((2 * a_x * V_y) - (a_x * D1L[2]))) / 2),
-R4SymLL[3] + trK * K_xz - 2 * KSqSymLL[3] + 4 * D12SymLL[3] + Gamma31SymLL[3] - Gamma11SymLL[3] + ADDSymLL[3] + ((((2 * a_z * V_x) - (a_z * D1L[1])) + ((2 * a_x * V_z) - (a_x * D1L[3]))) / 2),
-R4SymLL[4] + trK * K_yy - 2 * KSqSymLL[4] + 4 * D12SymLL[4] + Gamma31SymLL[4] - Gamma11SymLL[4] + ADDSymLL[4] + (a_y * ((2 * V_y) - D1L[2])),
-R4SymLL[5] + trK * K_yz - 2 * KSqSymLL[5] + 4 * D12SymLL[5] + Gamma31SymLL[5] - Gamma11SymLL[5] + ADDSymLL[5] + ((((2 * a_z * V_y) - (a_z * D1L[2])) + ((2 * a_y * V_z) - (a_y * D1L[3]))) / 2),
-R4SymLL[6] + trK * K_zz - 2 * KSqSymLL[6] + 4 * D12SymLL[6] + Gamma31SymLL[6] - Gamma11SymLL[6] + ADDSymLL[6] + (a_z * ((2 * V_z) - D1L[3])),
};
local GU0L = {
0,
0,
0,
};
local AKL = {
a_x * KUL[1][1] + a_y * KUL[2][1] + a_z * KUL[3][1],
a_x * KUL[1][2] + a_y * KUL[2][2] + a_z * KUL[3][2],
a_x * KUL[1][3] + a_y * KUL[2][3] + a_z * KUL[3][3],
};
local K12D23L = {
KUL[1][1] * DLUL[1][1][1] +KUL[1][2] * DLUL[1][2][1] +KUL[1][3] * DLUL[1][3][1] + KUL[2][1] * DLUL[1][1][2] +KUL[2][2] * DLUL[1][2][2] +KUL[2][3] * DLUL[1][3][2] + KUL[3][1] * DLUL[1][1][3] +KUL[3][2] * DLUL[1][2][3] +KUL[3][3] * DLUL[1][3][3],
KUL[1][1] * DLUL[2][1][1] +KUL[1][2] * DLUL[2][2][1] +KUL[1][3] * DLUL[2][3][1] + KUL[2][1] * DLUL[2][1][2] +KUL[2][2] * DLUL[2][2][2] +KUL[2][3] * DLUL[2][3][2] + KUL[3][1] * DLUL[2][1][3] +KUL[3][2] * DLUL[2][2][3] +KUL[3][3] * DLUL[2][3][3],
KUL[1][1] * DLUL[3][1][1] +KUL[1][2] * DLUL[3][2][1] +KUL[1][3] * DLUL[3][3][1] + KUL[2][1] * DLUL[3][1][2] +KUL[2][2] * DLUL[3][2][2] +KUL[2][3] * DLUL[3][3][2] + KUL[3][1] * DLUL[3][1][3] +KUL[3][2] * DLUL[3][2][3] +KUL[3][3] * DLUL[3][3][3],
};
local KD23L = {
KUL[1][1] * D1L[1] + KUL[2][1] * D1L[2] + KUL[3][1] * D1L[3],
KUL[1][2] * D1L[1] + KUL[2][2] * D1L[2] + KUL[3][2] * D1L[3],
KUL[1][3] * D1L[1] + KUL[2][3] * D1L[2] + KUL[3][3] * D1L[3],
};
local K12D12L = {
KUL[1][1] * DLUL[1][1][1] + KUL[1][2] * DLUL[1][2][1] + KUL[1][3] * DLUL[1][3][1] + KUL[2][1] * DLUL[2][1][1] + KUL[2][2] * DLUL[2][2][1] + KUL[2][3] * DLUL[2][3][1] + KUL[3][1] * DLUL[3][1][1] + KUL[3][2] * DLUL[3][2][1] + KUL[3][3] * DLUL[3][3][1],
KUL[1][1] * DLUL[1][1][2] + KUL[1][2] * DLUL[1][2][2] + KUL[1][3] * DLUL[1][3][2] + KUL[2][1] * DLUL[2][1][2] + KUL[2][2] * DLUL[2][2][2] + KUL[2][3] * DLUL[2][3][2] + KUL[3][1] * DLUL[3][1][2] + KUL[3][2] * DLUL[3][2][2] + KUL[3][3] * DLUL[3][3][2],
KUL[1][1] * DLUL[1][1][3] + KUL[1][2] * DLUL[1][2][3] + KUL[1][3] * DLUL[1][3][3] + KUL[2][1] * DLUL[2][1][3] + KUL[2][2] * DLUL[2][2][3] + KUL[2][3] * DLUL[2][3][3] + KUL[3][1] * DLUL[3][1][3] + KUL[3][2] * DLUL[3][2][3] + KUL[3][3] * DLUL[3][3][3],
};
local KD12L = {
KUL[1][1] * D3L[1] + KUL[2][1] * D3L[2] + KUL[3][1] * D3L[3],
KUL[1][2] * D3L[1] + KUL[2][2] * D3L[2] + KUL[3][2] * D3L[3],
KUL[1][3] * D3L[1] + KUL[2][3] * D3L[2] + KUL[3][3] * D3L[3],
};
local PL = {
GU0L[1] + AKL[1] - a_x * trK + K12D23L[1] + KD23L[1] - 2 * K12D12L[1] + 2 * KD12L[1],
GU0L[2] + AKL[2] - a_y * trK + K12D23L[2] + KD23L[2] - 2 * K12D12L[2] + 2 * KD12L[2],
GU0L[3] + AKL[3] - a_z * trK + K12D23L[3] + KD23L[3] - 2 * K12D12L[3] + 2 * KD12L[3],
};


		source[i][1] = -alpha * alpha * f * trK
		source[i][2] = -2 * alpha * K_xx
		source[i][3] = -2 * alpha * K_xy
		source[i][4] = -2 * alpha * K_xz
		source[i][5] = -2 * alpha * K_yy
		source[i][6] = -2 * alpha * K_yz
		source[i][7] = -2 * alpha * K_zz
		source[i][29] = alpha * SSymLL[1]
		source[i][30] = alpha * SSymLL[2]
		source[i][31] = alpha * SSymLL[3]
		source[i][32] = alpha * SSymLL[4]
		source[i][33] = alpha * SSymLL[5]
		source[i][34] = alpha * SSymLL[6]
		
		if self.useMomentumConstraints then
			source[i][35] = alpha * PL[1]
			source[i][36] = alpha * PL[2]
			source[i][37] = alpha * PL[3]
		else
			-- source terms of Gamma^i ?
		end
	
		-- [[ converge to constraints
		local eta = 1/dt	-- damping term / constraint enforcing of 1st order terms
		local _2dx = solver.xs[i+1] - solver.xs[i-1]
		-- a_i = alpha,i / alpha <=> a_i += eta (alpha,i / alpha - a_i)
		local dx_alpha = (solver.qs[i+1][1] - solver.qs[i-1][1]) / _2dx
		source[i][8] = source[i][8] + eta * (dx_alpha / alpha - a_x)
		-- alpha,y and alpha,z = 0
		source[i][9] = source[i][9] + eta * (-a_y)
		source[i][10] = source[i][10] + eta * (-a_z)
		
		for ij=1,6 do
			local gammaIndex = ij + 1
			local dIndex = ij-1 + 11
			-- d_xij = 1/2 gamma_ij,x <=> d_xij += eta (1/2 gamma_ij,x - d_kij)
			local dx_gamma_ij = (solver.qs[i+1][gammaIndex] - solver.qs[i-1][gammaIndex]) / _2dx
			source[i][dIndex] = source[i][dIndex] + eta * (.5 * dx_gamma_ij - solver.qs[i][dIndex])
			-- d_yij = d_zij = 0
			source[i][dIndex+6] = source[i][dIndex+6] + eta * (-solver.qs[i][dIndex+6])
			source[i][dIndex+12] = source[i][dIndex+12] + eta * (-solver.qs[i][dIndex+12])
		end
		
		-- V_i = d_ij^j - d^j_ji <=> V_i += eta (d_ij^j - d^j_ji - V_i)
		

		--]]
	end
	return source
end

--[[ enforce constraint V_i = d_im^m - d^m_mi
function ADM3D:postIterate(solver, qs)
	for i=1,solver.gridsize do
		local gamma_xx, gamma_xy, gamma_xz, gamma_yy, gamma_yz, gamma_zz = unpack(qs[i], 2, 7)
		local d_xxx, d_xxy, d_xxz, d_xyy, d_xyz, d_xzz = unpack(qs[i], 11, 16)
		local d_yxx, d_yxy, d_yxz, d_yyy, d_yyz, d_yzz = unpack(qs[i], 17, 22)
		local d_zxx, d_zxy, d_zxz, d_zyy, d_zyz, d_zzz = unpack(qs[i], 23, 28)
		local gammaUxx, gammaUxy, gammaUxz, gammaUyy, gammaUyz, gammaUzz = mat33sym.inv(gamma_xx, gamma_xy, gamma_xz, gamma_yy, gamma_yz, gamma_zz)

		local V_x, V_y, V_z = self:get_V_from_state(qs[i], gammaUxx, gammaUxy, gammaUxz, gammaUyy, gammaUyz, gammaUzz)

		if not self.useMomentumConstraints then
			-- Gamma_i = 2 d^m_mi - d_im^m
			-- but evolution should preserve this ... ?
		else
			-- direct assign (seems like this would be constantly overwriting the V_k source term contribution
			if self.momentumConstraintMethod == 'directAssign' then
				qs[i][35] = ((gammaUxy * d_xxy) + (gammaUxz * d_xxz) + (gammaUyy * d_xyy) + (2 * gammaUyz * d_xyz) + (((((((gammaUzz * d_xzz) - (gammaUxy * d_yxx)) - (gammaUxz * d_zxx)) - (gammaUyy * d_yxy)) - (gammaUyz * d_zxy)) - (gammaUyz * d_yxz)) - (gammaUzz * d_zxz)))
				qs[i][36] = ((gammaUxx * d_yxx) + (gammaUxy * d_yxy) + (2 * gammaUxz * d_yxz) + (gammaUyz * d_yyz) + (((((((gammaUzz * d_yzz) - (gammaUxx * d_xxy)) - (gammaUxz * d_zxy)) - (gammaUxy * d_xyy)) - (gammaUyz * d_zyy)) - (gammaUxz * d_xyz)) - (gammaUzz * d_zyz)))
				qs[i][37] = ((gammaUxx * d_zxx) + (2 * gammaUxy * d_zxy) + (gammaUxz * d_zxz) + (gammaUyy * d_zyy) + (((((((gammaUyz * d_zyz) - (gammaUxx * d_xxz)) - (gammaUxy * d_yxz)) - (gammaUxy * d_xyz)) - (gammaUyy * d_yyz)) - (gammaUxz * d_xzz)) - (gammaUyz * d_yzz)))
			
			-- linear projection ... would work fine if the D's weren't intermixed ... but because they are, this is 3 separate linear projections with intermixed terms ...
			elseif self.momentumConstraintMethod == 'linearProject' then

	-- x
local aDotA = (1 + (2 * (gammaUxy ^ 2)) + (2 * (gammaUxz ^ 2)) + (2 * (gammaUyy ^ 2)) + (6 * (gammaUyz ^ 2)) + (2 * (gammaUzz ^ 2)))
local vDotA = ((((((-(gammaUxy * d_xxy)) - (gammaUxz * d_xxz)) - (gammaUyy * d_xyy)) - (2 * gammaUyz * d_xyz)) - (gammaUzz * d_xzz)) + (gammaUxy * d_yxx) + (gammaUyy * d_yxy) + (gammaUyz * d_yxz) + (gammaUxz * d_zxx) + (gammaUyz * d_zxy) + (gammaUzz * d_zxz) + V_x)
local v_a = vDotA / aDotA
local epsilon = 1/100
qs[i][12] = qs[i][12] + (epsilon * v_a * gammaUxy)
qs[i][13] = qs[i][13] + (epsilon * v_a * gammaUxz)
qs[i][14] = qs[i][14] + (epsilon * v_a * gammaUyy)
qs[i][15] = qs[i][15] + (2 * epsilon * v_a * gammaUyz)
qs[i][16] = qs[i][16] + (epsilon * v_a * gammaUzz)
qs[i][17] = qs[i][17] + (-(epsilon * v_a * gammaUxy))
qs[i][18] = qs[i][18] + (-(epsilon * v_a * gammaUyy))
qs[i][19] = qs[i][19] + (-(epsilon * v_a * gammaUyz))
qs[i][23] = qs[i][23] + (-(epsilon * v_a * gammaUxz))
qs[i][24] = qs[i][24] + (-(epsilon * v_a * gammaUyz))
qs[i][25] = qs[i][25] + (-(epsilon * v_a * gammaUzz))
qs[i][35] = qs[i][35] + (-(epsilon * v_a))
	-- y
local aDotA = (1 + (2 * (gammaUxx ^ 2)) + (2 * (gammaUxy ^ 2)) + (6 * (gammaUxz ^ 2)) + (2 * (gammaUyz ^ 2)) + (2 * (gammaUzz ^ 2)))
local vDotA = ((gammaUxx * d_xxy) + (gammaUxy * d_xyy) + ((((((gammaUxz * d_xyz) - (gammaUxx * d_yxx)) - (gammaUxy * d_yxy)) - (2 * gammaUxz * d_yxz)) - (gammaUyz * d_yyz)) - (gammaUzz * d_yzz)) + (gammaUxz * d_zxy) + (gammaUyz * d_zyy) + (gammaUzz * d_zyz) + V_y)
local v_a = vDotA / aDotA
local epsilon = 1/100
qs[i][12] = qs[i][12] + (-(epsilon * v_a * gammaUxx))
qs[i][14] = qs[i][14] + (-(epsilon * v_a * gammaUxy))
qs[i][15] = qs[i][15] + (-(epsilon * v_a * gammaUxz))
qs[i][17] = qs[i][17] + (epsilon * v_a * gammaUxx)
qs[i][18] = qs[i][18] + (epsilon * v_a * gammaUxy)
qs[i][19] = qs[i][19] + (2 * epsilon * v_a * gammaUxz)
qs[i][21] = qs[i][21] + (epsilon * v_a * gammaUyz)
qs[i][22] = qs[i][22] + (epsilon * v_a * gammaUzz)
qs[i][24] = qs[i][24] + (-(epsilon * v_a * gammaUxz))
qs[i][26] = qs[i][26] + (-(epsilon * v_a * gammaUyz))
qs[i][27] = qs[i][27] + (-(epsilon * v_a * gammaUzz))
qs[i][36] = qs[i][36] + (-(epsilon * v_a))
	-- z
local aDotA = (1 + (2 * (gammaUxx ^ 2)) + (6 * (gammaUxy ^ 2)) + (2 * (gammaUxz ^ 2)) + (2 * (gammaUyy ^ 2)) + (2 * (gammaUyz ^ 2)))
local vDotA = ((gammaUxx * d_xxz) + (gammaUxy * d_xyz) + (gammaUxz * d_xzz) + (gammaUxy * d_yxz) + (gammaUyy * d_yyz) + ((((((gammaUyz * d_yzz) - (gammaUxx * d_zxx)) - (2 * gammaUxy * d_zxy)) - (gammaUxz * d_zxz)) - (gammaUyy * d_zyy)) - (gammaUyz * d_zyz)) + V_z)
local v_a = vDotA / aDotA
local epsilon = 1/100
qs[i][13] = qs[i][13] + (-(epsilon * v_a * gammaUxx))
qs[i][15] = qs[i][15] + (-(epsilon * v_a * gammaUxy))
qs[i][16] = qs[i][16] + (-(epsilon * v_a * gammaUxz))
qs[i][19] = qs[i][19] + (-(epsilon * v_a * gammaUxy))
qs[i][21] = qs[i][21] + (-(epsilon * v_a * gammaUyy))
qs[i][22] = qs[i][22] + (-(epsilon * v_a * gammaUyz))
qs[i][23] = qs[i][23] + (epsilon * v_a * gammaUxx)
qs[i][24] = qs[i][24] + (2 * epsilon * v_a * gammaUxy)
qs[i][25] = qs[i][25] + (epsilon * v_a * gammaUxz)
qs[i][26] = qs[i][26] + (epsilon * v_a * gammaUyy)
qs[i][27] = qs[i][27] + (epsilon * v_a * gammaUyz)
qs[i][37] = qs[i][37] + (-(epsilon * v_a))

			end
		end
	end
end
--]]

return ADM3D
