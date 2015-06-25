local class = require 'ext.class'
local table = require 'ext.table'

local Roe = require 'roe'
local mat33 = require 'mat33'

local ADM3DRoe = class(Roe)
ADM3DRoe.name = 'ADM 3D Roe'

function ADM3DRoe:init(args)
	args = table(args)
	args.equation = require 'adm3d'(args)
	ADM3DRoe.super.init(self, args)
end

function ADM3DRoe:fluxTransform(i, v)
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
		0,0,0			-- V_k
	}
end

function ADM3DRoe:eigenfields(i, v)

	-- interface eigenfield varialbes
	local avgQ = {}
	for j=1,self.numStates do 
		avgQ[j] = (self.qs[i-1][j] + self.qs[i][j]) / 2
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
	local f = self.equation.calc.f(alpha)

	-- cell variables
	-- what if, for the ADM equations, there is no distinction?
	-- they're used for Roe's scheme for computing deltas in eigenbasis coordinates by which to scale coordinates coinciding with the lambdas ...
	-- what about creating them solely from 'v' rather than using the average whatsoever?
	-- this would mean ensuring the inputs to the eigenfields() functions were always the state variables themselves (not differences or averages)
	-- 	and deferring differences or averages til after eigenfields() is called (assuming it is a linear function)
	-- this also has an issue with eigenfieldsInverse(), which is called on a flux vector, i.e. at cell interface, which would probably need the average of cells for that input

	return {
		((((-(2 * gUxz * v[37])) - (gUxx * v[8])) + (math.sqrt(f) * (gUxx ^ (3 / 2)) * v[29]) + (math.sqrt(f) * gUxy * v[30] * math.sqrt(gUxx)) + (math.sqrt(f) * gUxz * v[31] * math.sqrt(gUxx)) + (math.sqrt(f) * gUyy * v[32] * math.sqrt(gUxx)) + (math.sqrt(f) * gUyz * v[33] * math.sqrt(gUxx)) + (((math.sqrt(f) * gUzz * v[34] * math.sqrt(gUxx)) - (2 * gUxx * v[35])) - (2 * gUxy * v[36]))) / math.sqrt(gUxx)),
		(((-((gUxx ^ (3 / 2)) * v[12])) + ((v[30] * gUxx) - (v[36]))) / gUxx),
		(((-((gUxx ^ (3 / 2)) * v[13])) + ((v[31] * gUxx) - (v[37]))) / gUxx),
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
		(((((((v[8] - (f * gUxx * v[11])) - (f * gUxy * v[12])) - (f * gUxz * v[13])) - (f * gUyy * v[14])) - (f * gUyz * v[15])) - (f * gUzz * v[16]))),
		((((gUxx ^ (3 / 2)) * v[12]) + (v[30] * gUxx) + v[36]) / gUxx),
		((((gUxx ^ (3 / 2)) * v[13]) + (v[31] * gUxx) + v[37]) / gUxx),
		((math.sqrt(gUxx) * v[14]) + v[32]),
		((math.sqrt(gUxx) * v[15]) + v[33]),
		((math.sqrt(gUxx) * v[16]) + v[34]),
		((((gUxx ^ (3 / 2)) * v[8]) + (math.sqrt(f) * (gUxx ^ 2) * v[29]) + (math.sqrt(f) * gUxy * v[30] * gUxx) + (math.sqrt(f) * gUxz * v[31] * gUxx) + (math.sqrt(f) * gUyy * v[32] * gUxx) + (math.sqrt(f) * gUyz * v[33] * gUxx) + (math.sqrt(f) * gUzz * v[34] * gUxx) + (2 * v[35])) / gUxx)
	}
end

function ADM3DRoe:eigenfieldsInverse(i, v)
	
	-- interface eigenfield varialbes
	local avgQ = {}
	for j=1,self.numStates do 
		avgQ[j] = (self.qs[i-1][j] + self.qs[i][j]) / 2
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
	local f = self.equation.calc.f(alpha)

	return {
		v[7],
		v[8],
		v[9],
		v[10],
		v[11],
		v[12],
		v[13],
		(((-(v[37] * gUxx)) + (2 * gUxz * v[30] * math.sqrt(gUxx)) + (2 * gUxy * v[29] * math.sqrt(gUxx)) + (v[1] * gUxx) + (2 * v[28]) + (2 * (gUxx ^ (3 / 2)) * v[28])) / (-(2 * (gUxx ^ (3 / 2))))),
		v[14],
		v[15],
		(((-(v[37] * gUxx)) + (gUzz * v[36] * gUxx * f) + (gUyz * v[35] * gUxx * f) + (gUyy * v[34] * gUxx * f) + (gUxz * v[33] * gUxx * f) + (gUxy * v[32] * gUxx * f) + (2 * v[31] * (gUxx ^ (3 / 2))) + ((2 * gUxz * math.sqrt(gUxx) * v[30]) - (2 * gUxz * f * v[30])) + (((((((2 * gUxy * math.sqrt(gUxx) * v[29]) - (2 * gUxy * f * v[29])) - (gUzz * v[6] * f * gUxx)) - (gUyz * v[5] * f * gUxx)) - (gUyy * v[4] * f * gUxx)) - (gUxz * v[3] * f * gUxx)) - (gUxy * v[2] * f * gUxx)) + (v[1] * gUxx) + (2 * v[28]) + (2 * (gUxx ^ (3 / 2)) * v[28])) / (-(2 * (gUxx ^ (5 / 2)) * f))),
		(((-(v[32] * gUxx)) + (v[2] * gUxx) + (2 * v[29])) / (-(2 * (gUxx ^ (3 / 2))))),
		(((-(v[33] * gUxx)) + (v[3] * gUxx) + (2 * v[30])) / (-(2 * (gUxx ^ (3 / 2))))),
		(((-(v[34])) + v[4]) / (-(2 * math.sqrt(gUxx)))),
		(((-(v[35])) + v[5]) / (-(2 * math.sqrt(gUxx)))),
		(((-(v[36])) + v[6]) / (-(2 * math.sqrt(gUxx)))),
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
		((((((((v[37] * gUxx) - (gUzz * v[36] * gUxx * math.sqrt(f))) - (gUyz * v[35] * gUxx * math.sqrt(f))) - (gUyy * v[34] * gUxx * math.sqrt(f))) - (gUxz * v[33] * gUxx * math.sqrt(f))) - (gUxy * v[32] * gUxx * math.sqrt(f))) + (2 * gUxz * v[30] * math.sqrt(gUxx)) + ((((((2 * gUxy * v[29] * math.sqrt(gUxx)) - (gUzz * v[6] * math.sqrt(f) * gUxx)) - (gUyz * v[5] * math.sqrt(f) * gUxx)) - (gUyy * v[4] * math.sqrt(f) * gUxx)) - (gUxz * v[3] * math.sqrt(f) * gUxx)) - (gUxy * v[2] * math.sqrt(f) * gUxx)) + ((v[1] * gUxx) - (2 * v[28])) + (2 * (gUxx ^ (3 / 2)) * v[28])) / (2 * (gUxx ^ 2) * math.sqrt(f))),
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

function ADM3DRoe:addSourceToDerivCell(dq_dts, i)
	local alpha = self.qs[i][1]
	local g_xx, g_xy, g_xz, g_yy, g_yz, g_zz = unpack(self.qs[i], 2, 7)
	local A_x, A_y, A_z = unpack(self.qs[i], 8, 10)
	local D_xxx, D_xxy, D_xxz, D_xyy, D_xyz, D_xzz = unpack(self.qs[i], 11, 16)
	local D_yxx, D_yxy, D_yxz, D_yyy, D_yyz, D_yzz = unpack(self.qs[i], 17, 22)
	local D_zxx, D_zxy, D_zxz, D_zyy, D_zyz, D_zzz = unpack(self.qs[i], 23, 28)
	local K_xx, K_xy, K_xz, K_yy, K_yz, K_zz = unpack(self.qs[i], 29, 34)
	local V_x, V_y, V_z = unpack(self.qs[i], 35, 37)
	local gUxx, gUxy, gUxz, gUyy, gUyz, gUzz = mat33.inv(g_xx, g_xy, g_xz, g_yy, g_yz, g_zz)
	local f = self.equation.calc.f(alpha)

	local tr_K = 
		K_xx * gUxx
		+ K_yy * gUyy
		+ K_zz * gUzz
		+ 2 * K_xy * gUxy
		+ 2 * K_yz * gUyz
		+ 2 * K_xz * gUxz

	dq_dts[i][1] = dq_dts[i][1] - alpha * alpha * f * tr_K
	dq_dts[i][2] = dq_dts[i][2] - 2 * alpha * K_xx
	dq_dts[i][3] = dq_dts[i][3] - 2 * alpha * K_xy
	dq_dts[i][4] = dq_dts[i][4] - 2 * alpha * K_xz
	dq_dts[i][5] = dq_dts[i][5] - 2 * alpha * K_yy
	dq_dts[i][6] = dq_dts[i][6] - 2 * alpha * K_yz
	dq_dts[i][7] = dq_dts[i][7] - 2 * alpha * K_zz
	-- TODO dq_dts[i][29..34] == partial_t K_ij += alpha * S_ij
	-- partial_t V_k = alpha * P_k
end

-- enforce constraint V_k = (D_kmn - D_mnk) g^mn
function ADM3DRoe:iterate(...)
	ADM3DRoe.super.iterate(self, ...)
	for i=1,self.gridsize do
		-- [[ direct assign (seems like this would be constantly overwriting the V_k source term contribution
		local g_xx, g_xy, g_xz, g_yy, g_yz, g_zz = unpack(self.qs[i], 2, 7)
		local D_xxx, D_xxy, D_xxz, D_xyy, D_xyz, D_xzz = unpack(self.qs[i], 11, 16)
		local D_yxx, D_yxy, D_yxz, D_yyy, D_yyz, D_yzz = unpack(self.qs[i], 17, 22)
		local D_zxx, D_zxy, D_zxz, D_zyy, D_zyz, D_zzz = unpack(self.qs[i], 23, 28)
		local gUxx, gUxy, gUxz, gUyy, gUyz, gUzz = mat33.inv(g_xx, g_xy, g_xz, g_yy, g_yz, g_zz)
		self.qs[i][35] = 
			(D_xxy - D_yxx) * gUxy
			+ (D_xxz - D_zxx) * gUxz
			+ (D_xyy - D_yxy) * gUyy
			+ (D_xyz - D_yxz) * gUyz
			+ (D_xyz - D_zxy) * gUyz
			+ (D_xzz - D_zxz) * gUzz
		self.qs[i][36] = 
			(D_yxx - D_xxy) * gUxx
			+ (D_yxy - D_xyy) * gUxy
			+ (D_yxz - D_xyz) * gUxz
			+ (D_yxz - D_zxy) * gUxz
			+ (D_yyz - D_zyy) * gUyz
			+ (D_yzz - D_zyz) * gUzz
		self.qs[i][37] = 
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


return ADM3DRoe

