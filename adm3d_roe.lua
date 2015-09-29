local class = require 'ext.class'
local table = require 'ext.table'
local Roe = require 'roe'
local ADM3D = require 'adm3d'
local mat33 = require 'mat33'

local ADM3DRoe = class(Roe)
ADM3DRoe.name = 'ADM 3D Roe'

function ADM3DRoe:init(args)
	args = table(args)
	args.equation = ADM3D(args)
	ADM3DRoe.super.init(self, args)
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
