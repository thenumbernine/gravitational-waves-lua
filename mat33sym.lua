--[[
3x3 matrix routines
used by adm3d and adm3droe
--]]

local det = function(xx, xy, xz, yy, yz, zz)
	return xx * yy * zz
		+ xy * yz * xz
		+ xz * xy * yz
		- xz * yy * xz
		- yz * yz * xx
		- zz * xy * xy
end

local inv = function(xx, xy, xz, yy, yz, zz, d)
	if not d then d = det(xx, xy, xz, yy, yz, zz) end
	return 
		(yy * zz - yz * yz) / d,	-- xx
		(xz * yz - xy * zz) / d,	-- xy
		(xy * yz - xz * yy) / d,	-- xz
		(xx * zz - xz * xz) / d,	-- yy
		(xz * xy - xx * yz) / d,	-- yz
		(xx * yy - xy * xy) / d	-- zz
end

local mul = function(xx, xy, xz, yy, yz, zz, x, y, z)
	return 
		xx * x + xy * y + xz * z,
		xy * x + yy * y + yz * z,
		xz * x + yz * y + zz * z
end

return {
	det = det,
	inv = inv,
	mul = mul,
}
