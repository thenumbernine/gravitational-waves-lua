local class = require 'ext.class'
local Equation = require 'equation'
local symmath = require 'symmath'

local BSSNOK3D = class(Equation)
BSSNOK3D.name = 'BSSNOK 3D Hyperbolic'
BSSNOK3D.numStates = 30

function BSSNOK3D:init(args, ...)
	local function makesym(field)
		local v = args[field]
		if not v then return end
		return symmath.clone(v)() 
	end
	
	local x = assert(args.x)
	local y = args.y or symmath.var'y'
	local z = args.z or symmath.var'z'
	local vars = table{x,y,z}
	local symvars = table{'xx','xy','xz','yy','yz','zz'}
	local eta_ij = table{1,0,0,1,0,1}
	
	local gamma = symvars:map(function(ij_field, ij_index)
		local gamma_ij_field = 'gamma_'..ij_field
		local gamma_ij_expr = makesym(gamma_ij_field)
		if ij_field == 'xx' then assert(gamma_ij_expr) end
		return gamma_ij_expr or symmath.Constant(eta_ij[ij_index])
	end)

	error'TODO'
end
