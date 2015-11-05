local class = require 'ext.class'
local Equation = require 'equation'
local BSSNOK1DOriginal = class(Equation)

BSSNOK1DOriginal.name = 'BSSNOK 1D Original'

function BSSNOK1DOriginal:init(args, ...)
	local symmath = require 'symmath'
	local function makesym(field)
		return symmath.clone(assert(args[field], "expected to find field "..field)):simplify() 
	end
	
	local x = assert(args.x)

	local stateExprs = table{'alpha', 'g_xx', 'K_xx'}:map(function(name)
		return makesym(name), name
	end

	stateExprs.phi = -symmath.log(stateExprs.g_xx)/4
	stateExprs.K = (stateExprs.K_xx / stateExprs.g_xx):simplify()
	stateExprs.gammaTilde_xx = stateExprs.g_xx / stateExprs.phi^4
	stateExprs.ATilde_xx = symmath.exp(-4 * stateExprs.phi) * (stateExprs.K_xx - stateExprs.K/3 * stateExprs.g_xx)
	-- GammaTilde^x

	-- convert from symbolic functions to Lua functions
	self.calc = stateExprs:map(function(expr, name)
		return expr:compile{x}, name
	end)
	
	-- parameters that are symbolic functions -- based on f's alpha
	-- NOTICE: assign these *after* mapping stateExprs to calc
	local f = makesym'f'
	local f_param = assert(args.f_param)
	self.calc.f = f:compile{f_param}
	
	local dalpha_f = f:diff(args.f_param):simplify()
	self.calc.dalpha_f = dalpha_f:compile{args.f_param}
end

--phi = -1/(4*n) ln g_xx
-- exp(-4n phi) = g_xx for n=1
-- volume = sqrt(g_xx) = sqrt(exp(-4n phi)) = exp(-2n phi)
BSSNOK1DOriginal.graphInfos = table{
	{viewport={0/3, 0/3, 1/3, 1/3}, getter=function(self,i) return self.qs[i][1] end, name='alpha', color={1,0,1}},
	{viewport={1/3, 0/3, 1/3, 1/3}, getter=function(self,i) return self.qs[i][2] end, name='phi', color={.5,.5,1}},
	{viewport={2/3, 0/3, 1/3, 1/3}, getter=function(self,i) return self.qs[i][5] end, name='K', color={0,1,1}},
	{viewport={2/3, 1/3, 1/3, 1/3}, getter=function(self,i) return self.qs[i][6] end, name='ATilde_xx', color={0,1,1}},
	{viewport={0/3, 1/3, 1/3, 1/3}, getter=function(self,i) return self.qs[i][3] end, name='gammaTilde_xx', color={0,1,0}},
	{viewport={1/3, 1/3, 1/3, 1/3}, getter=function(self,i) return self.qs[i][4] end, name='GammaTilde^x', color={1,1,0}},
	{viewport={2/3, 2/3, 1/3, 1/3}, getter=function(self,i) return self.qs[i][1] * math.exp(-2 * self.qs[i][2]) end, name='volume', color={0,1,1}},
}
BSSNOK1DOriginal.graphInfoForNames = BSSNOK1DOriginal.graphInfos:map(function(info,i)
	return info, info.name
end)

function BSSNOK1DOriginal:initCell(sim,i)
	local x = sim.xs[i]
	local alpha = self.calc.alpha(x)
	local phi = self.calc.phi(x)
	local K = self.calc.K(x)
	local ATilde_xx = self.calc.ATilde_xx(x)
	local gammaTilde_xx = self.calc.gammaTilde_xx(x)
	--local ConnTilde_x = self.calc.ConnTilde_x(x)
	return {
		alpha = alpha,
		phi = phi,
		K = K,
		ATilde_xx = ATilde_xx,
		gammaTilde_xx = gammaTilde_xx,
	},
end

return BSSNOK1DOriginal
