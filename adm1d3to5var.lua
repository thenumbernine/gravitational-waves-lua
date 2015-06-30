--[[
3-var system extrapolated to 5-vars so the integration of alpha and g_xx are not separate
--]]

local class = require 'ext.class'
local table = require 'ext.table'

local ADM1D3to5Var = class()
ADM1D3to5Var.name = 'ADM1D3to5Var'

ADM1D3to5Var.numStates = 5

-- initial conditions
function ADM1D3to5Var:init(args, ...)
	
	local symmath = require 'symmath'
	local function makesym(field)
		return symmath.clone(assert(args[field], "expected to find field "..field)):simplify() 
	end
	
	-- parameters that are variables of symbolic functions
	local x = assert(args.x)
	
	-- parameters that are symbolic functions -- based on coordinates 
	local exprs = table{'alpha', 'g_xx', 'K_xx'}:map(function(name)
		return makesym(name), name
	end)
	
	-- derived functions
	exprs.dx_alpha = exprs.alpha:diff(x):simplify()
	exprs.dx_g_xx = exprs.g_xx:diff(x):simplify()
	
	-- convert from symbolic functions to Lua functions
	self.calc = exprs:map(function(expr, name)
		return expr:compile{x}, name
	end)
	
	-- parameters that are symbolic functions -- based on f's alpha
	-- NOTICE: assign these *after* mapping exprs to calc
	local f = makesym'f'
	local f_param = assert(args.f_param)
	self.calc.f = f:compile{f_param}
	
	local dalpha_f = f:diff(f_param):simplify()
	self.calc.dalpha_f = dalpha_f:compile{f_param}
end

ADM1D3to5Var.graphInfos = table{
	{viewport={0/3, 0/3, 1/3, 1/3}, getter=function(self,i) return self.qs[i][1] end, name='alpha', color={1,0,1}},
	{viewport={0/3, 1/3, 1/3, 1/3}, getter=function(self,i) return self.qs[i][3] end, name='A_x', color={0,1,0}},
	{viewport={1/3, 0/3, 1/3, 1/3}, getter=function(self,i) return self.qs[i][2] end, name='g_xx', color={.5,.5,1}},
	{viewport={1/3, 1/3, 1/3, 1/3}, getter=function(self,i) return self.qs[i][4] end, name='D_xxx', color={1,1,0}},
	{viewport={2/3, 0/3, 1/3, 1/3}, getter=function(self,i) return self.qs[i][5] / math.sqrt(self.qs[i][2]) end, name='K_xx', color={0,1,1}},
	{viewport={2/3, 1/3, 1/3, 1/3}, getter=function(self,i) return self.qs[i][1] * math.sqrt(self.qs[i][2]) end, name='volume', color={0,1,1}},
	{viewport={0/3, 2/3, 1/3, 1/3}, getter=function(self,i) return math.log(self.eigenbasisErrors[i]) end, name='log eigenbasis error', color={1,0,0}, range={-30, 30}},
	{viewport={1/3, 2/3, 1/3, 1/3}, getter=function(self,i) return math.log(self.fluxMatrixErrors[i]) end, name='log reconstuction error', color={1,0,0}, range={-30, 30}},
}
ADM1D3to5Var.graphInfoForNames = ADM1D3to5Var.graphInfos:map(function(info,i)
	return info, info.name
end)

function ADM1D3to5Var:initCell(sim,i)
	local x = sim.xs[i]
	local alpha = self.calc.alpha(x)
	local g_xx = self.calc.g_xx(x)
	local A_x = self.calc.dx_alpha(x) / self.calc.alpha(x)
	local D_xxx = 1/2 * self.calc.dx_g_xx(x)
	local K_xx = self.calc.K_xx(x)
	local KTilde_xx = K_xx / sqrt(g_xx)
	return {alpha, g_xx, A_x, D_xxx, KTilde_xx}
end

function ADM1D3to5Var:calcInterfaceEigenBasis(sim,i,qL,qR)
	local avgQ = {}
	for j=1,self.numStates do 
		avgQ[j] = (qL[j] + qR[j]) / 2
	end
	local alpha, g_xx, A_x, D_xxx, KTilde_xx = unpack(avgQ)
	local f = self.calc.f(alpha)
	local lambda = alpha * sqrt(f / g_xx)		
	sim.eigenvalues[i] = {-lambda, 0, 0, 0, lambda}
end

return ADM1D3to5Var

