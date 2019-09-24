--[[
1982 Teukolsky "Linearized quadrupole waves in general relativity and the motion of test particles"
l=2 m=0 even wave
--]]
local class = require 'ext.class'
local table = require 'ext.table'
local Solver = require 'solver'
local symmath = require 'symmath'

local TeukolskyExact = class(Solver)

TeukolskyExact.fixed_dt = .01
TeukolskyExact.name = 'Teukolsky Exact'

function TeukolskyExact:init(...)
	TeukolskyExact.super.init(self, ...)

	local oldenv = getfenv()
	local symenv = {}
	symmath.setup{env=symenv}
	local env = setmetatable({}, {__index=function(t,k)
		local v = symenv[k] if v ~= nil then return v end 
		local v = oldenv[k] if v ~= nil then return v end
	end})
	setfenv(1, env)
	
	local x = var'x'
	local lambda = var'lambda'
	local A = var'A'
	
	local F_def = A*x*exp(-(x/lambda)^2)
	print(F_def)
	local dF_dx_defs = table{F_def}
	local maxDerivs = 5
	for i=2,maxDerivs do
		dF_dx_defs[i] = dF_dx_defs[i-1]:diff(x)() 
		print(dF_dx_defs[i])
	end

	local AVal = 1
	local lambdaVal = 1
	local dF_dx_funcs = dF_dx_defs:mapi(function(dF_dx_def)
		return dF_dx_def:compile{lambda, A, x}:bind(lambdaVal, AVal)
	end)

	--local h = math.pi / 2	-- theta
	--local p = 0				-- phi
		-- 2013 Baumgarte et al
	local h = 1.61
	local p = 4.71

	-- eqn 9.54
	-- what is 'upper sign' vs 'lowre sign'? same as plus sign vs minus sign?
	-- enumerating -2..2
	local m_plus_3 = 3
	local f_rr = ({
		math.sin(h)^2 * math.sin(2 * p),
		math.sin(2 * h) * math.sin(p),
		2 - 3 * math.sin(h)^2,
		math.sin(2 * h) * math.cos(p),
		math.sin(h)^2 * math.cos(2 * p),
	})[m_plus_3]
	local f_rh = ({
		.5 * math.sin(2 * h) * math.sin(2 * p),
		math.cos(2 * h) * math.sin(p),
		-1.5 * math.sin(2 * h),
		math.cos(2 * h) * math.cos(p),
		.5 * math.sin(2 * h) * math.cos(2 * p),
	})[m_plus_3]
	local f_rp = ({
		math.sin(h) * math.cos(2 * p),
		math.cos(h) * math.cos(p),
		0,
		math.cos(h) * -math.sin(p),
		math.sin(h) * -math.sin(2 * p),
	})[m_plus_3]
	local f1_hh = ({
		(1 + math.cos(h)^2) * math.sin(2 * p),
		math.sin(2 * h) * -math.sin(p),
		3 * math.sin(h)^2,
		math.sin(2 * h) * -math.cos(p),
		(1 + math.cos(h)^2) * math.cos(2 * p),
	})[m_plus_3]
	local f2_hh = ({
		-math.sin(2 * p),
		0,
		-1,
		0,
		-math.cos(2 * p),
	})[m_plus_3]
	local f_hp = ({
		math.cos(h) * -math.cos(2 * p),
		math.sin(h) * math.cos(p),
		0,
		math.sin(h) * -math.sin(p),
		math.cos(h) * math.sin(2 * p),
	})[m_plus_3]
	local f1_pp = ({
		(1 + math.cos(h)^2) * -math.sin(2 * p),
		math.sin(2 * h) * math.sin(p),
		-3 * math.sin(h)^2,
		math.sin(2 * h) * math.cos(p),
		(1 + math.cos(h)^2) * -math.cos(2 * p),
	})[m_plus_3]
	local f2_pp = ({
		math.cos(h)^2 * math.sin(2 * p),
		math.sin(2 * h) * -math.sin(p),
		3 * math.sin(h)^2 - 1,
		math.sin(2 * h) * -math.cos(p),
		math.cos(h)^2 * math.cos(2 * p),
	})[m_plus_3]
	
	function self.genConsFunc(t)
		return function(r)
			local alpha = 1
			local a_x = 1
			local D_g = 1
			local KTilde = 1

			-- eqn 9.49
			-- F_1 and F_2 are defined in exercise 9.4 
			local F = dF_dx_funcs:mapi(function(dF_dx_func,i)
				local F_1 = dF_dx_func(t - r)
				local F_2 = dF_dx_func(t + r) * (i%2==0 and -1 or 1)
				return F_1 + F_2, i-1	-- make this table 0-based
			end)

			-- 2008 Rinne Appendix A, except A -> A/8, B -> B/4, C -> C/8
			-- 2010 B&S eqn 9.51
			local A = 3/r^3 * (F[2] + 1/r * (3 * F[1]  + 1/r * 3 * F[0])) 
--local A = F[0]			
--local A = (t-r)*math.exp(-(t-r)^2)/r^5
			-- 2010 B&S eqn 9.52
			local B = -(F[3] / r^2 + 3 * F[2] / r^2 + 6 * F[1] / r^4 + 6 * F[0] / r^5)
			-- 2010 B&S eqn 9.53
			local C = .25 * (F[4] / r + 2 * F[3] / r^2 + 9 * F[2] / r^3 + 21 * F[1] / r^4 + 21 * F[0] / r^5)

			-- eqn 9.48
			--[[
			local g_tt = -1
			local g_rr = 1 + A*f_rr
			local g_rh = B * f_rh * r
			local g_rp = B * f_rp * r * math.sin(h)
			local g_hh = (1 + C * f1_hh + A * f2_hh) * r^2
			local g_hp = (A - 2*C) * f_hp * r^2 * math.sin(h)
			local g_pp = (1 + C * f1_pp + A * f2_pp) * r^2 * math.sin(h)^2
			--]]
			-- [[ in orthonormal non-coordinates
			--[=[
			local h_rr = A*f_rr
			local h_rh = B * f_rh
			local h_rp = B * f_rp
			local h_hh = C * f1_hh + A * f2_hh
			local h_hp = (A - 2*C) * f_hp
			local h_pp = C * f1_pp + A * f2_pp
			--]=]
			--]]
			local h_rr = A*f_rr

			local gamma_xx = 1 + h_rr
			return alpha, gamma_xx, a_x, D_g, KTilde
		end
	end
end

function TeukolskyExact:step(dt)
	local consFunc = self.genConsFunc(self.t)
	for i=1,self.gridsize do
		local q = self.qs[i]
		local x = self.xs[i]
		fill(q, consFunc(x))
	end
end

return TeukolskyExact 
