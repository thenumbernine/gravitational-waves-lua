local class = require 'ext.class'
local Solver = require 'solver'
local symmath = require 'symmath'

-- source: https://en.wikipedia.org/wiki/Sod_shock_tube
-- http://www.sciencedirect.com/science/article/pii/0021999178900232
local SodExact = class(Solver)

SodExact.fixed_dt = .0001 
SodExact.name = 'Sod Exact'

--[[
args:
	x0 = initial guess
	epsilon = tolerance to stop
	maxiter = maximum iterations
	f = function
	if f is a symmath object:
		xvar = f dependent variable (for calculating the derivative)
	if f is a function:
		df = derivative function
--]]
local function newton(args)
	local f = assert(args.f)
	local df
	if symmath.Expression.is(f) then
		local xvar = assert(args.xvar)
		df = f:diff(xvar)():compile{xvar}
		f = f:compile{xvar}
	else
		df = assert(args.df)
	end
	local x = assert(args.x0)
	local epsilon = args.epsilon or 1e-50
	local maxiter = args.maxiter or 1000
	for i=1,maxiter do
		local dx = -f(x) / df(x)
		if math.abs(dx) < epsilon then break end
		x = x + dx
		if i == maxiter then
			print('warning reached maxiter='..maxiter)
		end
	end
	return x
end

function SodExact:init(args)
	SodExact.super.init(self, args)

	local gamma = self.equation.gamma
	-- TODO initial condition object to share these values with initialization
	local rhoL = 1
	local rhoR = .125
	local vL = 0
	local vR = 0
	local PL = 1
	local PR = .1

	local cL = math.sqrt(gamma * PL / rhoL)
	local cR = math.sqrt(gamma * PR / rhoR)

	-- find the zero of the sod function
	local Gamma = (gamma - 1) / (gamma + 1)
	local beta = (gamma - 1) / (2 * gamma)

	local PPostFunc, PPostDerivFunc
	do
		local P = symmath.var'P'
		local f = (P - PR) * ( (1-Gamma)*(1-Gamma) / (rhoR*(P + Gamma * PR)) )^.5
				- 2 * (gamma^.5/(gamma - 1)) * (1 - P^beta)
		PPostFunc = f:compile{P}
		local df_dp = f:diff(P)()
		PPostDerivFunc = df_dp:compile{P}
	end

	self.genConsFunc = function(t)
		
		-- 2 sqrt(gamma) / (gamma - 1) (1 - P^beta) = (P - PR) sqrt( (1 - Gamma)^2 / (rR (P + Gamma PR)) )
		local PPost = newton{
			f = PPostFunc,
			df = PPostDerivFunc,
			x0 = PR,	-- sod_exact uses pi as the initial guess.  why would't you use PR or PL or something?
		}
		
		local vPost = 2 * (math.sqrt(gamma) / (gamma - 1)) * (1 - PPost^beta)
		local rhoPost = rhoR * (PPost + Gamma * PR) / (PR + Gamma * PPost)
		local vShock = vPost * (rhoPost / rhoR) / (rhoPost / rhoR - 1)
		local rhoMiddle = rhoL * (PPost / PL) ^ (1 / gamma)
		
		local x0 = .5 * (self.domain.xmin + self.domain.xmax)
		local x1 = x0 - cL * t
		local x3 = x0 + vPost * t
		local x4 = x0 + vShock * t
		
		local c2 = cL - ((gamma - 1) / 2) * vPost
		local x2 = x0 + (vPost - c2) * t
		
		local consFunc = function(x)
			if x < x1 then	-- left
				return rhoL, vL, PL
			elseif x <= x2 then	-- rarefaction
				local c = Gamma * ((x0 - x) / t) + (1 - Gamma) * cL
				rho = rhoL * (c / cL) ^ (2 / (gamma - 1))
				v = (1 - Gamma) * ( -(x0 - x) / t + cL)
				P = PL * (rho / rhoL) ^ gamma
				return rho, v, P
			elseif x <= x3 then	-- middle
				return rhoMiddle, vPost, PPost
			elseif x <= x4 then	-- post 
				return rhoPost, vPost, PPost
			else	-- right
				return rhoR, vR, PR
			end
			error'here'
		end
	
		return consFunc
	end
end

function SodExact:step(dt)
	local consFunc = self.genConsFunc(self.t)
	for i=1,self.gridsize do
		local q = self.qs[i]
		local x = self.xs[i]
		fill(q, self.equation:calcConsFromPrim(consFunc(x)))
	end
end

return SodExact
