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
	if symmath.Expression:isa(f) then
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

	-- find the zero of the sod function
	local muSq = (gamma - 1) / (gamma + 1)
	local K = PL / rhoL^gamma
	--local beta = (gamma - 1) / (2 * gamma)
	
	local CsL = math.sqrt(gamma * PL / rhoL)
	local CsR = math.sqrt(gamma * PR / rhoR)

	local function solveP3()
		local P3 = .5 * (PL + PR)
		local epsilon = 1e-16
		while true do
			local f = ((((-2 * CsL) * (1 - ((P3 / PL) ^ ((-1 + gamma) / (2 * gamma))))) / (CsR * (-1 + gamma))) + ((-1 + (P3 / PR)) * ((0.75 / (gamma * (0.25 + (P3 / PR)))) ^ 0.5))) 
			local df_dP3 = ((-((((((1.5 * math.sqrt(0.75) * CsR * PR * (gamma ^ 1.5)) - ((0.75 ^ 1.5) * CsR * PR * math.sqrt(gamma))) - ((0.75 ^ 1.5) * CsR * (gamma ^ 2.5) * PR)) - (0.5 * P3 * math.sqrt(0.75) * CsR * math.sqrt(gamma))) - (0.5 * P3 * math.sqrt(0.75) * CsR * (gamma ^ 2.5))) + (((P3 * math.sqrt(0.75) * CsR * (gamma ^ 1.5)) - (0.25 * (PL ^ ((1 - gamma) / (2 * gamma))) * CsL * (P3 ^ (((-1) - gamma) / (2 * gamma))) * (PR ^ 1.5) * math.sqrt((P3 + (0.25 * PR))))) - ((PL ^ ((1 - gamma) / (2 * gamma))) * CsL * (P3 ^ ((-(1 - gamma)) / (2 * gamma))) * math.sqrt(PR) * math.sqrt((P3 + (0.25 * PR))))) + (0.5 * (PL ^ ((1 - gamma) / (2 * gamma))) * CsL * (P3 ^ (((-1) - gamma) / (2 * gamma))) * (PR ^ 1.5) * gamma * math.sqrt((P3 + (0.25 * PR)))) + (((2 * (PL ^ ((1 - gamma) / (2 * gamma))) * CsL * (P3 ^ ((-(1 - gamma)) / (2 * gamma))) * math.sqrt(PR) * gamma * math.sqrt((P3 + (0.25 * PR)))) - (0.25 * (PL ^ ((1 - gamma) / (2 * gamma))) * CsL * (P3 ^ (((-1) - gamma) / (2 * gamma))) * (gamma ^ 2) * (PR ^ 1.5) * math.sqrt((P3 + (0.25 * PR))))) - ((PL ^ ((1 - gamma) / (2 * gamma))) * CsL * (P3 ^ ((-(1 - gamma)) / (2 * gamma))) * (gamma ^ 2) * math.sqrt(PR) * math.sqrt((P3 + (0.25 * PR))))))) / (math.sqrt(PR) * CsR * ((P3 + (0.25 * PR)) ^ 1.5) * gamma * ((1 - (2 * gamma)) + (gamma ^ 2)))) 
			local dP3 = -f / df_dP3
			if math.abs(dP3) <= epsilon then break end
			if not math.isfinite(dP3) then error('delta is not finite! '..tostring(dP3)) end
			P3 = P3 + dP3 
		end
		return P3
	end
			
	local P3 = solveP3()
	local P4 = P3
	
	local rho3 = rhoL * (P3 / PL) ^ (1 / gamma)
	
	local v3 = vR + 2 * CsL / (gamma - 1) * (1 - (P3 / PL)^((gamma - 1)/(2*gamma)))
	local v4 = v3
	
	local rho4 = rhoR * (P4 + muSq * PR) / (PR + muSq * P4)
	
	local vshock = v4 * rho4 / (rho4 - rhoR)
	local vtail = CsL - v4 / (1 - muSq)
	
	-- between regions 1 and 2
	local s1 = -CsL	
	
	-- between regions 2 and 3
	-- http://www.itam.nsc.ru/flowlib/SRC/sod.f
	local s2 = -vtail

	local s3 = v3	-- between regions 3 and 4

	-- between regions 4 and 5 ...
	local s4 = vshock

	self.genConsFunc = function(t)
		local consFunc = function(x)
			--print('wavespeeds:',s1,s2,s3,s4)

			local rho, vx, P
			local xi = x / t
			if xi < s1 then
				rho = rhoL
				vx = vL
				P = PL
			elseif xi < s2 then
				vx = (1 - muSq) * (x/t + CsL)
				
				-- Dullemon:
				--rho = (rhoL^gamma / (gamma * PL) * (v2(x) - x/t)^2)^(1/(gamma-1))
				-- http://www.itam.nsc.ru/flowlib/SRC/sod.f
				rho = rhoL * (-muSq * (x / (CsL * t)) + (1 - muSq))^(2/(gamma-1))

				-- Dullemon:
				--P = K * rho2^gamma
				-- http://www.itam.nsc.ru/flowlib/SRC/sod.f
				P = PL * (-muSq * (x / (CsL * t)) + (1 - muSq)) ^ (2*gamma/(gamma-1))
			elseif xi < s3 then
				rho = rho3
				vx = v3
				P = P3
			elseif xi < s4 then
				rho = rho4
				vx = v4
				P = P4
			else
				rho = rhoR
				vx = vR
				P = PR
			end

			local EInt = P / (gamma - 1)
			local EKin = .5 * rho * vx*vx
			local ETotal = EKin + EInt
			return rho, rho * vx, ETotal	
		end
	
		return consFunc
	end
end

function SodExact:step(dt)
	local consFunc = self.genConsFunc(self.t)
	for i=1,self.gridsize do
		local q = self.qs[i]
		local x = self.xs[i]
		fill(q, consFunc(x))
	end
end

return SodExact
