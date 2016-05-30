--[[
from Toro, "Riemann Solvers and Numerical Methods for Fluid Dynamics" p.163
--]]

local class = require 'ext.class'
local Solver = require 'solver'

local Euler1DGodunov = class(Solver)
Euler1DGodunov.name = 'Euler 1D Godunov'
Euler1DGodunov.equation = require 'euler1d'()

function Euler1DGodunov:init(args, ...)
	Euler1DGodunov.super.init(self, args, ...)

	local godunovMethod = (args.godunovMethod or 'exact'):lower()
	local sampleMethod = (args.sampleMethod or 'toro'):lower()
	self.name = self.name .. ' ' .. godunovMethod .. ' ' .. sampleMethod

	self.calcPressureAndVelocity = ({
		exact = self.calcPressureAndVelocityExact,
		pvrs = self.calcPressureAndVelocityPVRS,
		twoshock = self.calcPressureAndVelocityTwoShock,
		adaptive = self.calcPressureAndVelocityAdaptive,
	})[godunovMethod]

	self.sample = ({
		toro = self.sampleToro,
		alt = self.sampleAlt,
	})[sampleMethod]
end

-- same as Euler1DBuregers:calcDT
function Euler1DGodunov:calcDT()
	local gamma = self.equation.gamma
	
	-- determine timestep based on cell velocity 
	local dt
	if self.fixed_dt then
		dt = self.fixed_dt
	else
		local result = huge
		for i=2,self.gridsize do
			local q = self.qs[i]
			local rho = q[1]
			local u = q[2] / rho
			local c = math.sqrt(gamma * (gamma - 1) * (q[3] / rho - .5 * u * u))

			local dx = self.ixs[i+1] - self.ixs[i]
			local dum = dx / (c + math.abs(u))
			result = math.min(result, dum)
		end
		dt = result * self.cfl
	end
	
	return dt
end
	
-- same as Euler1DBurgers:reset()
function Euler1DGodunov:reset()
	Euler1DGodunov.super.reset(self)

	-- also in Roe and HLL
	self.fluxes = {}
	for i=1,self.gridsize+1 do
		self.fluxes[i] = {}
		for j=1,self.numStates do
			self.fluxes[i][j] = 0
		end
	end
end

function Euler1DGodunov:calcFlux(dt)
	local gamma = self.equation.gamma
	local x_t = 0

	for i=2,self.gridsize do
		local qL = self.qs[i-1]
		local qR = self.qs[i]
		
		local rhoL = qL[1]
		local uL = qL[2] / rhoL
		local PL = (gamma - 1) * (qL[3] - .5 * rhoL * uL * uL)
		local cL = math.sqrt(gamma * PL / rhoL) 
		
		local rhoR = qR[1]
		local uR = qR[2] / rhoR
		local PR = (gamma - 1) * (qR[3] - .5 * rhoR * uR * uR)
		local cR = math.sqrt(gamma * PR / rhoR)

		-- both of these are modular ...
		-- and give the same results for either ...
		-- are they all wrong?
		local PM, uM = self:calcPressureAndVelocity(rhoL, uL, PL, cL, rhoR, uR, PR, cR)
		local rho, u, P = self:sample(rhoL, uL, PL, cL, rhoR, uR, PR, cR, PM, uM, x_t)

		self.fluxes[i][1] = rho * u
		self.fluxes[i][2] = rho * u * u + P
		self.fluxes[i][3] = (gamma / (gamma - 1) * P + .5 * rho * u * u) * u
	end
	
	local dq_dts = self:newState()
	for i=1,self.gridsize do
		local dx = self.ixs[i+1] - self.ixs[i]
		for j=1,self.numStates do
			dq_dts[i][j] = (self.fluxes[i][j] - self.fluxes[i+1][j]) / dx
		end
	end
	return dq_dts
end

function Euler1DGodunov:calcPressureAndVelocityPVRS(rhoL, uL, PL, cL, rhoR, uR, PR, cR)
--print'pvrs'
	-- purpose: to compute PM and uM  in the star region using
	--          the pvrs riemann solver.
	--          we use exact relations for density and exact
	--          solution for sonic flow in sampling routine
	-- declaration of variables

	-- compute PM and uM from pvrs riemann solver
	local cup = .25*(rhoL + rhoR)*(cL + cR)
	local PM  = .5*(PL + PR) + .5*(uL - uR)*cup
	-- reset pressure if negative
	PM  = math.max(0, PM)
	local uM  = .5*(uL + uR) + .5*(PL - PR)/cup
	return PM, uM
end

function Euler1DGodunov:calcPressureAndVelocityTwoShock(rhoL, uL, PL, cL, rhoR, uR, PR, cR)
--print'twoshock'
	local gamma = self.equation.gamma
--     purpose: to compute pm and um  in the star region using
--              the two-shock riemann solver.
--              we use exact relations for density and exact
--              solution for sonic flow in sampling routine
--     declaration of variables
--     compute guess pressure from pvrs riemann solver
	local cup = .25*(rhoL + rhoR)*(cL + cR)
	local ppv = .5*(PL + PR) + .5*(uL - uR)*cup
	local ppv = math.max(0, ppv)
--     two-shock riemann solver with pvrs as estimate
	local gel = math.sqrt((2/(gamma+1)/rhoL)/((gamma-1)/(gamma+1)*PL + ppv))
	local ger = math.sqrt((2/(gamma+1)/rhoR)/((gamma-1)/(gamma+1)*PR + ppv))
	local PM  = (gel*PL + ger*PR - (uR - uL))/(gel + ger)
	local uM  = .5*(uL + uR) + .5*(ger*(PM - PR) - gel*(PM - PL))
	return PM, uM
end

function Euler1DGodunov:calcPressureAndVelocityAdaptive(dl, ul, pl, cl, dr, ur, pr, cr)
--print'adaptive'
	local gamma = self.equation.gamma
	-- purpose: to compute PM and uM  in the star region using
	--    the adaptive riemann solver: pvrs, trrs, tsrs.
	--    we use exact relations for density and exact
	--    solution for sonic flow in sampling routine
	-- declaration of variables
	local quser = 2
	-- compute guess pressure from pvrs riemann solver
	local cup  = .25*(dl + dr)*(cl + cr)
	local ppv  = .5*(pl + pr) + .5*(ul - ur)*cup
	local ppv  = max(0, ppv)
	local pmin = min(pl,  pr)
	local pmax = max(pl,  pr)
	local qmax = pmax/pmin
	local PM, uM
	if qmax <= quser and pmin <= ppv and ppv <= pmax then
	-- select pvrs riemann solver
		PM = ppv
		uM = .5*(ul + ur) + .5*(pl - pr)/cup
	else
		if ppv < pmin then
	--		  select two-rarefaction riemann solver
			local pq  = (pl/pr)^((gamma-1)/(2*gamma))
			uM  = (pq*ul/cl + ur/cr + (2/(gamma-1))*(pq - 1))/(pq/cl + 1/cr)
			local ptl = 1 + (gamma-1)/2*(ul - uM)/cl
			local ptr = 1 + (gamma-1)/2*(uM - ur)/cr
			PM  = 0.5*(pl*ptl^((2*gamma)/(gamma-1)) + pr*ptr^((2*gamma)/(gamma-1)))
		else
	--		  select two-shock riemann solver with pvrs as estimate
			local gel = math.sqrt((2/(gamma+1)/dl)/((gamma-1)/(gamma+1)*pl + ppv))
			local ger = math.sqrt((2/(gamma+1)/dr)/((gamma-1)/(gamma+1)*pr + ppv))
			PM  = (gel*pl + ger*pr - (ur - ul))/(gel + ger)
			uM  = 0.5*(ul + ur) + 0.5*(ger*(PM - pr) - gel*(PM - pl))
		end
	end
	return PM, uM
end


function Euler1DGodunov:calcPressureAndVelocityExact(rhoL, uL, PL, cL, rhoR, uR, PR, cR)
--print'exact'
	local tolerance = 0
	local maxiters = 100
	local PStart = self:guessPressure(rhoL, uL, PL, cL, rhoR, uR, PR, cR)
	local P = PStart
	local uDiff = uR - uL
	local fR, fL, fdR, fdL
	for i=1,maxiters do
		fL, fdL = self:pressureFunction(P, rhoL, PL, cL)
		fR, fdR = self:pressureFunction(P, rhoR, PR, cR)
		local PNew = P - (fL + fR + uDiff) / (fdL + fdR)
		local change = 2 * math.abs((PNew - P) / (PNew + P))
		--print('change',change)
		if change <= tolerance then
			P = PNew
			--print('converged after '..i..' iterations')
			break
		end
		if PNew < 0 then PNew = tolerance end
		P = PNew
		if i == numIters then
			error'divergence in newton-raphson iteration'
		end
	end
	-- compute velocity in star region
	local u = .5 * (uL + uR + fR - fL)
	return P, u
end

function Euler1DGodunov:guessPressure(rhoL, uL, PL, cL, rhoR, uR, PR, cR)
	local gamma = self.equation.gamma
--
--     purpose: to provide a guessed value for pressure
--	        in the star region
--     declaration of variables
	local quser = 2
--     compute guess pressure from pvrs riemann solver
	local cup  = .25*(rhoL + rhoR)*(cL + cR)
	local ppv  = .5*(PL + PR) + .5*(uL - uR)*cup
	local ppv  = math.max(0, ppv)
	local pmin = math.min(PL,  PR)
	local pmax = math.max(PL,  PR)
	local qmax = pmax/pmin
	if qmax <= quser and  pmin <= ppv and ppv <= pmax then
--	  select pvrs riemann solver
		return ppv
	else
		if ppv < pmin then
--		  select two-rarefaction riemann solver
			local pq  = (PL / PR) ^ ((gamma - 1) / (2 * gamma))
			local um  = (pq * uL / cL + uR / cR + (2 / (gamma - 1)) * (pq - 1)) / (pq / cL + 1 / cR)
			local ptl = 1 + .5 * (gamma - 1) * (uL - um) / cL
			local ptr = 1 + .5 * (gamma - 1) * (um - uR) / cR
			return .5 * (PL * ptl ^ ((2 * gamma) / (gamma - 1)) + PR * ptr ^ ((2 * gamma) / (gamma - 1)))
		else
--		  select two-shock riemann solver with pvrs as estimate
			local gel = math.sqrt((2 / (gamma + 1) / rhoL) / ((gamma - 1) / (gamma + 1) * PL + ppv))
			local ger = math.sqrt((2 / (gamma + 1) / rhoR) / ((gamma - 1) / (gamma + 1) * PR + ppv))
			return (gel * PL + ger * PR - (uR - uL)) / (gel + ger)
		end
	end
end


function Euler1DGodunov:pressureFunction(P, rhoK, PK, cK)
	local gamma = self.equation.gamma
	if P <= PK then	-- rarefaction wave
		local PRatio = P / PK
        local f = (2 / (gamma - 1)) * cK * (PRatio ^ (.5 * (gamma - 1) / gamma) - 1)	-- Toro 5.85 says should be scaled by (1 - b * rhoK)
        local fd = (1 / (rhoK * cK)) * PRatio ^ (-(gamma + 1) / (2 * gamma))
		return f, fd
	else	-- shock wave
		local aK = 2 / (gamma + 1) / rhoK			-- Toro 4.87 says should be scaled by (1 - b * rhoK)
		local bK = (gamma - 1) / (gamma + 1) * PK	-- Toro 4.87
		local c = math.sqrt(aK / (bK + P))			-- Toro 4.85
		local f = (P - PK) * c						-- Toro 4.85
		local fd = (1 - .5 * (P - PK) / (bK + P)) * c
		return f, fd
	end
end

-- [[ from Toro's code on p.163
function Euler1DGodunov:sampleToro(rhoL, uL, PL, cL, rhoR, uR, PR, cR, PM, uM, slope)
	local gamma = self.equation.gamma
--     purpose: to sample the solution throughout the wave pattern
--     declaration of variables
	if slope <= uM then
--	  sampling point lies to the left of the contact discontinuity
		if PM <= PL then
--		  left rarefaction
			local shl = uL - cL
			if slope <= shl then
--			  sampled point is left data state
--print'left'
				return rhoL, uL, PL
			else
				local cml = cL * (PM / PL) ^ ((gamma - 1) / (2 * gamma))
				local stl = uM - cml
				if slope > stl then
--				  sampled point is star left state
					local rhoM = rhoL * (PM / PL) ^ (1 / gamma)
--print'leftstar'
					return rhoM, uM, PM
				else
--				  sampled point is inside left fan
					local u = 2 / (gamma + 1) * (cL + (gamma - 1) / 2 * uL + slope)
					local c = 2 / (gamma + 1) * (cL + (gamma - 1) / 2 * (uL - slope))
					local rho = rhoL * (c / cL) ^ (2 / (gamma - 1))
					local P = PL * (c / cL) ^ (2 * gamma / (gamma - 1))
--print'leftfan'
					return rho, u, P
				end
			end
		else
--		  left shock
			local pml = PM/PL
			local sl = uL - cL * math.sqrt((gamma + 1) / (2 * gamma) * pml + ((gamma - 1) / (2 * gamma)))
			if slope <= sl then
--			  sampled point is left data state
--print'left2'
				return rhoL, uL, PL
			else
--			  sampled point is star left state
				local rhoM = rhoL * (pml + (gamma - 1) / (gamma + 1)) / (pml * (gamma - 1) / (gamma + 1) + 1)
--print'leftstar2'
				return rhoM, uM, PM
			end
		end
	else
--	  sampling point lies to the right of the contact discontinuity
		if PM > PR then
--		  right shock
			local pmr = PM/PR
			local sr = uR + cR * math.sqrt((gamma + 1) / (2 * gamma) * pmr + ((gamma - 1) / (2 * gamma)))
			if slope >= sr then
--			  sampled point is right data state
--print'right'
				return rhoR, uR, PR
			else
--			  sampled point is star right state
				local rhoM = rhoR * (pmr + (gamma - 1) / (gamma + 1)) / (pmr * (gamma - 1) / (gamma + 1) + 1)
--print'rightstar'
				return rhoM, uM, PM
			end
		else
--		  right rarefaction
			local shr = uR + cR
			if slope >= shr then
--			  sampled point is right data state
--print'right2'
				return rhoR, uR, PR
			else
				local cmr = cR * (PM / PR) ^ ((gamma - 1) / (2 * gamma))
				local str = uM + cmr
				if slope <= str then
--				  sampled point is star right state
					local rhoM = rhoR * (PM / PR) ^ (1 / gamma)
--print'rightstar2'
					return rhoM, uM, PM
				else
--				  sampled point is inside right fan
					local u = 2 / (gamma + 1) * (-cR + (gamma - 1) / 2 * uR + slope)
					local c = 2 / (gamma + 1) * (cR - (gamma - 1) / 2 * (uR - slope))
					local rho = rhoR * (c / cR) ^ (2 / (gamma - 1))
					local P = PR * (c / cR) ^ (2 * gamma / (gamma - 1))
--print'rightfan'
					return rho, u, P
				end
			end
		end
	end
end
--]]

-- [[ from http://coffeecfd.blogspot.com/2014/03/some-compressible-cfd-codes-for.html -- verifying the accuracy of the Toro sample method
function Euler1DGodunov:sampleAlt(rhoL, uL, PL, cL, rhoR, uR, PR, cR, PM, uM, slope)

	--asshat who wrote this baked all this numbers ...
	local K = 1.4
	local K_1_2 = 0.2
	local K_1 = 0.4
	local K_K = 0.142857142857
	local _2K_K_1 = 7.0
	local K__1_2 = 1.2
	local K__1_2K = 0.857142857143
	local K_1_2K = 0.142857142857
	local K_1_K__1 = 0.166666666667
	local K_K_1 = 3.5
	local _2_K__1 = 0.833333333333
	local K__1 = 2.4
	local _3K_1 = 3.2
	local _2_K_1 = 5.0
	local _4K = 5.6

    -- left wave is SW
	local aL
	if PM >= PL then
		aL = math.sqrt(rhoL * (K__1_2 * PM + K_1_2 * PL))
    -- left wave is RW
	elseif PL - PM > 300 then
		aL = K_1_2K * rhoL * cL * (1 - PM / PL) / (1 - (PM/PL)^K_K)
	else
		aL = rhoL * cL
	end

	-- rigt wave is SW
	local aR
	if PM > PR then
		aR = math.sqrt(rhoR * (K__1_2 * PM + K_1_2 * PR))
	-- rigt wave is RW
	elseif PR - PM > 300 then
		aR = K_1_2K * rhoR * cR * (1 - PM / PR) / (1 - (PM/PR)^K_K)
	else
		aR = rhoR * cR
	end

	-- velocity of CD
	uM = (aL * uL + aR * uR + PL - PR) / (aL + aR)

	-- characteristic velocities
	-- left wave is SW
	if PM > PL then
		dL = uL - aL / rhoL
		dL_ = dL
		RhoL = rhoL * aL / (aL - rhoL * (uL - uM))
	-- left wave is RW
	else
		dL = uL - cL
		cL_ = cL + K_1_2 * (uL - uM)
		dL_ = uM - cL_
		RhoL = K * PM / (cL_ * cL_)
	end
	
	-- right wave is SW
	if PM > PR then
		dR = uR + aR / rhoR
		dR_ = dR
		RhoR = rhoR * aR / (aR + rhoR * (uR - uM))
	-- right wave is RW
	else
		dR = uR + cR
		cR_ = cR - K_1_2 * (uR - uM)
		dR_ = uM + cR_
		RhoR = K * PM / (cR_ * cR_)
	end

	-- parameters from zone 3 or 4
	if dL_ < 0 and dR_ > 0 then
		if uM >= 0 then
			rhoM = rhoL
		else
			rhoM = rhoR
		end
		return rhoM, uM, PM
	end

	-- parameters from zone 2
	if dL < 0 and dL_ >= 0 then
		cL_ = _2_K__1 * cL + K_1_K__1 * uL
		uM = cL_
		PM = PL * (cL_/cL)^_2K_K_1
		rhoM = K * PM / (cL_ * cL)
		return rhoM, uM, PM
	end

	-- parameters from zone 5
	if dR_ <= 0 and dR > 0 then
		cR_ = _2_K__1 * cR - K_1_K__1 * uR
		uM = -cR_
		PM = PR * (cR_/cR)^_2K_K_1
		rhoM = K * PM / (cR_ * cR_)
		return rhoM, uM, PM
	end

	if dL >= 0 then
		return rhoL, uL, PL
	end

	return rhoR, uR, PR
end
--]]

return Euler1DGodunov
