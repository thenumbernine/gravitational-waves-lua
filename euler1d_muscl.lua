local class = require 'ext.class'
local fluxLimiters = require 'limiter' 
local Solver = require 'solver'

local Euler1DMUSCL = class(Solver)

Euler1DMUSCL.equation = require 'euler1d'()

function Euler1DMUSCL:init(args)
	-- baseScheme can be Roe or HLL
	-- TODO should MUSCL-Roe be using eigenbasis-transformed states somewhere in there? 
	self.baseScheme = args.baseScheme or require 'roe'()

	-- limiter of ratio
	-- popular limiters: Fromm, Beam-Warming, Lax-Wendroff, minmod
	-- notice that winded-ness needs to vary per-scheme
	self.slopeLimiter = args.slopeLimiter or fluxLimiters.Fromm
end

-- first calculate new state interface left and right
-- then call the old calcDT which also calculates eigen basis stuff
function Euler1DMUSCL:calcDT(sim, getLeft, getRight, getLeft2, getRight2)
	
	local sigma = sim:newState()
	for i=3,sim.gridsize-1 do
		local qL2 = getLeft2 and getLeft2(i) or sim.qs[i-2]
		local qL = getLeft and getLeft(i) or sim.qs[i-1]
		local qR = getRight and getRight(i) or sim.qs[i]
		local qR2 = getRight2 and getRight2(i) or sim.qs[i+1]
		for j=1,sim.numStates do
			local dq = qR[j] - qL[j]

			-- since we're doing flux/subgrid per-state (un-transformed)
			--  then, for determining the next-cell of advection,
			--   I'm going to use the velocity
			--   ...which means this is a Euler-only MUSCL solver
			local iu = .5*(qL[2]/qL[1] + qR[2]/qR[1])
			
			-- ratio of slope to next cell slope
			local r = dq == 0 and 0 or 
				(iu >= 0 
					and ((qL[j] - qL2[j]) / dq) 
					or ((qR2[j] - qR[j]) / dq))
			
			local phi = self.slopeLimiter(r)
			
			-- slope limiter:
			-- 	sigma = minmod(dq[i],dq[i']) for i current cell and i' next advection cell
			-- flux limiter:
			--  phi = minmod(1,r) = minmod(dq[i]/dq[i], dq[i']/dq[i])
			-- ... so scaling phi by dq[i] ... or simply not dividing r by dq[i] ... should get us the slope
			-- ... but if the slope limiter is constrained to [0,2] then we want it divided at first, then multiply after 
			sigma[i][j] = phi * dq
		end
	end

	-- subgrid values
	local iqLs = sim:newState()
	local iqRs = sim:newState()
	for j=1,sim.numStates do
		iqLs[1][j] = sim.qs[1][j]
		iqRs[1][j] = sim.qs[1][j]
	end
	for i=2,sim.gridsize do
		local qL = getLeft and getLeft(i) or sim.qs[i-1]
		local qR = getRight and getRight(i) or sim.qs[i]
		local dx = sim.ixs[i] - sim.ixs[i-1]
		for j=1,sim.numStates do
			iqLs[i][j] = qL[j] + .5 * dx * sigma[i-1][j]
			iqRs[i][j] = qR[j] - .5 * dx * sigma[i][j]
		end
	end

	-- flux based on subgrid values
	local ifLs = sim:newState()
	local ifRs = sim:newState()
	for i=1,sim.gridsize do
		ifLs[i] = sim:calcFluxForState(iqLs[i])
		ifRs[i] = sim:calcFluxForState(iqRs[i])
	end

	-- how should we get this dt?
	-- based on Roe eigenvector without MUSCL?  or based on Burgers?  or fixed + implicit (optimistically)
	-- should this be the dt that we consistently use even after using MUSCL to adjust the state?
	local dt = self.baseScheme.calcDT(self, sim)

	-- half step in time
	self.iqhLs = sim:newState()
	self.iqhRs = sim:newState()
	for j=1,sim.numStates do
		self.iqhLs[1][j] = iqLs[1][j]
		self.iqhRs[1][j] = iqRs[1][j]
		self.iqhLs[sim.gridsize][j] = iqLs[sim.gridsize][j]
		self.iqhRs[sim.gridsize][j] = iqRs[sim.gridsize][j]
	end
	for i=2,sim.gridsize-1 do
		local dx = sim.ixs[i] - sim.ixs[i-1]
		for j=1,sim.numStates do
			self.iqhLs[i][j] = iqLs[i][j] + .5 * dt/dx * (ifLs[i][j] - ifRs[i-1][j])
			self.iqhRs[i][j] = iqRs[i][j] + .5 * dt/dx * (ifLs[i+1][j] - ifRs[i][j])
		end
	end
	
	-- once we have *this* collection of subgrid L & R states,
	--  we use them for whatever method you want ...
	-- ... be it HLL or Roe, etc 

	self.getLeft = function(i) return self.iqhLs[i] end
	self.getRight = function(i) return self.iqhRs[i] end
	self.getLeft2 = function(i) return self.iqhLs[i-1] end
	self.getRight2 = function(i) return self.iqhRs[i+1] end

	--[[
	TODO now the paper has officially gave me a circular dependency:
	to do the MUSCL subgrid half-step in time it needs to know dt
	however, to calculate dt, it needs the eigenvalues
		eigenvalues need Roe matrix at interface
		interface needs the left and right half-step values
		... which need the dt
	--]]
	local dt = self.baseScheme.calcDT(self, sim, self.getLeft, self.getRight)

	return dt
end

function Euler1DMUSCL:calcFlux(sim, dt)
	return self.baseScheme:calcFlux(sim, dt, self.getLeft, self.getRight, self.getLeft2, self.getRight2)
end

function Euler1DMUSCL:postIterate(...)
	if self.baseScheme.postIterate then
		return self.baseScheme:postIterate(...)
	end
end

function Euler1DMUSCL:addSourceToDeriv(...)
	if self.baseScheme.addSourceToDeriv then
		self.baseScheme:addSourceToDeriv(...)
	end
end

return Euler1DMUSCL

