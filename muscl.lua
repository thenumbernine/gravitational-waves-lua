local class = require 'ext.class'
local limiter = require 'limiter' 

local function MUSCLBehavior(parentClass)

	local MUSCLTemplate = class(parentClass)

	function MUSCLTemplate:init(args)
		MUSCLTemplate.super.init(self, args)
		self.name = self.name .. ' MUSCL'

		-- limiter of ratio
		-- popular limiters: Fromm, Beam-Warming, Lax-Wendroff, minmod
		-- notice that winded-ness needs to vary per-scheme
		self.slopeLimiter = args.slopeLimiter or limiter.Fromm

		-- slopes
		self.sigma = self:newState()
		-- subgrid values
		self.iqLs = self:newState()
		self.iqRs = self:newState()
		-- flux based on subgrid values
		self.ifLs = self:newState()
		self.ifRs = self:newState()
		-- half step in time
		self.iqhLs = self:newState()
		self.iqhRs = self:newState()
	end

	-- first calculate new state interface left and right
	-- then call the old calcDT which also calculates eigen basis stuff
	function MUSCLTemplate:calcFluxes(dt)

		-- calculate the slopes
		for i=3,self.gridsize-1 do
			local qL2 = MUSCLTemplate.super.get_qL(self, i-1)
			local qL = MUSCLTemplate.super.get_qL(self, i)
			local qR = MUSCLTemplate.super.get_qR(self, i) 
			local qR2 = MUSCLTemplate.super.get_qR(self, i+1)
			for j=1,self.numStates do
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
				
				local phi = self.slopeLimiter.func(r)
				
				-- slope limiter:
				-- 	self.sigma = minmod(dq[i],dq[i']) for i current cell and i' next advection cell
				-- flux limiter:
				--  phi = minmod(1,r) = minmod(dq[i]/dq[i], dq[i']/dq[i])
				-- ... so scaling phi by dq[i] ... or simply not dividing r by dq[i] ... should get us the slope
				-- ... but if the slope limiter is constrained to [0,2] then we want it divided at first, then multiply after 
				self.sigma[i][j] = phi * dq
			end
		end

		-- subgrid values
		for j=1,self.numStates do
			self.iqLs[1][j] = self.qs[1][j]
			self.iqRs[1][j] = self.qs[1][j]
		end
		for i=2,self.gridsize do
			local qL = MUSCLTemplate.super.get_qL(self, i)
			local qR = MUSCLTemplate.super.get_qR(self, i)
			local dx = self.ixs[i] - self.ixs[i-1]
			for j=1,self.numStates do
				self.iqLs[i][j] = qL[j] + .5 * dx * self.sigma[i-1][j]
				self.iqRs[i][j] = qR[j] - .5 * dx * self.sigma[i][j]
			end
		end

		-- flux based on subgrid values
		for i=1,self.gridsize do
			fill(self.ifLs[i], self.equation:calcFluxForState(self.iqLs[i]))
			fill(self.ifRs[i], self.equation:calcFluxForState(self.iqRs[i]))
		end

		-- half step in time
		for j=1,self.numStates do
			self.iqhLs[1][j] = self.iqLs[1][j]
			self.iqhRs[1][j] = self.iqRs[1][j]
			self.iqhLs[self.gridsize][j] = self.iqLs[self.gridsize][j]
			self.iqhRs[self.gridsize][j] = self.iqRs[self.gridsize][j]
		end
		for i=2,self.gridsize-1 do
			local dx = self.ixs[i] - self.ixs[i-1]
			for j=1,self.numStates do
				self.iqhLs[i][j] = self.iqLs[i][j] + .5 * dt/dx * (self.ifLs[i][j] - self.ifRs[i-1][j])
				self.iqhRs[i][j] = self.iqRs[i][j] + .5 * dt/dx * (self.ifLs[i+1][j] - self.ifRs[i][j])
			end
		end
		
		return MUSCLTemplate.super.calcFluxes(self, dt)
	end

	function MUSCLTemplate:get_qL(i)
		return self.iqhLs[i]
	end

	function MUSCLTemplate:get_qR(i)
		return self.iqhRs[i]
	end

	return MUSCLTemplate
end

return MUSCLBehavior
