--[[

  d [ rho ]    d    [ rho ]     d    [ 0 ]
  - [rho u] +  - (u [rho u]) +  - (P [ 1 ]) = 0
 dt [rho e]   dx    [rho e]    dx    [ u ]

the 1st and 2nd terms are integrated via the flux integration
the 1st and 3rd terms are integrated via the pressure integration
	that is split into first the momentum and then the work diffusion 

--]]

local class = require 'ext.class'
local Euler1D = require 'euler1d'
local Solver = require 'solver'

local Euler1DBurgers = class(Solver)

Euler1DBurgers.equation = Euler1D()
Euler1DBurgers.name = 'Euler 1D Burgers'

function Euler1DBurgers:calcDT(getLeft, getRight)
	local gamma = self.equation.gamma
	
	-- determine timestep based on cell velocity 
	local dt
	if self.fixed_dt then
		dt = self.fixed_dt
	else
		local result = huge
		for i=1,self.gridsize do
			local rho = self.qs[i][1]
			local u = self.qs[i][2] / rho
			local eTotal = self.qs[i][3] / rho
			local eInt = eTotal - .5 * u * u
			local Cs = sqrt(gamma * (gamma - 1) * eInt)

			local dx = self.ixs[i+1] - self.ixs[i]
			local dum = dx / (Cs + abs(u))
			result = min(result, dum)
		end
		dt = result * self.cfl
	end

	return dt
end

function Euler1DBurgers:reset()
	Euler1DBurgers.super.reset(self)

	-- also in Roe and HLL
	self.fluxes = {}
	for i=1,self.gridsize+1 do
		self.fluxes[i] = {}
		for j=1,self.numStates do
			self.fluxes[i][j] = 0
		end
	end
end

function Euler1DBurgers:calcFlux(dt, getLeft, getRight, getLeft2, getRight2)
	local dq_dts = self:newState()
	
	local gamma = self.equation.gamma

	for i=3,self.gridsize-1 do
		local qL2 = getLeft2 and getLeft2(self,i) or self.qs[i-2]
		local qL = getLeft and getLeft(self,i) or self.qs[i-1]
		local qR = getRight and getRight(self,i) or self.qs[i]
		local qR2 = getRight2 and getRight2(self,i) or self.qs[i+1]
		local uL = qL[2] / qL[1]
		local uR = qR[2] / qR[1]
		local iu = .5 * (uL + uR)
		local dx = self.xs[i] - self.xs[i-1]
		for j=1,self.numStates do
			local dq = qR[j] - qL[j]
			local r = dq == 0 and 0 or 
				(iu >= 0 
					and ((qL[j] - qL2[j]) / dq) 
					or ((qR2[j] - qR[j]) / dq))
			local phi = self.fluxLimiter(r)
			local theta = iu >= 0 and 1 or -1
			self.fluxes[i][j] =
				.5 * iu * ((1 + theta) * qL[j] + (1 - theta) * qR[j])
				-- why does a + make things worse, and a - make things *much* smoother?  (original: http://www.mpia.de/homes/dullemon/lectures/fluiddynamics/Chapter_4.pdf says + )
				+ .5 * dq * phi * abs(iu) * (1 - abs(iu * dt/dx))
		end
	end

	for i=1,self.gridsize do
		local dx = self.ixs[i+1] - self.ixs[i]
		for j=1,self.numStates do
			dq_dts[i][j] = dq_dts[i][j] - (self.fluxes[i+1][j] - self.fluxes[i][j]) / dx
		end
	end

	return dq_dts
end
	
function Euler1DBurgers:postIterate(dt)
	self.Ps = self.Ps or {}
	
	local gamma = self.equation.gamma

	for i=1,self.gridsize do
		local rho = self.qs[i][1]
		local u = self.qs[i][2] / rho
		local eTotal = self.qs[i][3] / rho
		local eInt = eTotal - .5 * u * u
		self.Ps[i] = (gamma - 1) * rho * eInt
	end

	--[[ artificial viscosity
	for i=2,self.gridsize-1 do
		local rho = self.qs[i][1]
		local uL = self.qs[i-1][2] / self.qs[i-1][1]
		local uR = self.qs[i+1][2] / self.qs[i+1][1]
		local zeta = 2
		local du = zeta * .5 * (uR - uL)
		self.Ps[i] = self.Ps[i] + du * du * rho
	end
	--]]

	-- diffuse momentum
	for i=2,self.gridsize-1 do
		local dP = self.Ps[i+1] - self.Ps[i-1]
		local dx = self.xs[i+1] - self.xs[i-1]
		self.qs[i][2] = self.qs[i][2] - dP * dt / dx
	end

	-- diffuse work
	for i=2,self.gridsize-1 do
		local WR = self.Ps[i+1] * self.qs[i+1][2] / self.qs[i+1][1]
		local WL = self.Ps[i-1] * self.qs[i-1][2] / self.qs[i-1][1]
		local dW = WR - WL
		local dx = self.xs[i+1] - self.xs[i-1]
		self.qs[i][3] = self.qs[i][3] - dW * dt / dx
	end
end

return Euler1DBurgers

