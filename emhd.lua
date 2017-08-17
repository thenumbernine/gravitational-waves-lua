--[[
primitive variables
rho, vx, vy, vz, P
Ex, Ey, Ez
Bx, By, Bz

state variables
rho, mx, my, mz, ETotal
Ex, Ey, Ez
Bx, By, Bz

time advance equations:
rho,t + (rho u^j),j = 0		(continuity equation)
(rho u^i),t + (rho u^i u^j + delta^ij P),j = n q (E + u x B)^i		(momentum, lorentz force law as source term)
ETotal,t + ((ETotal + P) u^j),j = n q (u dot E) 	(energy total)
E^i,t - 1/(mu0 eps0) (e_ijk B^k),j = -1/eps0 n q u^i
B^i,t + (e_ijk E^k),j = 0

relations:
rho = n m	(density = particle number x particle mass)
ETotal = P / (gamma - 1) + rho u^2 / 2	(total energy = internal ideal gas energy + kinetic energy)
P = n kB T 	(pressure related to temperature)

... wow, this is just like Euler fluid + Maxwell 
... but with each eqns vars in the others source terms

so here's the idea:

make a Euler1D and a Maxwell eqn, then spend the rest of the routines picking them apart and putting them together
...and make sure your units are correctly scaled, so the cfl is correct 

--]]
local class = require 'ext.class'
local Equation = require 'equation'

local Euler3D = require 'euler3d'
local Maxwell = require 'maxwell'

local EMHD = class(Equation)
EMHD.name = 'EMHD'

EMHD.species = {
	{mass=1, charge=1},
}

EMHD.numStates = #EMHD.species * Euler3D.numStates + Maxwell.numStates
EMHD.gamma = Euler3D.gamma
EMHD.mu0 = 1
EMHD.epsilon0 = 1

do
	local U = function(solver, i) return solver.qs[i] end
	local gamma = function(solver) return solver.equation.gamma end
	local mu0 = function(solver) return solver.equation.mu0 end
	local epsilon0 = function(solver) return solver.equation.epsilon0 end
	local rho = U:_(1)
	local mx, my, mz = U:_(2), U:_(3), U:_(4)
	local EnergyHydro = U:_(5)
	local vx, vy, vz = mx/rho, my/rho, mz/rho
	local EnergyKinetic = .5 * rho * (vx*vx + vy*vy + vz*vz)
	local EnergyInternal = EnergyHydro - EnergyKinetic
	local P = (gamma - 1) * EnergyInternal
	local S = P / rho^gamma
	local H = EnergyInternal + P 
	local h = H / rho 
	local HTotal = EnergyHydro + P 
	local hTotal = HTotal / rho
	local Ex = U:_(6) / epsilon0
	local Ey = U:_(7) / epsilon0
	local Ez = U:_(8) / epsilon0
	local ESq = Ex^2 + Ey^2 + Ez^2
	local E = math.sqrt:o(ESq)
	local Bx = U:_(9)
	local By = U:_(10)
	local Bz = U:_(11)
	local BSq = Bx^2 + By^2 + Bz^2
	local B = math.sqrt:o(BSq)
	local EnergyEM = .5 * (ESq * epsilon0 + BSq / mu0)

	EMHD:buildGraphInfos{
		-- Euler 3D
		-- prims
		{rho = rho},
		{vx = vx},
		{vy = vy},
		{vz = vz},
		{P = P},
		-- other vars
		{EnergyInternal = EnergyInternal},
		{EnergyKinetic = EnergyKinetic},
		{H = H},
		{S = S},
		-- conservative
		--{mx = mx},
		{EnergyHydro = EnergyHydro},

		-- Maxwell
		{Ex=Ex}, {Ey=Ey}, {Ez=Ez}, {E=E},
		{Bx=Bx}, {By=By}, {Bz=Bz}, {B=B},
		{EnergyEM = EnergyEM},
	}
end

function EMHD:init(...)
	self.maxwell = Maxwell(...)
	self.euler = Euler3D(...)
end

function EMHD:initCell(solver, i)
	local x = solver.xs[i]
	-- Brio & Wu
	self.gamma = 2
	local rho = x < 0 and 1 or .125
	local vx, vy, vz = 0, 0, 0
	local P = x < 0 and 1 or .1
	local Ex, Ey, Ez = 0, 0, 0
	local Bx, By, Bz = .75, (x < 0 and 1 or -1), 0
	return {self:calcConsFromPrim(rho, vx, vy, vz, P, Ex, Ey, Ez, Bx, By, Bz)}
end

function EMHD:calcConsFromPrim(rho, vx, vy, vz, P, Ex, Ey, Ez, Bx, By, Bz)
	local rho, mx, my, mz, EnergyTotal = self.euler:calcConsFromPrim(rho, vx, vy, vz, P)
	return rho, mx, my, mz, EnergyTotal, Ex, Ey, Ez, Bx, By, Bz
end

function EMHD:calcPrimFromCons(rho, mx, my, mz, EnergyTotal, Ex, Ey, Ez, Bx, By, Bz)
	local rho, vx, vy, vz, P = self.euler:calcPrimFromCons(rho, mx, my, mz, EnergyTotal)
	return rho, vx, vy, vz, P, Ex, Ey, Ez, Bx, By, Bz
end

function EMHD:calcRoeValues(qL, qR)
	local rho, vx, vy, vz, hTotal, Cs = self.euler:calcRoeValues(qL, qR)
	local Ex, Ey, Ez, Bx, By, Bz = self.maxwell:calcRoeValues({table.unpack(qL, 6)}, {table.unpack(qR, 6)})
	return rho, vx, vy, vz, hTotal, Cs, Ex, Ey, Ez, Bx, By, Bz
end

function EMHD:calcEigenvalues(vx, Cs)
	local lambdaEuler = {self.euler:calcEigenvalues(vx, Cs)}
	local lambdaMaxwell = {self.maxwell:calcEigenvalues()}
	-- hmm, usually I keep these arranged from smallest to largest .. but that would mean mixing up the eigenvector elements as well...
	-- so I'll just keep them separate: euler, then maxwell
	return table(lambdaEuler):append(lambdaMaxwell):unpack()
end

local function minAndMax(...) return math.min(...), math.max(...) end

function EMHD:calcMinMaxEigenvaluesFromCons(...)
	-- it *should* be the maxwell eigenvalues, because they're the speed of light
	return minAndMax(self:calcEigenvaluesFromCons(...))
end

function EMHD:calcEigenvaluesFromCons(rho, mx, my, mz, EnergyTotal, Ex, Ey, Ez, Bx, By, Bz)
	local rho, vx, vy, vz, P = self.euler:calcPrimFromCons(rho, mx, my, mz, EnergyTotal)
	local hTotal = self.euler:calc_hTotal(rho, P, EnergyTotal)
	local Cs = self.euler:calcSpeedOfSound(vx, vy, vz, hTotal)
	return self:calcEigenvalues(vx, Cs)
end

function EMHD:calcFluxForState(U)
	local fluxEuler = {self.euler:calcFluxForState{table.unpack(U, 1, 5)}}
	local fluxMaxwell = {self.maxwell:calcFluxForState{table.unpack(U, 6)}}
	return table(fluxEuler):append(fluxMaxwell):unpack()
end

function EMHD:fillSubMatrices(dest,a,b)
	for i=1,self.euler.numStates do
		for j=1,self.euler.numStates do
			dest[i][j] = assert(a[i][j])
		end
		for j=self.euler.numStates+1,self.numStates do
			dest[i][j] = 0
		end
	end
	for i=self.euler.numStates+1,self.numStates do
		for j=1,self.euler.numStates do
			dest[i][j] = 0
		end
		for j=self.euler.numStates+1,self.numStates do
			dest[i][j] = assert(b[i-self.euler.numStates][j-self.euler.numStates])
		end
	end
	assert(#dest == self.numStates)
	for i=1,self.numStates do
		assert(#dest[i] == self.numStates, "failed for row "..i.." is "..#dest[i].." should be "..self.numStates)
	end
end

function EMHD:calcEigenBasis(lambda, evR, evL, dF_dU, rho, vx, vy, vz, hTotal, Cs, Ex, Ey, Ez, Bx, By, Bz)
	local lambdaEuler, evrEuler, evlEuler, dFdUEuler = {}, {}, {}, {}
	for i=1,self.euler.numStates do
		evrEuler[i], evlEuler[i], dFdUEuler[i] = {}, {}, {}
	end
	self.euler:calcEigenBasis(lambdaEuler, evrEuler, evlEuler, dFdUEuler, rho, vx, vy, vz, hTotal, Cs)
	
	local lambdaMaxwell, evrMaxwell, evlMaxwell, dFdUMaxwell = {}, {}, {}, {}
	for i=1,self.maxwell.numStates do
		evrMaxwell[i], evlMaxwell[i], dFdUMaxwell[i] = {}, {}, {}
	end
	self.maxwell:calcEigenBasis(lambdaMaxwell, evrMaxwell, evlMaxwell, dFdUMaxwell, Ex, Ey, Ez, Bx, By, Bz)
	
	local lambda = table():append(lambdaEuler):append(lambdaMaxwell)
	if dF_dU then self:fillSubMatrices(dF_dU, dFdUEuler, dFdUMaxwell) end
	if evR then self:fillSubMatrices(evR, evrEuler, evrMaxwell) end
	if evL then self:fillSubMatrices(evL, evlEuler, evlMaxwell) end
end

function EMHD:calcInterfaceEigenvalues(sim, i, qL, qR)
	local rho, vx, vy, vz, hTotal, Cs = self:calcRoeValues(qL, qR)
	fill(sim.eigenvalues[i], self:calcEigenvalues(vx, Cs))
end

function EMHD:calcCellCenterRoeValues(solver, i)
	local rho, vx, vy, vz, hTotal = self.euler:calcCellCenterRoeValues(solver, i)
	local Ex, Ey, Ez, Bx, By, Bz = table.unpack(self.qs[i], 6)
	local Cs = self.euler:calcSpeedOfSound(vx, vy, vz, hTotal)
	return rho, vx, vy, vz, hTotal, Cs, Ex, Ey, Ez, Bx, By, Bz
end

function EMHD:sourceTerm(solver, Us)
do return solver:newState() end
	local m = self.species[1].mass
	local q = self.species[1].charge

	local m0 = 1
	local n0 = 1
	local q0 = 1
	local B0 = 1
	local rL = m0 * self.maxwell.mu0 / (q0 * B0)
	local lambda = math.sqrt(self.maxwell.epsilon0 * self.maxwell.mu0^2 * m0 / (n0 * q0^2))
	local src = solver:newState()
	for i=1,solver.gridsize do
		local rho, mx, my, mz, EnergyTotal, Ex, Ey, Ez, Bx, By, Bz = table.unpack(Us[i])
		local vx, vy, vz = mx/rho, my/rho, mz/rho
		local n = rho / m
		fill(src[i],
			-- Euler 3D
			-- TODO, one of these for each kind of fluid 
			0,	-- rho
			n * q * (Ex + vy * Bz - vz * By),	-- mx
			n * q * (Ey + vz * Bx - vx * Bz),	-- my
			n * q * (Ez + vx * By - vy * Bx), -- mz
			n * q * (vx * Ex + vy * Ey + vz * Ez),	-- EnergyTotal
			-- Maxwell	
			-rL / lambda^2 * n * q * vx,	-- Ex
			-rL / lambda^2 * n * q * vy,	-- Ey
			-rL / lambda^2 * n * q * vz,	-- Ez
			0,	-- Bx
			0,	-- By
			0)	-- Bz
	end
	return src
end

function EMHD:postIterate(solver, Us)
	-- enforce div B = 0 and div E = 1/eps0 sum_a n_a q_a
	-- ... not necessary in 1D?  I know discrete 1D poisson solvers create 1st order solutions .. i.e. solutions that look like c*|x - x0| 
end

return EMHD
