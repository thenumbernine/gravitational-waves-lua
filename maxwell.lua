-- going by Trangenstein
local class = require 'ext.class'
local table = require 'ext.table'
local Equation = require 'equation'

local Maxwell = class(Equation)

Maxwell.numStates = 6	--E,B xyz

Maxwell.epsilon0 = 1	-- permittivity
Maxwell.mu0 = 1	-- permissivity
Maxwell.sigma = 1	-- conductivity

do
	local q = function(self,i) return self.qs[i] end
	local Ex = q:_(1) / Maxwell.epsilon0
	local Ey = q:_(2) / Maxwell.epsilon0
	local Ez = q:_(3) / Maxwell.epsilon0
	local ESq = Ex^2 + Ey^2 + Ez^2
	local E = math.sqrt:o(ESq)
	local Bx = q:_(4)
	local By = q:_(5)
	local Bz = q:_(6)
	local BSq = Bx^2 + By^2 + Bz^2
	local B = math.sqrt:o(BSq)
	local energy = .5 * (ESq * Maxwell.epsilon0 + BSq / Maxwell.mu0)
	Maxwell:buildGraphInfos{
		{Ex=Ex}, {Ey=Ey}, {Ez=Ez}, {E=E},
		{Bx=Bx}, {By=By}, {Bz=Bz}, {B=B},
		{energy = energy},
	}
end
Maxwell.graphInfoForNames = Maxwell.graphInfos:map(function(info,i)
	return info, info.name
end)

function Maxwell:initCell(sim, i)
	local x = sim.xs[i]
	local Ex = 0
	local Ey = 0
	local Ez = 1
	local Bx = 1
	local By = x < 0 and 1 or -1
	local Bz = 0
	return {Ex * self.epsilon0, Ey * self.epsilon0, Ez * self.epsilon0, Bx, By, Bz}
end

function Maxwell:calcFluxForState(q)
	local seu = math.sqrt(self.epsilon0/self.mu0)/self.mu0
	return
		0,
		seu * q[6],
		-seu * q[5],
		0,
		-q[3]/seu,
		q[2]/seu
end

-- return nothing.  this is passed on to 'calcEigenBasis', which does nothing and stores nothing.
-- eigenvalues/vectors
function Maxwell:calcInterfaceRoeValues(solver, i)
end

-- typically stores eigenvalues, left and right eigenvectors, flux matrix
-- but no information is needed to recreate any of those, so don't store any
-- except for the values, which are used externally
function Maxwell:calcEigenBasis(lambda, evr, evl, dF_dU)
	fill(lambda, self:calcEigenvalues())

	-- .. but for orthogonality error's sake I'll fill them anyways
	-- (that and, now, EMHD, which references this)
	-- If you wanted to bureaucratize this even further, you could make a 'calcEigenBasis' and a 'getEigenBasis' 
	-- one to calculate what's needed upon update (which for Maxwell is nothing) and one for getting the actual eigenbasis (which EMHD needs)
	if dF_dU then
		fill(dF_dU[1], 0, 0, 0, 0, 0, 0)
		fill(dF_dU[2], 0, 0, 0, 0, 0, 1/self.mu0)
		fill(dF_dU[3], 0, 0, 0, 0, -1/self.mu0, 0)
		fill(dF_dU[4], 0, 0, 0, 0, 0, 0)
		fill(dF_dU[5], 0, 0, -1/self.epsilon0, 0, 0, 0)
		fill(dF_dU[6], 0, 1/self.epsilon0, 0, 0, 0, 0)
	end

	local _2 = 1/math.sqrt(2)
	local ue = math.sqrt(self.mu0 / (2 * self.epsilon0))
	local eu = math.sqrt(self.epsilon0 / (2 * self.mu0))

	fill(evr[1],	0,	0,	1,	0,	0,	0)
	fill(evr[2],	_2,	0,	0,	0,	_2,	0)
	fill(evr[3],	0,	_2,	0,	0,	0,	_2)
	fill(evr[4],	0,	0,	0,	1,	0,	0)
	fill(evr[5],	0,	ue,	0,	0,	0,	-ue)
	fill(evr[6],	-ue,0,	0,	0,	ue,	0)

	fill(evl[1],	0,	_2,	0,	0,	0,	-eu)
	fill(evl[2],	0,	0,	_2,	0,	eu,	0)
	fill(evl[3],	1,	0,	0,	0,	0,	0)
	fill(evl[4],	0,	0,	0,	1,	0,	0)
	fill(evl[5],	0,	_2,	0,	0,	0,	eu)
	fill(evl[6],	0,	0,	_2,	0,	-eu,0)
end

function Maxwell:calcInterfaceEigenvalues(sim,i,qL,qR)
	fill(sim.eigenvalues[i], self:calcEigenvalues())
end

function Maxwell:applyFluxMatrix(sim, i, v)
	return {self:calcFluxForState(v)}
end

function Maxwell:eigenLeftTransform(solver, m, v)
	local se = math.sqrt(self.epsilon0/2)
	local ise = 1/se
	local su = math.sqrt(self.mu0/2)
	local isu = 1/su
	return {
		ise * v[3] + isu * v[5],
		-ise * v[2] + isu * v[6],
		-ise * v[1] + isu * v[4],
		ise * v[1] + isu * v[4],
		ise * v[2] + isu * v[6],
		-ise * v[3] + isu * v[5]
	}
end

function Maxwell:eigenRightTransform(solver, m, v)
	local se = math.sqrt(self.epsilon0/2)
	local su = math.sqrt(self.mu0/2)
	return {
		-se * v[3] + se * v[4],
		-se * v[2] + se * v[5],
		se * v[1] - se * v[6],
		su * v[3] + su * v[4],
		su * v[1] + su * v[6],
		su * v[2] + su * v[5]
	}
end

function Maxwell:calcMaxEigenvalue()
	local lambda = 1/math.sqrt(self.epsilon0 * self.mu0)
	return lambda
end

function Maxwell:calcEigenvalues(...)
	local lambda = self:calcMaxEigenvalue()
	return -lambda, -lambda, 0, 0, lambda, lambda
end
Maxwell.calcEigenvaluesFromCons = Maxwell.calcEigenvalues

function Maxwell:calcMinMaxEigenvaluesFromCons(...)
	local lambda = self:calcMaxEigenvalue()
	return -lambda, lambda
end

function Maxwell:sourceTerm(sim, qs)
	local source = sim:newState()
	for i=1,sim.gridsize do
		local eEx, eEy, eEz, Bx, By, Bz = unpack(qs[i])
		local Ex = eEx / self.epsilon0
		local Ey = eEy / self.epsilon0
		local Ez = eEz / self.epsilon0
		source[i][1] = -Ex * self.sigma
		source[i][2] = -Ey * self.sigma
		source[i][3] = -Ez * self.sigma
	end
	return source
end

return Maxwell
