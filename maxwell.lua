-- going by Trangenstein
local class = require 'ext.class'
local table = require 'ext.table'
local Equation = require 'equation'

local Maxwell = class(Equation)

Maxwell.numStates = 6	--E,B xyz

local e0 = 1	-- permittivity
local u0 = 1	-- permissivity

do
	local q = function(self,i) return self.qs[i] end
	local Ex = q:_(1) / e0
	local Ey = q:_(2) / e0
	local Ez = q:_(3) / e0
	local ESq = Ex^2 + Ey^2 + Ez^2
	local E = math.sqrt:o(ESq)
	local Bx = q:_(4)
	local By = q:_(5)
	local Bz = q:_(6)
	local BSq = Bx^2 + By^2 + Bz^2
	local B = math.sqrt:o(BSq)
	local energy = .5 * (ESq * e0 + BSq / u0)
	Maxwell:buildGraphInfos{
		{Ex=Ex}, {Ey=Ey}, {Ez=Ez}, {E=E},
		{Bx=Bx}, {By=By}, {Bz=Bz}, {B=B},
		{['log eigenbasis error'] = function(self,i) return self.eigenbasisErrors and math.log(self.eigenbasisErrors[i]) end},
		{['log reconstruction error'] = function(self,i) return self.fluxMatrixErrors and math.log(self.fluxMatrixErrors[i]) end},
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
	return {Ex * e0, Ey * e0, Ez * e0, Bx, By, Bz}
end

function Maxwell:calcInterfaceEigenBasis(sim,i,qL,qR)
	local avgQ = {}
	for j=1,sim.numStates do
		avgQ[j] = (qL[j] + qR[j]) / 2
	end
	local lambda = 1/sqrt(e0 * u0)
	sim.eigenvalues[i] = {-lambda, -lambda, 0, 0, lambda, lambda}
	local se = sqrt(e0/2)
	local su = sqrt(u0/2)
	local seu = sqrt(e0/u0)/u0
	local ise = 1/se
	local isu = 1/su
	self.fluxTransform = function(self, sim, i, v)
		return {
			0,
			seu * v[6],
			-seu * v[5],
			0,
			-v[3]/seu,
			v[2]/seu
		}
	end
	self.eigenfields = function(self, sim, i, v)
		return {
			ise * v[3] + isu * v[5],
			-ise * v[2] + isu * v[6],
			-ise * v[1] + isu * v[4],
			ise * v[1] + isu * v[4],
			ise * v[2] + isu * v[6],
			-ise * v[3] + isu * v[5]
		}
	end
	self.eigenfieldsInverse = function(self, sim, i, v)
		return {
			-se * v[3] + se * v[4],
			-se * v[2] + se * v[5],
			se * v[1] - se * v[6],
			su * v[3] + su * v[4],
			su * v[1] + su * v[6],
			su * v[2] + su * v[5]
		}
	end
end

local sigma = 1	-- conductivity
function Maxwell:sourceTerm(sim, qs)
	local source = sim:newState()
	for i=1,sim.gridsize do
		local eEx, eEy, eEz, Bx, By, Bz = unpack(qs[i])
		local Ex = eEx / e0
		local Ey = eEy / e0
		local Ez = eEz / e0
		source[i][1] = -Ex * sigma
		source[i][2] = -Ey * sigma
		source[i][3] = -Ez * sigma
	end
	return source
end

return Maxwell
