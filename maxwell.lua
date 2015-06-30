-- going by Trangenstein
local class = require 'ext.class'
local table = require 'ext.table'

local Maxwell = class()

Maxwell.numStates = 6	--E,B xyz

local e0 = 1	-- permittivity
local u0 = 1	-- permissivity

function Maxwell:init(...)
	Maxwell.super.init(self, ...)
	local get_state = function(self,i) return self.qs[i] end
	local Ex = function(self,i) return self.qs[i][1] / e0 end
	local Ey = function(self,i) return self.qs[i][2] / e0 end
	local Ez = function(self,i) return self.qs[i][3] / e0 end
	local ESq = Ex^2 + Ey^2 + Ez^2
	local Bx = function(self,i) return self.qs[i][4] end
	local By = function(self,i) return self.qs[i][5] end
	local Bz = function(self,i) return self.qs[i][6] end
	local BSq = Bx^2 + By^2 + Bz^2
	self.graphInfos = table{
		{viewport={0/4, 0/3, 1/4, 1/3}, getter=Ex, name='Ex', color={1,0,0}},
		{viewport={1/4, 0/3, 1/4, 1/3}, getter=Ey, name='Ey', color={0,1,0}},
		{viewport={2/4, 0/3, 1/4, 1/3}, getter=Ez, name='Ez', color={0,0,1}},
		{viewport={3/4, 0/3, 1/4, 1/3}, getter=sqrt:compose(ESq), name='E', color={0,0,1}},
		{viewport={0/4, 1/3, 1/4, 1/3}, getter=Bx, name='Bx', color={0,1,1}},
		{viewport={1/4, 1/3, 1/4, 1/3}, getter=By, name='By', color={1,0,1}},
		{viewport={2/4, 1/3, 1/4, 1/3}, getter=Bz, name='Bz', color={1,1,0}},
		{viewport={3/4, 1/3, 1/4, 1/3}, getter=sqrt:compose(BSq), name='B', color={1,1,0}},
		{viewport={0/4, 2/3, 1/4, 1/3}, getter=function(self,i) return math.log(self.eigenbasisErrors[i]) end, name='log eigenbasis error', color={1,0,0}, range={-30, 30}},
		{viewport={1/4, 2/3, 1/4, 1/3}, getter=function(self,i) return math.log(self.fluxMatrixErrors[i]) end, name='log reconstruction error', color={1,0,0}, range={-30, 30}},
		{viewport={3/4, 2/3, 1/4, 1/3}, getter=.5 * (ESq * e0 + BSq / u0), name='energy', color={0,.5,1}},
	}
	self.graphInfoForNames = self.graphInfos:map(function(info,i)
		return info, info.name
	end)
end

function Maxwell:initCell(i)
	local x = self.xs[i]
	local Ex = 0
	local Ey = 0
	local Ez = 1
	local Bx = 1
	local By = x < 0 and 1 or -1
	local Bz = 0
	return {Ex * e0, Ey * e0, Ez * e0, Bx, By, Bz}
end

function Maxwell:calcInterfaceEigenBasis(i,qL,qR)
	local avgQ = {}
	for j=1,self.numStates do
		avgQ[j] = (qL[j] + qR[j]) / 2
	end
	local lambda = 1/sqrt(e0 * u0)
	self.eigenvalues[i] = {-lambda, -lambda, 0, 0, lambda, lambda}
	local se = sqrt(e0/2)
	local su = sqrt(u0/2)
	local seu = sqrt(e0/u0)/u0
	local ise = 1/se
	local isu = 1/su
	self.fluxTransform = function(self, i, v)
		return {
			0,
			seu * v[6],
			-seu * v[5],
			0,
			-v[3]/seu,
			v[2]/seu
		}
	end
	self.eigenfields = function(self, i, v)
		return {
			ise * v[3] + isu * v[5],
			-ise * v[2] + isu * v[6],
			-ise * v[1] + isu * v[4],
			ise * v[1] + isu * v[4],
			ise * v[2] + isu * v[6],
			-ise * v[3] + isu * v[5]
		}
	end
	self.eigenfieldsInverse = function(self, i, v)
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
function Maxwell:addSourceToDerivCell(dq_dts, i)
	local eEx, eEy, eEz, Bx, By, Bz = unpack(self.qs[i])
	local Ex = eEx / e0
	local Ey = eEy / e0
	local Ez = eEz / e0
	dq_dts[i][1] = dq_dts[i][1] - Ex * sigma
	dq_dts[i][2] = dq_dts[i][2] - Ey * sigma
	dq_dts[i][3] = dq_dts[i][3] - Ez * sigma
end

return Maxwell

