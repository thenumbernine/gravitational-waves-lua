local class = require 'ext.class'
local Equation = require 'equation'

local Euler3D = class(Equation)
Euler3D.name = 'Euler 3D'
Euler3D.numStates = 5
Euler3D.gamma = 5/3

do
	local q = function(self, i) return self.qs[i] end
	local gamma = function(self) return self.equation.gamma end
	local rho = q:_(1)
	local mx, my, mz = q:_(2), q:_(3), q:_(4)
	local ETotal = q:_(5)
	local vx, vy, vz = mx/rho, my/rho, mz/rho
	local eKin = .5 * (vx*vx + vy*vy + vz*vz)
	local EKin = rho * eKin
	local EInt = ETotal - EKin
	local eInt = EInt / rho
	local P = (gamma - 1) * EInt
	local S = P / rho^gamma
	local H = EInt + P 
	local h = H / rho 
	local HTotal = ETotal + P 
	local hTotal = HTotal / rho

	Euler3D:buildGraphInfos{
		-- prims
		{rho = rho},
		{vx = vx},
		{vy = vy},
		{vz = vz},
		{P = P},
		-- other vars
		{EInt = EInt},
		{EKin = EKin},
		{H = H},
		{S = S},
		-- conservative
		--{mx = mx},
		{ETotal = ETotal},
	}
end

function Euler3D:initCell(sim,i)
	local x = sim.xs[i]
	local rho = x < 0 and 1 or .125
	local vx, vy, vz = 0, 0, 0
	local P = x < 0 and 1 or .1
	local q = {self:calcConsFromPrim(rho, vx, vy, vz, P)}
	return q
end

function Euler3D:calcConsFromPrim(rho, vx, vy, vz, P)
	local EInt = P/(self.gamma-1)
	local EKin = .5 * rho * (vx*vx + vy*vy + vz*vz)
	local ETotal = EInt + EKin
	return 
		rho,
		rho * vx,
		rho * vy,
		rho * vz,
		ETotal
end

function Euler3D:calcPrimFromCons(rho, mx, my, mz, ETotal)
	local vx, vy, vz = mx/rho, my/rho, mz/rho
	local EKin = .5 * rho * (vx*vx + vy*vy + vz*vz)
	local EInt = ETotal - EKin
	local P = (self.gamma - 1) * EInt
	return rho, vx, vy, vz, P
end

function Euler3D:calcRoeValues(qL, qR)
	local gamma = self.gamma

	local ETotalL = qL[5]
	local rhoL, vxL, vyL, vzL, PL = self:calcPrimFromCons(table.unpack(qL))
	local hTotalL = self:calc_hTotal(rhoL, PL, ETotalL)
	local sqrtRhoL = math.sqrt(rhoL)
	
	local ETotalR = qR[5]
	local rhoR, vxR, vyR, vzR, PR = self:calcPrimFromCons(table.unpack(qR))
	local hTotalR = self:calc_hTotal(rhoR, PR, ETotalR)
	local sqrtRhoR = math.sqrt(rhoR)
	
	local rho = sqrtRhoL * sqrtRhoR
	local vx = (sqrtRhoL * vxL + sqrtRhoR * vxR) / (sqrtRhoL + sqrtRhoR)
	local vy = (sqrtRhoL * vyL + sqrtRhoR * vyR) / (sqrtRhoL + sqrtRhoR)
	local vz = (sqrtRhoL * vzL + sqrtRhoR * vzR) / (sqrtRhoL + sqrtRhoR)
	local hTotal = (sqrtRhoL * hTotalL + sqrtRhoR * hTotalR) / (sqrtRhoL + sqrtRhoR)
	local Cs = self:calcSpeedOfSound(vx, vy, vz, hTotal)

	return rho, vx, vy, vz, hTotal, Cs
end

function Euler3D:calcEigenvalues(vx, Cs)
	return vx - Cs, vx, vx, vx, vx + Cs
end

function Euler3D:calcSpeedOfSound(vx, vy, vz, hTotal)
	local eKin = .5 * (vx*vx + vy*vy + vz*vz)
	return math.sqrt((self.gamma - 1) * (hTotal - eKin))
end

function Euler3D:calc_hTotal(rho, P, ETotal)
	return (ETotal + P) / rho
end

function Euler3D:calcEigenvaluesFromCons(rho, mx, my, mz, ETotal)
	local rho, vx, vy, vz, P = self:calcPrimFromCons(rho, mx, my, mz, ETotal)
	local hTotal = self:calc_hTotal(rho, P, ETotal)
	local Cs = self:calcSpeedOfSound(vx, vy, vz, hTotal)
	return self:calcEigenvalues(vx, Cs)
end

-- x-direction flux
function Euler3D:calcFluxForState(q)
	local gamma = self.gamma
	local rho, mx, my, mz, ETotal = table.unpack(q)
	local rho, vx, vy, vz, P = self:calcPrimFromCons(table.unpack(q))
	return 
		mx,
		vx * mx + P,
		vy * mx,
		vz * mx,
		(ETotal + P) * vx
end

function Euler3D:calcEigenBasis(lambda, evR, evL, dF_dU, rho, vx, vy, vz, hTotal, Cs)
	Cs = Cs or self:calcSpeedOfSound(vx, vy, vz, hTotal)
	
	local gamma = self.gamma
	local gamma_1 = gamma - 1
	local gamma_3 = gamma - 3
	local CsSq = Cs * Cs
	local vSq = vx*vx + vy*vy + vz*vz

	fill(lambda, self:calcEigenvalues(vx, Cs))

	if dF_dU then
		fill(dF_dU[1], 0, 									1, 								0, 					0, 					0			)
		fill(dF_dU[2], -vx * vx + .5 * gamma_1 * vSq,		-vx * gamma_3,					-vy * gamma_1,		-vz * gamma_1,		gamma - 1	)
		fill(dF_dU[3], -vx * vy,							vy, 							vx, 				0, 					0			)
		fill(dF_dU[4], -vx * vz, 							vz, 							0, 					vx, 				0			)
		fill(dF_dU[5], 	vx * (.5 * vSq * gamma_1 - hTotal),	-gamma_1 * vx * vx + hTotal,	-gamma_1 * vx * vy,	-gamma_1 * vx * vz,	gamma * vx	)
	end

	fill(evR[1], 1, 				1, 			0, 	0, 	1				)
	fill(evR[2], vx - Cs, 			vx, 		0, 	0, 	vx + Cs			)
	fill(evR[3], vy, 				vy, 		1, 	0, 	vy				)
	fill(evR[4], vz, 				vz, 		0, 	1, 	vz				)
	fill(evR[5], hTotal - Cs * vx, .5 * vSq, 	vy, vz, hTotal + Cs * vx)

	local invDenom = .5 / CsSq
	fill(evL[1], (.5 * gamma_1 * vSq + Cs * vx) * invDenom,	-(Cs + gamma_1 * vx) * invDenom,	-gamma_1 * vy * invDenom,		-gamma_1 * vz * invDenom,		gamma_1 * invDenom		)
	fill(evL[2], 1 - gamma_1 * vSq * invDenom,				gamma_1 * vx * 2 * invDenom,		gamma_1 * vy * 2 * invDenom,	gamma_1 * vz * 2 * invDenom,	-gamma_1 * 2 * invDenom	)
	fill(evL[3], -vy, 										0, 									1, 								0, 								0						)
	fill(evL[4], -vz, 										0, 									0, 								1, 								0						)
	fill(evL[5], (.5 * gamma_1 * vSq - Cs * vx) * invDenom,	(Cs - gamma_1 * vx) * invDenom,		-gamma_1 * vy * invDenom,		-gamma_1 * vz * invDenom,		gamma_1 * invDenom		)
end

function Euler3D:calcInterfaceEigenvalues(sim, i, qL, qR)
	local rho, vx, vy, vz, hTotal, Cs = self:calcRoeValues(qL, qR)
	fill(sim.eigenvalues[i], self:calcEigenvalues(vx, Cs))
end

function Euler3D:calcCellCenterRoeValues(solver, i)
	local rho, mx, my, mz, ETotal = table.unpack(solver.qs[i])
	local rho, vx, vy, vz, P = self:calcPrimFromCons(rho, mx, my, mz, ETotal)
	local hTotal = self:calc_hTotal(rho, P, ETotal)
	return rho, vx, vy, vz, hTotal
end

return Euler3D
