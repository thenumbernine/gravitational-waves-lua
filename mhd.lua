require 'ext'
local Simulation = require 'simulation'
local MHDSimulation = class(Simulation)

MHDSimulation.numStates = 8
MHDSimulation.gamma = 5/3	

function MHDSimulation:init(...)
	Simulation.init(self, ...)
	
	local getState = index:bind(self.qs)
	self.graphInfos = {
		{viewport={0/3, 0/2, 1/3, 1/2}, getter=getState:index(1), name='rho', color={1,0,1}},
		{viewport={1/3, 0/2, 1/3, 1/2}, getter=getState:index(2) / getState:index(1), name='u', color={0,1,0}},
		{viewport={2/3, 0/2, 1/3, 1/2}, getter=getState:index(3) / getState:index(1), name='E', color={.5,.5,1}},
		{viewport={0/3, 1/2, 1/3, 1/2}, getter=log:compose(index:bind(self.eigenbasisErrors)), name='log eigenbasis error', color={1,0,0}, range={-30, 30}},
		{viewport={1/3, 1/2, 1/3, 1/2}, getter=log:compose(index:bind(self.fluxMatrixErrors)), name='log reconstruction error', color={1,0,0}, range={-30, 30}},
	}
end

MHDSimulation.mu = 1

function MHDSimulation:initCell(i)
	local rho = self.xs[i] < 0 and .1 or 1
	local ux, uy, uz = 0, 0, 0
	local Bx, By, Bz = 0, 0, 0
	local eInt = 1
	local eKin = .5 * (ux^2 + uy^2 + uz^2)
	local EMag = .5 * (Bx^2 + By^2 + Bz^2) / self.mu
	local ETotal = rho * (eInt + eKin) + EMag 
	return {rho, rho * ux, rho * uy, rho * uz, Bx, By, Bz, ETotal}
end

function MHDSimulation:stateToPrims(rho, mx, my, mz, Bx, By, Bz, ETotal)
	local ux, uy, uz = mx / rho, my / rho, mz / rho
	local EMag = .5 * (Bx^2 + By^2 + Bz^2) / self.mu
	local EHydro = ETotal - EMag
	local eHydro = EHydro / rho
	local eKin = .5 * (ux^2 + uy^2 + uz^2)
	local eInt = eHydro - eKin
	local p = eInt / (self.gamma - 1)
	return rho, ux, uy, uz, Bx, By, Bz, p
end

function MHDSimulation:calcInterfaceEigenBasis(i)
	local gamma = self.gamma
	local gammaMinusOne = gamma - 1

	local rhoL, uxL, uyL, uzL, BxL, ByL, BzL, pL = self:stateToPrims(unpack(self.qs[i-1]))
	local rhoR, uxR, uyR, uzR, BxR, ByR, BzR, pR = self:stateToPrims(unpack(self.qs[i]))

	local rho = .5 * (rhoL + rhoR)
	local ux = .5 * (uxL + uxR)
	local uy = .5 * (uyL + uyR)
	local uz = .5 * (uzL + uzR)
	local Bx = .5 * (BxL + BxR)
	local By = .5 * (ByL + ByR)
	local Bz = .5 * (BzL + BzR)
	local p = .5 * (pL + pR)

	local mu = self.mu
	local bSq = Bx*Bx + By*By + Bz*Bz
	local vaxSq = Bx*Bx / (mu * rho)
	local vax = sqrt(vaxSq)
	local CsSq = gamma * p / rho
	local Cs= sqrt(CsSq)
	local vaSq = bSq / (mu * rho)
	local va = sqrt(vaSq)
	local cStarSq = .5 * (vaSq + CsSq)
	local discr = sqrt(cStarSq * cStarSq - vaxSq * CsSq)
	local vfSq = cStarSq + discr
	local vsSq = cStarSq - discr
	local vf = sqrt(vfSq)
	local vs = sqrt(vsSq)
	local sgnBx = Bx >= 0 and 1 or -1

	self.fluxMatrix[i] = {
		{0,0,0,0,0,0,0,0},
		{0,0,0,0,0,0,0,0},
		{0,0,0,0,0,0,0,0},
		{0,0,0,0,0,0,0,0},
		{0,0,0,0,0,0,0,0},
		{0,0,0,0,0,0,0,0},
		{0,0,0,0,0,0,0,0},
		{0,0,0,0,0,0,0,0}
	}
	
	local S = self.eigenvalues[i]
	S[1] = ux - vf
	S[2] = ux - vax
	S[3] = ux - vs
	S[4] = ux
	S[5] = ux
	S[6] = ux + vs
	S[7] = ux + vax
	S[8] = ux + vf

	local tau = 1 / rho
	
	local afSq = (vfSq - vaxSq) / (vfSq - vsSq)
	local af
	if afSq > 0 then
		af = sqrt(afSq)
	else
		af, afSq = 1, 1
	end
	
	local asSq = (vfSq - CsSq) / (vfSq - vsSq)
	local as
	if asSq > 0 then
		as = sqrt(asSq)
	else
		as, asSq = 1, 1
	end

	local BTSq = By * By + Bz * Bz
	local betay = 1
	local betaz = 0
	if BTSq > 0 then
		local BT = sqrt(BTSq)
		betay = By / BT
		betaz = Bz / BT
	end
		
	local RFast = vf / sqrt(afSq * (vfSq + CsSq) + asSq * (vfSq + vaxSq))
	local RSlow = vfSq / sqrt(afSq * CsSq * (vfSq + CsSq) + asSq * vfSq * (vsSq + CsSq))

	local R_A = {	--right eigenvectors
	--fast -
		{
			af * tau * RFast,
			af * vf * RFast,
			-as * betay * vax * sgnBx * RFast,
			-as * betaz * vax * sgnBx * RFast,
			0,
			-as * betay * vf * sqrt(mu * rho) * RFast,
			-as * betaz * vf * sqrt(mu * rho) * RFast,
			-af * gamma * p * RFast
		},
	--alfven -
		{
			0,
			0,
			-betaz * vf / sqrt(2),
			betay * vf / sqrt(2),
			0,
			-sgnBx * sqrt(mu * rho) * betaz * vf / sqrt(2),
			sgnBx * sqrt(mu * rho) * betay * vf / sqrt(2),
			0
		},
	--slow -
		{
			as * tau,
			as * vs,
			af * betay * Cs * sgnBx,
			af * betaz * Cs * sgnBx,
			0,
			af * betay * CsSq / vf * sqrt(mu * rho),
			af * betaz * CsSq / vf * sqrt(mu * rho),
			-as * gamma * p
		},
	--entropy
		{tau, 0, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 1, 0, 0, 0},
	--zero
	--slow +
		{
			-as * tau,
			as * vs,
			af * betay * Cs * sgnBx,
			af * betaz * Cs * sgnBx,
			0,
			-af * betay * CsSq / vf * sqrt(mu * rho),
			-af * betaz * CsSq / vf * sqrt(mu * rho),
			as * gamma * p
		},
	--alfven +
		{
			0,
			0,
			betaz * vf / sqrt(2),
			-betay * vf / sqrt(2),
			0,
			-sgnBx * sqrt(mu * rho) * betaz * vf / sqrt(2),
			sgnBx * sqrt(mu * rho) * betay * vf / sqrt(2),
			0
		},
	--fast +
		{
			-af * tau * RFast,
			af * vf * RFast,
			-as * betay * vax * sgnBx * RFast,
			-as * betaz * vax * sgnBx * RFast,
			0,
			as * betay * vf * sqrt(mu * rho) * RFast,
			as * betaz * vf * sqrt(mu * rho) * RFast,
			af * gamma * p * RFast
		}
	}

	--specify in rows (so it will look as-is)
	--so L_ij = L_A[i][j]
	local L_A = {
	--fast -
		{
			0,
			 af * vf * RFast / vfSq,
			 -as * betay * vax * sgnBx * RFast / vfSq,
			 -as * betaz * vax * sgnBx * RFast / vfSq,
			 0,
			 -as * betay * vf / sqrt(mu * rho) * RFast / vfSq,
			 -as * betaz * vf / sqrt(mu * rho) * RFast / vfSq,
			 -af * tau * RFast / vfSq
		},
	--alfven -
		{
			0,
			0,
			-betaz / (vf * sqrt(2)),
			betay / (vf * sqrt(2)),
			0,
			-sgnBx * betaz / sqrt(mu * rho) / (vf * sqrt(2)),
			sgnBx * betay / sqrt(mu * rho) / (vf * sqrt(2)),
			0
		},
	--slow -
		{
			0,
			as * vs * RSlow / vfSq,
			af * betay * Cs * sgnBx * RSlow / vfSq,
			af * betaz * Cs * sgnBx * RSlow / vfSq,
			0,
			af * betay * CsSq / (sqrt(mu * rho) * vf) * RSlow / vfSq,
			af * betaz * CsSq / (sqrt(mu * rho) * vf) * RSlow / vfSq,
			-as * gamma * p * RSlow / vfSq
		},
	--entropy
		{rho, 0, 0, 0, 0, 0, 0, 1 / (gamma * p)},
	--zero
		{0, 0, 0, 0, 1, 0, 0, 0},
	--slow +
		{
			0,
			as * vs * RSlow / vfSq,
			af * betay * Cs * sgnBx * RSlow / vfSq,
			af * betaz * Cs * sgnBx * RSlow / vfSq,
			0,
			-af * betay * CsSq / (sqrt(mu * rho) * vf) * RSlow / vfSq,
			-af * betaz * CsSq / (sqrt(mu * rho) * vf) * RSlow / vfSq,
			as * gamma * p * RSlow / vfSq
		},
	--alfven +
		{
			0,
			0, 
			betaz / (vf * sqrt(2)),
			-betay / (vf * sqrt(2)),
			0, 
			-sgnBx * betaz / sqrt(mu * rho) / (vf * sqrt(2)),
			sgnBx * betay / sqrt(mu * rho) / (vf * sqrt(2)),
			0
		},
	--fast +
		{
			0,
			 af * vf * RFast / vfSq,
			 -as * betay * vax * sgnBx * RFast / vfSq,
			 -as * betaz * vax * sgnBx * RFast / vfSq,
			 0,
			 as * betay * vf / sqrt(mu * rho) * RFast / vfSq,
			 as * betaz * vf / sqrt(mu * rho) * RFast / vfSq,
			 af * tau * RFast / vfSq
		}
	}

	local tauSq = tau * tau
	local uSq = ux * ux + uy * uy + uz * uz

	--specified by row, so dU_dW[i][j] == (dU/dW)_ij
	local dU_dW = {
		{-rho*rho,	0,	0,	0,	0,	0,	0,	0},
		{-ux*rho*rho,	rho,	0,	0,	0,	0,	0,	0},
		{-uy*rho*rho,	0,	rho,	0,	0,	0,	0,	0},
		{-uz*rho*rho,	0,	0,	rho,	0,	0,	0,	0},
		{0,	0,	0,	0,	1,	0,	0,	0},
		{0,	0,	0,	0,	0,	1,	0,	0},
		{0,	0,	0,	0,	0,	0,	1,	0},
		{-.5*uSq*rho*rho,	ux*rho,	uy*rho,	uz*rho,	Bx/mu,	By/mu,	Bz/mu,	1/gammaMinusOne}
	}

	--specified by row, so dW_dU[i][j] == (dW/dU)_ij
	local dW_dU = {
		{-tauSq,	0,	0,	0,	0,	0,	0,	0},
		{-tau*ux,	tau,	0,	0,	0,	0,	0,	0},
		{-tau*uy,	0,	tau,	0,	0,	0,	0,	0},
		{-tau*uz,	0,	0,	tau,	0,	0,	0,	0},
		{0,	0,	0,	0,	1,	0,	0,	0},
		{0,	0,	0,	0,	0,	1,	0,	0},
		{0,	0,	0,	0,	0,	0,	1,	0},
		{.5*uSq * gammaMinusOne,	ux * gammaMinusOne,	uy * gammaMinusOne,	uz * gammaMinusOne,	-Bx/mu * gammaMinusOne,	-By/mu * gammaMinusOne,	-Bz/mu * gammaMinusOne,	gammaMinusOne}
	}

	--now transform these to the left and right eigenvectors of the flux ...
	--with transformations: R_U = dU/dW * R_A and L_U = L_A * dW/dU
	--don't forget indexing is A_ij == A[i][j] except R_A is transposed
	for j=1,8 do
		for k=1,8 do
			local sum = 0
			for m=1,8 do
				sum = sum + dU_dW[j][m] * R_A[k][m]
			end
			self.eigenvectors[i][j][k] = sum
		end
		for k=1,8 do
			local sum = 0
			for m=1,8 do
				sum = sum + L_A[j][m] * dW_dU[m][k]
			end
			self.eigenvectorsInverse[i][j][k] = sum
		end
	end
end

return MHDSimulation

