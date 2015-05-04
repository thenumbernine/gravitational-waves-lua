require 'ext'
local Simulation = require 'simulation'
local MHDSimulation = class(Simulation)

MHDSimulation.numStates = 8
MHDSimulation.gamma = 5/3	
MHDSimulation.mu = 1

function MHDSimulation:init(...)
	Simulation.init(self, ...)

	local function initArray()
		local t = {}
		for i=1,self.gridsize do
			t[i] = 0
		end
		return t
	end
	self.primOrthoError = initArray()	-- ortho error of dU/dW * dW/dU
	self.consOrthoError = initArray()	-- ortho error of L_A * R_A

	local mu = self.mu
	local gamma = self.gamma
	local getState = index:bind(self.qs)
	local rho = getState:index(1)
	local ux = getState:index(2) / rho
	local uy = getState:index(3) / rho
	local uz = getState:index(4) / rho
	local Bx = getState:index(5) / mu
	local By = getState:index(6) / mu
	local Bz = getState:index(7) / mu
	local ETotal = getState:index(8)
	local eTotal = ETotal / rho
	local EMag = .5*(Bx*Bx + By*By + Bz*Bz) / mu
	local EHydro = ETotal - EMag
	local eHydro = EHydro / rho
	local eKin = .5*(ux*ux + uy*uy + uz*uz)
	local eInt = eHydro - eKin
	local p = eInt / (gamma - 1)
	self.graphInfos = {
		{viewport={0/5, 0/4, 1/5, 1/4}, getter=(index:bind(self.eigenbasisErrors)), name='eigenbasis error', color={1,0,0}, range={-30, 30}},
		{viewport={0/5, 1/4, 1/5, 1/4}, getter=(index:bind(self.fluxMatrixErrors)), name='reconstruction error', color={1,0,0}, range={-30, 30}},
		{viewport={0/5, 2/4, 1/5, 1/4}, getter=(index:bind(self.primOrthoError)), name='dU/dW ortho error', color={1,0,0}, range={-30, 30}},
		{viewport={0/5, 3/4, 1/5, 1/4}, getter=(index:bind(self.consOrthoError)), name='dF/dU eigenbasis error', color={1,0,0}, range={-30, 30}},
		{viewport={1/5, 0/4, 1/5, 1/4}, getter=rho, name='rho', color={1,0,1}},
		{viewport={1/5, 1/4, 1/5, 1/4}, getter=rho, name='p', color={1,0,1}},
		{viewport={2/5, 0/4, 1/5, 1/4}, getter=ux, name='ux', color={0,1,0}},
		{viewport={2/5, 1/4, 1/5, 1/4}, getter=uy, name='uy', color={0,1,0}},
		{viewport={2/5, 2/4, 1/5, 1/4}, getter=uz, name='uz', color={0,1,0}},
		{viewport={3/5, 0/4, 1/5, 1/4}, getter=Bx, name='Bx', color={.5,.5,1}},
		{viewport={3/5, 1/4, 1/5, 1/4}, getter=By, name='By', color={.5,.5,1}},
		{viewport={3/5, 2/4, 1/5, 1/4}, getter=Bz, name='Bz', color={.5,.5,1}},
		{viewport={4/5, 0/4, 1/5, 1/4}, getter=eTotal, name='eTotal', color={1,1,0}},
		{viewport={4/5, 1/4, 1/5, 1/4}, getter=eKin, name='eKin', color={1,1,0}},
		{viewport={4/5, 2/4, 1/5, 1/4}, getter=eInt, name='eInt', color={1,1,0}},
		{viewport={4/5, 3/4, 1/5, 1/4}, getter=eInt, name='eHydro', color={1,1,0}},
		{viewport={3/5, 3/4, 1/5, 1/4}, getter=eInt, name='EMag', color={1,1,0}},
	}
end

function MHDSimulation:initCell(i)
	local x = self.xs[i]
	local rho = x < 0 and .1 or 1
	local ux, uy, uz = 0, 0, 0
	local Bx = 0	-- .75
	local By = 0	-- x < 0 and 1 or -1
	local Bz = 0
	local p = 1		-- x < 0 and 1 or .1
	local eInt = p / (self.gamma - 1)
	local eKin = .5*(ux*ux + uy*uy + uz*uz)
	local EMag = .5*(Bx*Bx + By*By + Bz*Bz) / self.mu
	local ETotal = rho*(eInt + eKin) + EMag 
	return {rho, rho*ux, rho*uy, rho*uz, Bx, By, Bz, ETotal}
end

function MHDSimulation:stateToPrims(rho, mx, my, mz, Bx, By, Bz, ETotal)
	local ux, uy, uz = mx / rho, my / rho, mz / rho
	local EMag = .5*(Bx*Bx + By*By + Bz*Bz) / self.mu
	local EHydro = ETotal - EMag
	local eHydro = EHydro / rho
	local eKin = .5*(ux*ux + uy*uy + uz*uz)
	local eInt = eHydro - eKin
	local p = eInt / (self.gamma - 1)
	return rho, ux, uy, uz, Bx, By, Bz, p
end

function MHDSimulation:calcInterfaceEigenBasis(i)
	local gamma = self.gamma
	local gammaMinusOne = gamma - 1

	local rhoL, uxL, uyL, uzL, BxL, ByL, BzL, pL = self:stateToPrims(unpack(self.qs[i-1]))
	local rhoR, uxR, uyR, uzR, BxR, ByR, BzR, pR = self:stateToPrims(unpack(self.qs[i]))

	local rho = .5*(rhoL + rhoR)
	local ux = .5*(uxL + uxR)
	local uy = .5*(uyL + uyR)
	local uz = .5*(uzL + uzR)
	local Bx = .5*(BxL + BxR)
	local By = .5*(ByL + ByR)
	local Bz = .5*(BzL + BzR)
	local p = .5*(pL + pR)

	local mu = self.mu
	local BSq = Bx*Bx + By*By + Bz*Bz
	local vaxSq = Bx*Bx / (mu*rho)
	local vax = sqrt(vaxSq)
	local CsSq = gamma*p / rho
	local Cs= sqrt(CsSq)
	local vaSq = BSq / (mu*rho)
	local va = sqrt(vaSq)
	local cStarSq = .5*(vaSq + CsSq)
	local discr = sqrt(cStarSq*cStarSq - vaxSq*CsSq)
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
	local epsilon = 1e-20

	local afSq = (vfSq - vaxSq) / (vfSq - vsSq)
	local af
	if afSq > epsilon then
		af = sqrt(afSq)
	else
		af, afSq = 1, 1
	end
	
	local asSq = (vfSq - CsSq) / (vfSq - vsSq)
	local as
	if asSq > epsilon then
		as = sqrt(asSq)
	else
		as, asSq = 1, 1
	end

	local BTSq = By*By + Bz*Bz
	local betay = 0
	local betaz = 0
	if BTSq > epsilon then
		local BT = sqrt(BTSq)
		betay = By / BT
		betaz = Bz / BT
	end
		
	local RFast = vf / sqrt(afSq*(vfSq + CsSq) + asSq*(vfSq + vaxSq))
	local RSlow = vfSq / sqrt(afSq*CsSq*(vfSq + CsSq) + asSq*vfSq*(vsSq + CsSq))

	local sqrt1_2 = sqrt(.5)

	local R_A = {	--right eigenvectors
	--fast -
		{
			af*tau*RFast,
			af*vf*RFast,
			-as*betay*vax*sgnBx*RFast,
			-as*betaz*vax*sgnBx*RFast,
			0,
			-as*betay*vf*sqrt(mu*rho)*RFast,
			-as*betaz*vf*sqrt(mu*rho)*RFast,
			-af*gamma*p*RFast
		},
	--alfven -
		{
			0,
			0,
			-betaz*vf*sqrt1_2,
			betay*vf*sqrt1_2,
			0,
			-sgnBx*sqrt(mu*rho)*betaz*vf*sqrt1_2,
			sgnBx*sqrt(mu*rho)*betay*vf*sqrt1_2,
			0
		},
	--slow -
		{
			as*tau,
			as*vs,
			af*betay*Cs*sgnBx,
			af*betaz*Cs*sgnBx,
			0,
			af*betay*CsSq / vf*sqrt(mu*rho),
			af*betaz*CsSq / vf*sqrt(mu*rho),
			-as*gamma*p
		},
	--entropy
		{tau, 0, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 1, 0, 0, 0},
	--zero
	--slow +
		{
			-as*tau,
			as*vs,
			af*betay*Cs*sgnBx,
			af*betaz*Cs*sgnBx,
			0,
			-af*betay*CsSq / vf*sqrt(mu*rho),
			-af*betaz*CsSq / vf*sqrt(mu*rho),
			as*gamma*p
		},
	--alfven +
		{
			0,
			0,
			betaz*vf*sqrt1_2,
			-betay*vf*sqrt1_2,
			0,
			-sgnBx*sqrt(mu*rho)*betaz*vf*sqrt1_2,
			sgnBx*sqrt(mu*rho)*betay*vf*sqrt1_2,
			0
		},
	--fast +
		{
			-af*tau*RFast,
			af*vf*RFast,
			-as*betay*vax*sgnBx*RFast,
			-as*betaz*vax*sgnBx*RFast,
			0,
			as*betay*vf*sqrt(mu*rho)*RFast,
			as*betaz*vf*sqrt(mu*rho)*RFast,
			af*gamma*p*RFast
		}
	}

	--specify in rows (so it will look as-is)
	--so L_ij = L_A[i][j]
	local L_A = {
	--fast -
		{
			0,
			 af*vf*RFast / vfSq,
			 -as*betay*vax*sgnBx*RFast / vfSq,
			 -as*betaz*vax*sgnBx*RFast / vfSq,
			 0,
			 -as*betay*vf / sqrt(mu*rho)*RFast / vfSq,
			 -as*betaz*vf / sqrt(mu*rho)*RFast / vfSq,
			 -af*tau*RFast / vfSq
		},
	--alfven -
		{
			0,
			0,
			-betaz*sqrt1_2 / vf,
			betay*sqrt1_2 / vf,
			0,
			-sgnBx*betaz / sqrt(mu*rho)*sqrt1_2 / vf,
			sgnBx*betay / sqrt(mu*rho)*sqrt1_2 / vf,
			0
		},
	--slow -
		{
			0,
			as*vs*RSlow / vfSq,
			af*betay*Cs*sgnBx*RSlow / vfSq,
			af*betaz*Cs*sgnBx*RSlow / vfSq,
			0,
			af*betay*CsSq / (sqrt(mu*rho)*vf)*RSlow / vfSq,
			af*betaz*CsSq / (sqrt(mu*rho)*vf)*RSlow / vfSq,
			-as*gamma*p*RSlow / vfSq
		},
	--entropy
		{rho, 0, 0, 0, 0, 0, 0, 1 / (gamma*p)},
	--zero
		{0, 0, 0, 0, 1, 0, 0, 0},
	--slow +
		{
			0,
			as*vs*RSlow / vfSq,
			af*betay*Cs*sgnBx*RSlow / vfSq,
			af*betaz*Cs*sgnBx*RSlow / vfSq,
			0,
			-af*betay*CsSq / (sqrt(mu*rho)*vf)*RSlow / vfSq,
			-af*betaz*CsSq / (sqrt(mu*rho)*vf)*RSlow / vfSq,
			as*gamma*p*RSlow / vfSq
		},
	--alfven +
		{
			0,
			0, 
			betaz*sqrt1_2 / vf,
			-betay*sqrt1_2 / vf,
			0, 
			-sgnBx*betaz / sqrt(mu*rho)*sqrt1_2 / vf,
			sgnBx*betay / sqrt(mu*rho)*sqrt1_2 / vf,
			0
		},
	--fast +
		{
			0,
			 af*vf*RFast / vfSq,
			 -as*betay*vax*sgnBx*RFast / vfSq,
			 -as*betaz*vax*sgnBx*RFast / vfSq,
			 0,
			 as*betay*vf / sqrt(mu*rho)*RFast / vfSq,
			 as*betaz*vf / sqrt(mu*rho)*RFast / vfSq,
			 af*tau*RFast / vfSq
		}
	}

	local rhoSq = rho*rho
	local tauSq = tau*tau
	local uSq = ux*ux + uy*uy + uz*uz

	--specified by row, so dU_dW[i][j] == (dU/dW)_ij
	local dU_dW = {
		{-rhoSq,	0,	0,	0,	0,	0,	0,	0},
		{-ux*rhoSq,	rho,	0,	0,	0,	0,	0,	0},
		{-uy*rhoSq,	0,	rho,	0,	0,	0,	0,	0},
		{-uz*rhoSq,	0,	0,	rho,	0,	0,	0,	0},
		{0,	0,	0,	0,	1,	0,	0,	0},
		{0,	0,	0,	0,	0,	1,	0,	0},
		{0,	0,	0,	0,	0,	0,	1,	0},
		{.5*BSq/mu,	ux,	uy,	uz,	Bx*tau/mu,	By*tau/mu,	Bz*tau/mu,	1/gammaMinusOne}
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
		{gammaMinusOne*tau*(uSq + .5*BSq*tau/mu),	-gammaMinusOne*tau*ux,	-gammaMinusOne*tau*uy,	-gammaMinusOne*tau*uz,	-gammaMinusOne*Bx*tau/mu,	-gammaMinusOne*By*tau/mu,	-gammaMinusOne*Bz*tau/mu,	gammaMinusOne}
	}
	--now transform these to the left and right eigenvectors of the flux ...
	--with transformations: R_U = dU/dW*R_A and L_U = L_A*dW/dU
	--don't forget indexing is A_ij == A[i][j] except R_A is transposed
	for j=1,8 do
		for k=1,8 do
			local sum = 0
			for m=1,8 do
				sum = sum + dU_dW[j][m]*R_A[k][m]
			end
			self.eigenvectors[i][j][k] = sum
		end
	end
	for j=1,8 do
		for k=1,8 do
			local sum = 0
			for m=1,8 do
				sum = sum + L_A[j][m]*dW_dU[m][k]
			end
			self.eigenvectorsInverse[i][j][k] = sum
		end
	end
	do
		local totalError = 0
		for j=1,8 do
			for k=1,8 do
				local sum = 0
				for m=1,8 do
					sum = sum + dU_dW[j][m] * dW_dU[m][k]
				end
				totalError = totalError + abs((j==k and 1 or 0) - sum)
			end
		end
		self.primOrthoError[i] = totalError
	end
	do
		local totalError = 0
		for j=1,8 do
			for k=1,8 do
				local sum = 0
				for m=1,8 do
					sum = sum + L_A[j][m] * R_A[k][m]
				end
				totalError = totalError + abs((j==k and 1 or 0) - sum)
			end
		end
		self.consOrthoError[i] = totalError
	end	
end

return MHDSimulation