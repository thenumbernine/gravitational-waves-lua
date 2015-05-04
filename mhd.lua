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
	--[[ brio-wu
	local Bx = .75
	local By = x < 0 and 1 or -1
	local Bz = 0
	local p = x < 0 and 1 or .1
	--]]
	-- [[ sod
	local Bx = 0
	local By = 0
	local Bz = 0
	local p = 1
	--]]
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
	return rho, ux, uy, uz, Bx, By, Bz, p, ETotal
end

function MHDSimulation:calcInterfaceEigenBasis(i)
	local gamma = self.gamma
	local gammaMinusOne = gamma - 1
	local gammaMinusTwo = gamma - 2
	local mu = self.mu

	local rhoL, uxL, uyL, uzL, BxL, ByL, BzL, pL, ETotalL = self:stateToPrims(unpack(self.qs[i-1]))
	
	local rhoR, uxR, uyR, uzR, BxR, ByR, BzR, pR, ETotalR = self:stateToPrims(unpack(self.qs[i]))

	local sqrtRhoL = sqrt(rhoL)
	local sqrtRhoR = sqrt(rhoR)
	local sumSqrtRho = sqrtRhoL + sqrtRhoR
	local invDenom = 1 / sumSqrtRho

	local rho = sqrtRhoL*sqrtRhoR
	local ux = (sqrtRhoL*uxL + sqrtRhoR*uxR)*invDenom
	local uy = (sqrtRhoL*uyL + sqrtRhoR*uyR)*invDenom
	local uz = (sqrtRhoL*uzL + sqrtRhoR*uzR)*invDenom
	local Bx = (sqrtRhoL*BxL + sqrtRhoR*BxR)*invDenom
	local By = (sqrtRhoL*ByL + sqrtRhoR*ByR)*invDenom
	local Bz = (sqrtRhoL*BzL + sqrtRhoR*BzR)*invDenom

	local x = .5*((ByL - ByR)^2 + (BzL - BzR)^2)*invDenom^2
	local y = .5*sumSqrtRho/rho;
	local pbL = .5*(Bx*Bx + ByL*ByL + BzL*BzL)	
	local pbR = .5*(Bx*Bx + ByR*ByR + BzR*BzR)

	local h = ((ETotalL + pL + pbL)/sqrtRhoL + (ETotalR + pR + pbR)/sqrtRhoR)*invDenom

	local invRho = 1/rho
	local uSq = ux*ux + uy*uy + uz*uz
	local BTSq = By*By + Bz*Bz
 	local bt_starsq = (gammaMinusOne - gammaMinusTwo*y)*BTSq
	local vaxSq = Bx*Bx*invRho
	local hp = h - (vaxSq + BTSq*invRho)
	local epsilon = 1e-20
	local twid_asq = max((gammaMinusOne*(hp-0.5*uSq)-gammaMinusTwo*x), epsilon)

	local ct2 = bt_starsq*invRho
	local tsum = vaxSq + ct2 + twid_asq
	local tdif = vaxSq + ct2 - twid_asq
	local cf2_cs2 = sqrt(tdif*tdif + 4*twid_asq*ct2)

	local cfsq = .5*(tsum + cf2_cs2)
	local cf = sqrt(cfsq)

	local cssq = twid_asq*vaxSq/cfsq
	local cs = sqrt(cssq)

	local BT = sqrt(BTSq)
	local BTStar = sqrt(bt_starsq)
	local beta_y, beta_z
	if BT == 0 then
		beta_y = 1.0
		beta_z = 0.0
	else
		beta_y = By/BT
		beta_z = Bz/BT
	end
	local betaStar_y = beta_y/sqrt(gammaMinusOne - gammaMinusTwo*y)
	local betaStar_z = beta_z/sqrt(gammaMinusOne - gammaMinusTwo*y)
	local betaStarSq = betaStar_y*betaStar_y + betaStar_z*betaStar_z
	local vbet = uy*betaStar_y + uz*betaStar_z

	local alpha_f, alpha_s
	if cfsq - cssq == 0 then
		alpha_f = 1
		alpha_s = 0
	elseif twid_asq - cssq <= 0 then
		alpha_f = 0
		alpha_s = 1
	elseif cfsq - twid_asq <= 0 then
		alpha_f = 1
		alpha_s = 0
	else
		alpha_f = sqrt((twid_asq - cssq)/(cfsq - cssq))
		alpha_s = sqrt((cfsq - twid_asq)/(cfsq - cssq))
	end

	local sqrtd = sqrt(rho)
	local isqrtd = 1.0/sqrtd
	local s = Bx < 0 and -1 or 1
	local twid_a = sqrt(twid_asq);
	local qf = cf*alpha_f*s
	local qs = cs*alpha_s*s
	local af_prime = twid_a*alpha_f*isqrtd
	local as_prime = twid_a*alpha_s*isqrtd
	local afpbb = af_prime*BTStar*betaStarSq
	local aspbb = as_prime*BTStar*betaStarSq

-- Compute eigenvalues (eq. B17) 

	local vax = sqrt(vaxSq)
	local eigenvalues = self.eigenvalues[i]
	eigenvalues[1] = ux - cf
	eigenvalues[2] = ux - vax
	eigenvalues[3] = ux - cs
	eigenvalues[4] = ux
	eigenvalues[5] = ux
	eigenvalues[6] = ux + cs
	eigenvalues[7] = ux + vax
	eigenvalues[8] = ux + cf

-- Right-eigenvectors, stored as COLUMNS (eq. B21) 
-- Note statements are grouped in ROWS for optimization, even though rem[*][n]
-- is the nth right eigenvector

	local right_eigenmatrix = self.eigenvectors[i]
	right_eigenmatrix[1][1] = alpha_f
	right_eigenmatrix[1][2] = 0
	right_eigenmatrix[1][3] = alpha_s
	right_eigenmatrix[1][4] = 1
	right_eigenmatrix[1][5] = 0
	right_eigenmatrix[1][6] = alpha_s
	right_eigenmatrix[1][7] = 0
	right_eigenmatrix[1][8] = alpha_f

	right_eigenmatrix[2][1] = alpha_f*eigenvalues[1]
	right_eigenmatrix[2][2] = 0
	right_eigenmatrix[2][3] = alpha_s*eigenvalues[3]
	right_eigenmatrix[2][4] = ux
	right_eigenmatrix[2][5] = 0
	right_eigenmatrix[2][6] = alpha_s*eigenvalues[6]
	right_eigenmatrix[2][7] = 0
	right_eigenmatrix[2][8] = alpha_f*eigenvalues[8]

	local qa = alpha_f*uy
	local qb = alpha_s*uy
	local qc = qs*betaStar_y
	local qd = qf*betaStar_y
	right_eigenmatrix[3][1] = qa + qc
	right_eigenmatrix[3][2] = -beta_z
	right_eigenmatrix[3][3] = qb - qd
	right_eigenmatrix[3][4] = uy
	right_eigenmatrix[3][5] = 0
	right_eigenmatrix[3][6] = qb + qd
	right_eigenmatrix[3][7] = beta_z
	right_eigenmatrix[3][8] = qa - qc

	local qa = alpha_f*uz
	local qb = alpha_s*uz
	local qc = qs*betaStar_z
	local qd = qf*betaStar_z
	right_eigenmatrix[4][1] = qa + qc
	right_eigenmatrix[4][2] = beta_y
	right_eigenmatrix[4][3] = qb - qd
	right_eigenmatrix[4][4] = uz
	right_eigenmatrix[4][5] = 0
	right_eigenmatrix[4][6] = qb + qd
	right_eigenmatrix[4][7] = -beta_y
	right_eigenmatrix[4][8] = qa - qc

	right_eigenmatrix[5][1] = 0
	right_eigenmatrix[5][2] = 0
	right_eigenmatrix[5][3] = 0
	right_eigenmatrix[5][4] = 0
	right_eigenmatrix[5][5] = 1
	right_eigenmatrix[5][6] = 0
	right_eigenmatrix[5][7] = 0
	right_eigenmatrix[5][8] = 0

	right_eigenmatrix[6][1] = alpha_f*(hp - ux*cf) + qs*vbet + aspbb
	right_eigenmatrix[6][2] = -(uy*beta_z - uz*beta_y)
	right_eigenmatrix[6][3] = alpha_s*(hp - ux*cs) - qf*vbet - afpbb
	right_eigenmatrix[6][4] = 0.5*uSq + gammaMinusTwo*x/gammaMinusOne
	right_eigenmatrix[6][5] = alpha_s*(hp + ux*cs) + qf*vbet - afpbb
	right_eigenmatrix[6][6] = 0
	right_eigenmatrix[6][7] = -right_eigenmatrix[6][2]
	right_eigenmatrix[6][8] = alpha_f*(hp + ux*cf) - qs*vbet + aspbb

	right_eigenmatrix[7][1] = as_prime*betaStar_y;
	right_eigenmatrix[7][2] = -beta_z*s*isqrtd;
	right_eigenmatrix[7][3] = -af_prime*betaStar_y;
	right_eigenmatrix[7][4] = 0
	right_eigenmatrix[7][5] = 0
	right_eigenmatrix[7][6] = right_eigenmatrix[7][3]
	right_eigenmatrix[7][7] = right_eigenmatrix[7][2]
	right_eigenmatrix[7][8] = right_eigenmatrix[7][1]

	right_eigenmatrix[8][1] = as_prime*betaStar_z
	right_eigenmatrix[8][2] = beta_y*s*isqrtd
	right_eigenmatrix[8][3] = -af_prime*betaStar_z
	right_eigenmatrix[8][4] = 0
	right_eigenmatrix[8][5] = 0
	right_eigenmatrix[8][6] = right_eigenmatrix[8][3]
	right_eigenmatrix[8][7] = right_eigenmatrix[8][2]
	right_eigenmatrix[8][8] = right_eigenmatrix[8][1]

	-- Left-eigenvectors, stored as ROWS (eq. B29)

	-- Normalize by 1/2a^{2}: quantities denoted by \hat{f}
	local norm = .5/twid_asq
	local cff = norm*alpha_f*cf
	local css = norm*alpha_s*cs
	qf = qf * norm
	qs = qs * norm
	local af = norm*af_prime*rho
	local as = norm*as_prime*rho
	local afpb = norm*af_prime*BTStar
	local aspb = norm*as_prime*BTStar

	-- Normalize by (gamma-1)/2a^{2}: quantities denoted by \bar{f}
	norm = norm * gammaMinusOne
	local alpha_f = alpha_f * norm
	local alpha_s = alpha_s * norm
	local q2_star = betaStar_y/betaStarSq
	local q3_star = betaStar_z/betaStarSq
	local vqstr = uy*q2_star + uz*q3_star
	norm = norm * 2

	local left_eigenmatrix = self.eigenvectorsInverse[i]
	left_eigenmatrix[1][1] = alpha_f*(uSq-hp) + cff*(cf+ux) - qs*vqstr - aspb
	left_eigenmatrix[1][2] = -alpha_f*ux - cff
	left_eigenmatrix[1][3] = -alpha_f*uy + qs*q2_star
	left_eigenmatrix[1][4] = -alpha_f*uz + qs*q3_star
	left_eigenmatrix[1][5] = 0
	left_eigenmatrix[1][6] = alpha_f
	left_eigenmatrix[1][7] = as*q2_star - alpha_f*By
	left_eigenmatrix[1][8] = as*q3_star - alpha_f*Bz

	left_eigenmatrix[2][1] = .5*(uy*beta_z - uz*beta_y)
	left_eigenmatrix[2][2] = 0
	left_eigenmatrix[2][3] = -0.5*beta_z
	left_eigenmatrix[2][4] = 0.5*beta_y
	left_eigenmatrix[2][5] = 0
	left_eigenmatrix[2][6] = 0
	left_eigenmatrix[2][7] = -0.5*sqrtd*beta_z*s
	left_eigenmatrix[2][8] = 0.5*sqrtd*beta_y*s

	left_eigenmatrix[3][1] = alpha_s*(uSq-hp) + css*(cs+ux) + qf*vqstr + afpb
	left_eigenmatrix[3][2] = -alpha_s*ux - css
	left_eigenmatrix[3][3] = -alpha_s*uy - qf*q2_star
	left_eigenmatrix[3][4] = -alpha_s*uz - qf*q3_star
	left_eigenmatrix[3][5] = 0
	left_eigenmatrix[3][6] = alpha_s
	left_eigenmatrix[3][7] = -af*q2_star - alpha_s*By
	left_eigenmatrix[3][8] = -af*q3_star - alpha_s*Bz

	left_eigenmatrix[4][1] = 1.0 - norm*(0.5*uSq - gammaMinusTwo*x/gammaMinusOne)
	left_eigenmatrix[4][2] = norm*ux
	left_eigenmatrix[4][3] = norm*uy
	left_eigenmatrix[4][4] = norm*uz
	left_eigenmatrix[4][5] = 0
	left_eigenmatrix[4][6] = -norm
	left_eigenmatrix[4][7] = norm*By
	left_eigenmatrix[4][8] = norm*Bz

	left_eigenmatrix[5][1] = 0
	left_eigenmatrix[5][2] = 0
	left_eigenmatrix[5][3] = 0
	left_eigenmatrix[5][4] = 0
	left_eigenmatrix[5][5] = 1
	left_eigenmatrix[5][6] = 0
	left_eigenmatrix[5][7] = 0
	left_eigenmatrix[5][8] = 0

	left_eigenmatrix[6][1] = alpha_s*(uSq-hp) + css*(cs-ux) - qf*vqstr + afpb
	left_eigenmatrix[6][2] = -alpha_s*ux + css
	left_eigenmatrix[6][3] = -alpha_s*uy + qf*q2_star
	left_eigenmatrix[6][4] = -alpha_s*uz + qf*q3_star
	left_eigenmatrix[6][5] = 0
	left_eigenmatrix[6][6] = alpha_s
	left_eigenmatrix[6][7] = left_eigenmatrix[3][7]
	left_eigenmatrix[6][8] = left_eigenmatrix[3][8]

	left_eigenmatrix[7][1] = -left_eigenmatrix[2][1]
	left_eigenmatrix[7][2] = 0
	left_eigenmatrix[7][3] = -left_eigenmatrix[2][3]
	left_eigenmatrix[7][4] = -left_eigenmatrix[2][4]
	left_eigenmatrix[7][5] = 0
	left_eigenmatrix[7][6] = 0
	left_eigenmatrix[7][7] = left_eigenmatrix[2][7]
	left_eigenmatrix[7][8] = left_eigenmatrix[2][8]

	left_eigenmatrix[8][1] = alpha_f*(uSq-hp) + cff*(cf-ux) + qs*vqstr - aspb
	left_eigenmatrix[8][2] = -alpha_f*ux + cff
	left_eigenmatrix[8][3] = -alpha_f*uy - qs*q2_star
	left_eigenmatrix[8][4] = -alpha_f*uz - qs*q3_star
	left_eigenmatrix[8][5] = 0
	left_eigenmatrix[8][6] = alpha_f
	left_eigenmatrix[8][7] = left_eigenmatrix[1][7]
	left_eigenmatrix[8][8] = left_eigenmatrix[1][8]
	
	-- TODO
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


--[[ TODO does athena multiply by dU/dW ?

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
--]]
end

return MHDSimulation

