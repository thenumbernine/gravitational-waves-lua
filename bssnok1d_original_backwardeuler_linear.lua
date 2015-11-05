--[[
Alcubierre p.84 / Baumgarte & Shapiro p.394 or so

d/dt = partial_t - Lie_beta

d/dt alpha = -alpha^2 f K 		<- Bona-Masso family slicing
alpha,t = -alpha^2 f K + beta^k (-alpha^2 f K),k
alpha,t =
	- 2 alpha beta^k f K alpha,k
	+ beta^k alpha,k

d/dt phi = -1/6 alpha K
phi,t =
	- 1/6 alpha K
	+ beta^k phi,k

d/dt gammaTilde_ij = -2 alpha ATilde_ij
gammaTilde_ij,t = 
	- 2 alpha ATilde_ij
	+ beta^k gammaTilde_ij,k

d/dt K = -D_i D^i alpha + alpha (ATilde_ij ATilde^ij + 1/3 K^2) + 4 pi alpha (rho + S)
d/dt K = -D_i (alpha,i) + alpha (ATilde_ij ATilde^ij + 1/3 K^2) + 4 pi alpha (rho + S)

d/dt ATilde_ij = exp(-4 phi) (-D_i D_j alpha + alpha R_ij + 4 pi alpha (gamma_ij (S - rho) - 2 S_ij))^TF	(TF=trace-free)
d/dt GammaTilde^i = - 2 ATilde^ij alpha_,j + 2 alpha (GammaTilde^i_jk ATilde^jk + 6 ATilde^ij phi_,j - 2/3 gammaTilde^ij K,j - 8 pi jTilde^i) + gammaTilde^jk beta^i_,kj + 1/3 gammaTilde^ij beta^k_,kj
--]]
local class = require 'ext.class'
local table = require 'ext.table'
local SolverFV = require 'solverfv'
local BSSNOK1D = require 'bssnok1d'

local BSSNOK1DOriginalBackwardEuler = class(SolverFV)
BSSNOK1DOriginalBackwardEuler.name = 'BSSNOK Original 1D Backward Euler Linear'

function BSSNOK1DOriginalBackwardEuler:init(args)
	args = table(args)
	args.equation = BSSNOK1D(args)
	self.linearSolver = args.linearSolver or require 'linearsolvers'.conjres
	self.fluxMatrix = {}
	self.eigenvectors = {}
	self.eigenvectorsInverse = {}
	self.eigenbasisErrors = {}
	self.fluxMatrixErrors = {}
	BSSNOK1DOriginalBackwardEuler.super.init(self, args)
end

function BSSNOK1DOriginalBackwardEuler:calcDT(getLeft, getRight)
	for i=2,self.gridsize do
		local qL = getLeft and getLeft(self,i) or self.qs[i-1]
		local qR = getRight and getRight(self,i) or self.qs[i]
		self.equation:calcInterfaceEigenBasis(self,i,qL,qR)
	end
	BSSNOK1DOriginalBackwardEuler.super.calcDT(self, getLeft, getRight)
end

function BSSNOK1DOriginalBackwardEuler:step(dt, getLeft, getRight)
	local function calc_dq_dt(qs)
		local dq_dt = self:newState() 
		-- [=[
		for i=1,self.gridsize do
			local alpha, phi, A_x, Phi_x, K, ATilde_xx = table.unpack(qs[i])
			local f = self.equation.calc.f(alpha)
			local dalpha_f = self.equation.calc.dalpha_f(alpha)
		
			-- TODO compute derivative 
		
		end
		--]=]
		return dq_dt
	end

	--[=[ forward euler
	self.qs = self.qs + dt * calc_dq_dt(self.qs)
	--]=]
	-- [=[ backward euler
	local qs = self.qs
	self.qs = self.linearSolver{
		b = qs:clone(),
		x0 = qs:clone(),
		A = function(qs)
			qs = qs - dt * calc_dq_dt(qs)
			--self.boundaryMethod(qs)
			return qs
		end,
	}
	--]=]
end

return BSSNOK1DOriginalBackwardEuler
