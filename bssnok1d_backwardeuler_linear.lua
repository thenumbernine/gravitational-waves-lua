--[[

diff_t alpha = -alpha^2 f K
diff_t gammaTilde_ij = -2 alpha ATilde_ij
diff_t phi = -1/6 alpha K

diff_t ATilde_ij =
	exp(-4 phi) (
		-diff_i diff_j alpha
		+ 1/3 gammaTilde_ij gammaTilde^kl diff_k diff_l alpha
		
		+ GammaTilde^k_ji diff_k alpha
		- 1/3 gammaTilde_ij GammaTilde^k diff_alpha
		
		+ 2 delta^k_i diff_j phi diff_k alpha
		- 2/3 gamma_ij gamma^mn delta^k_m diff_n phi diff_k alpha
		
		+ 2 delta^k_j diff_i phi diff_k alpha
		- 2/3 gamma_ij gamma^mn delta^k_n diff_m phi diff_k alpha
		
		+ alpha R_ij
		- 1/3 alpha gamma_ij R

		- 8 pi alpha S_ij
		+ 8/3 pi alpha gamma_ij gamma^kl S_kl
	)
	+ alpha K ATilde_ij
	- 2 alpha ATilde_ik gammaTilde^kl ATilde_lj

diff_t K =
	exp(-4 phi) (
		-gammaTilde^ij diff_i diff_j alpha
		+ GammaTilde^k diff_k alpha
		- 2 gammaTilde^ij diff_j phi diff_i alpha
	)
	+ alpha ATilde_ij gammaTilde^im gammaTilde^jn ATilde_mn
	+ 1/3 alpha K^2
	+ 4 pi alpha (rho + S)

diff_t GammaTilde^i =
	-2 gammaTilde^im ATilde^mn gammaTilde^nj diff_j alpha
	+ 2 alpha GammaTilde^i_jk gammaTilde^jm ATilde^mn gammaTilde^nk
	+ 12 alpha gammaTilde^im ATilde_mn gammaTilde^nj diff_j phi
	- 4/3 alpha gammaTilde^ij diff_j K
	- 16 pi alpha jTilde^i

--]]

local class = require 'ext.class'
local table = require 'ext.table'
local SolverFV = require 'solverfv'
local BSSNOK1D = require 'bssnok1d'

local BSSNOK1DBackwardEulerLinear = class(SolverFV)
BSSNOK1DBackwardEulerLinear.name = 'BSSNOK 1D Backward Euler Linear'

function BSSNOK1DBackwardEulerLinear:init(args)
	args = table(args)
	args.equation = BSSNOK1D(args)
	self.linearSolver = args.linearSolver or require 'linearsolvers'.conjres
	self.fluxMatrix = {}
	self.eigenvectors = {}
	self.eigenvectorsInverse = {}
	self.eigenbasisErrors = {}
	self.fluxMatrixErrors = {}
	BSSNOK1DBackwardEulerLinear.super.init(self, args)
end

function BSSNOK1DBackwardEulerLinear:calcDT(getLeft, getRight)
	for i=2,self.gridsize do
		local qL = getLeft and getLeft(self,i) or self.qs[i-1]
		local qR = getRight and getRight(self,i) or self.qs[i]
		self.equation:calcInterfaceEigenBasis(self,i,qL,qR)
	end
	BSSNOK1DBackwardEulerLinear.super.calcDT(self, getLeft, getRight)
end

function BSSNOK1DBackwardEulerLinear:step(dt, getLeft, getRight)
	local function calc_dq_dt(qs)
		local dq_dt = self:newState() 
		-- [=[
		for i=1,self.gridsize do
			local alpha, phi, A_x, Phi_x, K, ATilde_xx = table.unpack(qs[i])
			local f = self.equation.calc.f(alpha)
			local dalpha_f = self.equation.calc.dalpha_f(alpha)
			
			--alpha,t = -alpha^2 f K
			dq_dt[i][1] = dq_dt[i][1] - alpha * alpha * f * K
			--phi,t = -alpha K / 6
			dq_dt[i][2] = dq_dt[i][2] - alpha * K / 6
			
			local dx = self.ixs[i+1] - self.ixs[i]
			if i>1 then
				--A_x,t = -alpha K f,alpha partial_x alpha - alpha f partial_x K
				dq_dt[i][3] = dq_dt[i][3] + alpha * K * dalpha_f * qs[i-1][1] / (2 * dx)
				dq_dt[i][3] = dq_dt[i][3] + alpha * f * qs[i-1][5] / (2 * dx)
				--Phi_x,t = -alpha/6 partial_x K
				dq_dt[i][4] = dq_dt[i][4] + alpha / 6 * qs[i-1][5] / (2 * dx)
				--K,t = -alpha exp(-4 phi) partial_x A_x
				dq_dt[i][5] = dq_dt[i][5] + alpha * math.exp(-4 * phi) * qs[i-1][3] / (2 * dx)
				--ATilde_xx,t = -alpha exp(-4 phi) partial_x A_x - 2 alpha exp(-4 phi) partial_x Phi_x
				dq_dt[i][6] = dq_dt[i][6] + alpha * math.exp(-4 * phi) * qs[i-1][3] / (2 * dx)
				dq_dt[i][6] = dq_dt[i][6] + 2 * alpha * math.exp(-4 * phi) * qs[i-1][4] / (2 * dx)
			end
			if i<self.gridsize then
				--A_x,t = -alpha K f,alpha partial_x alpha - alpha f partial_x K
				dq_dt[i][3] = dq_dt[i][3] - alpha * K * dalpha_f * qs[i+1][1] / (2 * dx)
				dq_dt[i][3] = dq_dt[i][3] - alpha * f * qs[i+1][5] / (2 * dx)
				--Phi_x,t = -alpha/6 partial_x K
				dq_dt[i][4] = dq_dt[i][4] - alpha / 6 * qs[i+1][5] / (2 * dx)
				--K,t = -alpha exp(-4 phi) partial_x A_x
				dq_dt[i][5] = dq_dt[i][5] - alpha * math.exp(-4 * phi) * qs[i+1][3] / (2 * dx)
				--ATilde_xx,t = -alpha exp(-4 phi) partial_x A_x - 2 alpha exp(-4 phi) partial_x Phi_x
				dq_dt[i][6] = dq_dt[i][6] - alpha * math.exp(-4 * phi) * qs[i+1][3] / (2 * dx)
				dq_dt[i][6] = dq_dt[i][6] - 2 * alpha * math.exp(-4 * phi) * qs[i+1][4] / (2 * dx)
			end
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

return BSSNOK1DBackwardEulerLinear
