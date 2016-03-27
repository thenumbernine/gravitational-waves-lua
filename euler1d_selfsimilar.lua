local class = require 'ext.class'
local Euler1DEqn = require 'euler1d'
local Euler1DSSEqn = class(Euler1DEqn)

-- used by HLL
function Euler1DSSEqn:calcFluxForState(sim, i, q, flux)
	local flux = Euler1DSSEqn.super.calcFluxForState(self, sim, i, q, flux)
	for j=1,3 do
		flux[j] = flux[j] - sim.ixs[i] * q[j]
	end
	return flux
end

-- used by HLL
-- offset self-similar eigenvalues by coordinate
function Euler1DSSEqn:calcInterfaceEigenvalues(sim, i, qL, qR, S)
	local S = Euler1DSSEqn.super.calcInterfaceEigenvalues(self, sim, i, qL, qR, S)
	for j=1,3 do
		S[j] = S[j] - sim.ixs[i]
	end
end

-- used by Roe
function Euler1DSSEqn:calcInterfaceEigenBasis(sim,i,qL,qR)
	Euler1DSSEqn.super.calcInterfaceEigenBasis(self,sim,i,qL,qR)
	local S = sim.eigenvalues[i]
	for j=1,3 do
		S[j] = S[j] - sim.ixs[i]
	end
	local F = sim.fluxMatrix[i]
	for j=1,3 do
		F[j][j] = F[j][j] - sim.ixs[i]
	end
end


local class = require 'ext.class'
local table = require 'ext.table'
local Roe = require 'roe'
local HLL = require 'hll'
local matrix = require 'matrix'
local solveQR = require 'LinearSolvers.solveQR'
local solveLUP = require 'LinearSolvers.solveLUP'
local Euler1DSS = class(Roe)
--local Euler1DSS = class(HLL)

-- solves x for A x = b
solveLinear = solveQR
--solveLinear = solveLUP

function Euler1DSS:init(args)
	args.equation = Euler1DSSEqn()
	Euler1DSS.super.init(self, args)
end

function Euler1DSS:iterate()
	-- instead of integrating, we're going to do a newton descent
	-- don't apply boundary. don't touch boundary cells.
	self:applyBoundary()
	
	local n = 3*(self.gridsize - 4)

	-- push self.qs into U
	local function pushState(qs)
		local U = matrix.zeros(n)
		for i=3,self.gridsize-2 do	-- exclude 2 ghost cells
			for j=1,3 do
				U[j+3*(i-3)] = qs[i][j]
			end
		end
		return U
	end
	
	local function popState(U)
		for i=3,self.gridsize-2 do	-- exclude 2 ghost cells
			for j=1,3 do
				self.qs[i][j] = U[j+3*(i-3)]
			end
		end
		self:applyBoundary()
	end
	
	local function calcG(U)
		-- G = U + dF/dxi
		-- store U in self.qs
		popState(U)
		local dt = self:calcDT()	-- calc eigenbasis for flux to use
		local dt = 0				-- dt used for slope limiter
		-- calculate dq_dt = -df/dx based on self.qs
		local dF_dx = -pushState(self:calcFlux(dt))
		return U + dF_dx
	end
	
	local U = pushState(self.qs)
	
	local dG_dU_T = matrix.zeros(n)
	local epsilon = 1e-2
	for j=1,n do	-- j varies across d/dU
		local Uplus = matrix(U)
		Uplus[j] = Uplus[j] + epsilon
		local Gplus = calcG(Uplus)
		
		local Uminus = matrix(U)
		Uminus[j] = Uminus[j] - epsilon
		local Gminus = calcG(Uminus)
	
		dG_dU_T[j] = (Gplus - Gminus) / (2 * epsilon)
	end
	local dG_dU = dG_dU_T:transpose()

	local G = calcG(U)
	local dU = matrix(solveLinear(dG_dU, G))
	
	local alpha = .0001	-- alpha is the line trace coefficient
	U = U - alpha * dU
	-- in absense of jacobian (J=I) and alpha=1 we should get ...
	-- U = U - dU = U - (U + dF/dx) = -dF/dx	... of course it'll reach a steady state of zero
	-- with jacobian ...

	popState(U)
end

return Euler1DSS 
