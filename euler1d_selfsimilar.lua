-- offset lambdas is one way to offset Roe flux, but turns out to be buggy ...
local offsetLambdasByXi = false
-- offset flux by u xi is the other way. works better.
local offsetFluxByUXi = true

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

--[[ I've found this doesn't help
-- used by HLL
-- offset self-similar eigenvalues by coordinate
function Euler1DSSEqn:calcInterfaceEigenvalues(sim, i, qL, qR)
	local lambda = self.eigenvalues[i]
	Euler1DSSEqn.super.calcInterfaceEigenvalues(self, sim, i, qL, qR)
	for j=1,3 do
		lambda[j] = lambda[j] - sim.ixs[i]
	end
end
--]]

-- used by Roe
if offsetLambdasByXi then
	function Euler1DSSEqn:calcInterfaceEigenBasis(sim,i,qL,qR)
		Euler1DSSEqn.super.calcInterfaceEigenBasis(self,sim,i,qL,qR)
		--[[ doesn't help:
		-- one option is offsetting the eigenvalues ... 
		-- mathematically correct, but don't do this
		-- it forms errors in discontinuities
		local S = sim.eigenvalues[i]
		for j=1,3 do
			S[j] = S[j] - sim.ixs[i]
		end
		--]]
		--[[ isn't used
		local F = sim.fluxMatrix[i]
		-- TODO get back the primitives from the interface, and calculate conservative values from them
		for j=1,3 do
			F[j][j] = F[j][j] - sim.ixs[i] * self.iqs[i][j]
		end
		--]]
		error"TODO: get the actual flux, computed at interface, before finite volume calculation, and offset that by -xi * q"
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

if Euler1DSS.super == Roe
and offsetFluxByUXi
then
	function Euler1DSS:calcFluxesAtInterface(...)
		Euler1DSS.super.calcFluxesAtInterface(self, ...)

		-- now offset all fluxes by U * xi
		for i=2,self.gridsize do
			local xi = self.ixs[i]
			-- U is based on the Roe averaged variables
			-- but i don't store them or expose them
			-- so I'd have to recalculate them here ...
			error"TODO"
			--local U = 
			for j=1,sefl.numStates do
				self.fluxes[i] = self.fluxes[i][j] - U[j] * xi
			end
		end
	end
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
	local epsilon = 1e-7
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
	local alpha = .01	-- alpha is the line trace coefficient
	
	U = U - alpha * dU
	-- in absense of jacobian (J=I) and alpha=1 we should get ...
	-- U = U - dU = U - (U + dF/dx) = -dF/dx	... of course it'll reach a steady state of zero
	-- with jacobian ...

	popState(U)
end

return Euler1DSS 
