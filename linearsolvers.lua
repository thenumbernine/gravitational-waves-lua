local table = require 'ext.table'
local Jacobi = require 'LinearSolvers.Jacobi'
local ConjugateGradient = require 'LinearSolvers.ConjugateGradient'
local ConjugateResidual = require 'LinearSolvers.ConjugateResidual'
local BiconjugateGradientStabilized = require 'LinearSolvers.BiconjugateGradientStabilized'
local GeneralizedMinimalResidual = require 'LinearSolvers.GeneralizedMinimalResidual'

local function prepArgs(args)
	args = table(args)
	assert(args.A)
	assert(args.b)
	-- args.x0 is optional, defaults to b
	-- args.ADiag is only used by jacobi
	if not args.maxiter then args.maxiter = 100 end
	if not args.epsilon then args.epsilon = 1e-10 end
	if not args.dot then args.dot = args.b.dot end
	if not args.clone then args.clone = args.b.clone end
	if not args.scale then args.scale = args.b.perElementMultiply end
	if not args.invScale then args.invScale = args.b.perElementDivide end
	return args
end

return {
	jacobi = function(args)
		return Jacobi(prepArgs(args))
	end,
	conjgrad = function(args)
		return ConjugateGradient(prepArgs(args))
	end,
	conjres = function(args)
		return ConjugateResidual(prepArgs(args))
	end,
	bicgstab = function(args)
		return BiconjugateGradientStabilized(prepArgs(args))
	end,
	gmres = function(args)
		return GeneralizedMinimalResidual(prepArgs(args))
	end,
}

