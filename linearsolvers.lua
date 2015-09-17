local function conjgrad(args)
	return require 'LinearSolvers.ConjugateGradient'{
		A = args.A,
		b = args.b,
		x0 = args.x0,
		maxiter = args.maxiter or 100,
		epsilon = args.epsilon or 1e-10,
		dot = args.b.dot,
		clone = args.b.clone,
	}
end

local function conjres(args)
	return require 'LinearSolvers.ConjugateResidual'{
		A = args.A,
		b = args.b,
		x0 = args.x0,
		maxiter = args.maxiter or 100,
		epsilon = args.epsilon or 1e-10,
		dot = args.b.dot,
		clone = args.b.clone,
	}
end

local function jacobi(args)
	return require 'LinearSolvers.Jacobi'{
		A = args.A,
		ADiag = args.ADiag,
		b = args.b,
		x0 = args.x0,
		maxiter = args.maxiter or 100,
		epsilon = args.epsilon or 1e-10,
		scale = b.perElementMultiply,
		invScale = b.perElementDivide,
	}
end

return {
	conjgrad = conjgrad,
	conjres = conjres,
	jacobi = jacobi,
}

