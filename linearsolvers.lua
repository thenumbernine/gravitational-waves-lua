local function conjgrad(args)
	local A = args.A
	local b = args.b
	local x = args.x0
	local maxiter = args.maxiter or 100
	local epsilon = args.epsilon or 1e-10
	
	local r = b - A(x)
	local r2 = r:dot(r)
	if r2 < epsilon then return x end
	local p = r:clone()
	for iter=1,maxiter do
		local Ap = A(p)
		local alpha = r2 / p:dot(Ap)
		x = x + alpha * p
		local nr = r - alpha * Ap
		local nr2 = nr:dot(nr)
		local beta = nr2 / r2
		if nr2 < epsilon then break end
		r = nr
		r2 = nr2
		p = r + beta * p
	end
	return x
end

local function conjres(args)
	local A = args.A
	local b = args.b
	local x = args.x0
	local maxiter = args.maxiter or 100
	local epsilon = args.epsilon or 1e-10

	local r = b - A(x)
	local r2 = r:dot(r)
	if r2 <= epsilon then return x end
	
	local Ar = A(r)
	local rAr = r:dot(Ar)
	local p = r:clone()
	local Ap = A(p)
	for iter=1,maxiter do
		local alpha = rAr / Ap:dot(Ap)
		x = x + alpha * p
		local nr = r - alpha * Ap
		local Anr = A(nr)
		local nrAr = nr:dot(Anr)
		local beta = nrAr / rAr
		local nr2 = nr:dot(nr)
		if nr2 < epsilon then break end
		r = nr
		rAr = nrAr
		Ar = Anr
		p = r + beta * p
		Ap = Ar + beta * Ap
	end
	return x
end

local function jacobi(args)
	local A = args.A
	local A_diag = args.A_diag
	local b = args.b
	local x = args.x0
	local maxiter = args.maxiter or 100
	local epsilon = args.epsilon or 1e-10
	
	local x = b:clone()
	for iter=1,maxiter do
		local nx = b - A(x)
		nx = nx + x:perElementMultiply(A_diag)	-- remove diagonal
		nx = nx:perElementDivide(A_diag)		-- divide by diagonal
		local r = b - A(nx)
		local r2 = r:dot(r)
		if r2 < epsilon then break end
		x = nx
	end
	return x
end

return {
	conjgrad = conjgrad,
	conjres = conjres,
	jacobi = jacobi,
}

