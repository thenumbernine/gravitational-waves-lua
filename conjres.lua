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

return conjres

