local table = require 'ext.table'

local function prepArgs(args)
	args = table(args)
	assert(args.A)
	assert(args.b)
	-- args.x is optional, defaults to b
	-- args.ADiag is only used by jacobi
	if not args.maxiter then args.maxiter = 10 * #args.b end
	if not args.restart then args.restart = #args.b end	-- used by gmres
	if not args.epsilon then args.epsilon = 1e-10 end
	if not args.dot then args.dot = assert(args.b.dot) end
	if not args.clone then args.clone = assert(args.b.clone) end
	if not args.scale then args.scale = args.b.emul end
	if not args.invScale then args.invScale = args.b.ediv end
	return args
end

local solvers = {}
for _,name in ipairs{'jacobi', 'conjgrad', 'conjres', 'bicgstab', 'gmres'} do
	local solver = require('solver.'..name)
	solvers[name] = function(args)
		return solver(prepArgs(args)) 
	end
end
return solvers
