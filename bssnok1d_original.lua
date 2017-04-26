local class = require 'ext.class'
local Equation = require 'equation'
local BSSNOK1DOriginal = class(Equation)

BSSNOK1DOriginal.name = 'BSSNOK 1D Original'
BSSNOK1DOriginal.numStates = 6

function BSSNOK1DOriginal:init(args, ...)
	
	local symmath = require 'symmath'
	symmath.tostring = require 'symmath.tostring.SingleLine'
	local function makesym(field)
		return symmath.clone(assert(args[field], "expected to find field "..field))() 
	end
	
	-- parameters that are variables of symbolic functions
	local x = assert(args.x)

	-- parameters that are symbolic functions -- based on coordinates 
	local exprs = table{'alpha', 'gamma_xx', 'K_xx'}:map(function(name)
		return makesym(name), name
	end)

	-- Alcubierre uses tildes for both kinds of rescaling
	-- Baumgarte & Shapiro uses bars vs tildes
	-- gamma_ij = exp(4 phi) gammaBar_ij
	-- det gammaBar_ij = 1
	-- det gamma_ij = det exp(4 phi) gammaBar_ij = exp(12 phi) det gammaBar_ij = exp(12 phi)
	-- phi = 1/12 log det gamma_ij
	exprs.phi = (symmath.log(exprs.gamma_xx)/12)()
	exprs.K = (exprs.K_xx / exprs.gamma_xx)()
	exprs.gammaTilde_xx = (exprs.gamma_xx / exprs.phi^4)()
	exprs.gammaTildeUxx = (1 / exprs.gammaTilde_xx)()
	exprs.ATilde_xx = (symmath.exp(-4 * exprs.phi) * (exprs.K_xx - exprs.K/3 * exprs.gamma_xx))()
	
	-- GammaTilde^i = exp(4 phi) Gamma^i + 2 gammaTilde^ij phi_,j
	--	= -gammaTilde^ij_,j
	exprs.GammaTildeUx = (-exprs.gammaTildeUxx:diff(x))()

	-- convert from symbolic functions to Lua functions
	self.calc = exprs:map(function(expr, name)
		local f, code = expr:compile{x}
print()
print(name)
print(code)
print(expr)
		return f, name
	end)
	
	-- parameters that are symbolic functions -- based on f's alpha
	-- NOTICE: assign these *after* mapping exprs to calc
	exprs.f = makesym'f'
	local f_param = assert(args.f_param)
	self.calc.f = exprs.f:compile{f_param}
	
	local dalpha_f = exprs.f:diff(args.f_param)()
	self.calc.dalpha_f = dalpha_f:compile{args.f_param}
end

--phi = -1/(4*n) ln gamma_xx
-- exp(-4n phi) = gamma_xx for n=1
-- volume = sqrt(gamma_xx) = sqrt(exp(-4n phi)) = exp(-2n phi)
do
	local q = function(self,i) return self.qs[i] end
	local alpha = q:_(1)
	local phi = q:_(2)
	local gammaTilde_xx = q:_(3)
	local GammaTildeUx = q:_(4)
	local gamma_xx = math.exp:o(-12 * phi)
	local K = q:_(5)
	local ATilde_xx = q:_(6)
	local volume = alpha * math.sqrt:o(gamma_xx)
	BSSNOK1DOriginal:buildGraphInfos{
		{alpha = alpha},
		{phi = phi},
		{K = K},
		{ATilde_xx = ATilde_xx},
		{gammaTilde_xx = gammaTilde_xx},
		{GammaTildeUx = GammaTildeUx},
		{volume = volume},
	}
end

function BSSNOK1DOriginal:initCell(sim,i)
	local x = sim.xs[i]
	local alpha = self.calc.alpha(x)
assert(math.isfinite(alpha))	
	local phi = self.calc.phi(x)
assert(math.isfinite(phi))	
	local gammaTilde_xx = self.calc.gammaTilde_xx(x)
assert(math.isfinite(gammaTilde_xx))	
	local GammaTildeUx = self.calc.GammaTildeUx(x)
assert(math.isfinite(GammaTildeUx))	
	local K = self.calc.K(x)
assert(math.isfinite(K))	
	local ATilde_xx = self.calc.ATilde_xx(x)
assert(math.isfinite(ATilde_xx))	
	return {alpha, phi, gammaTilde_xx, GammaTildeUx, K, ATilde_xx}
end

return BSSNOK1DOriginal
