--[[
going by Alcubierre's book, with some help from the 5-var system in his paper at http://arxiv.org/pdf/gr-qc/9609015v2.pdf 

p.85 non-hyperbolic description:
partial_t gammaTilde_ij = -2 alpha ATilde_ij
partial_t phi = -1/6 alpha K
partial_t ATilde_ij = exp(-4 phi) ( -D_i D_j alpha + alpha R_ij + 4 pi alpha (gamma_ij (S - rho) - 2 S_ij))^TF + alpha (K ATilde_ij - 2 ATilde_ik ATilde^k_j)
partial_t K = -D_i D^i alpha + alpha (ATilde_ij ATilde^ij + 1/3 K^2) + 4 pi alpha (rho + S)
partial_t Gamma^i = -2 ATilde^ij partial_j alpha + 2 alpha (GammaTilde^i_jk ATilde^jk + 6 ATilde^ij partial_j phi - 2/3 gammaTilde^ij partial_j K - 8 pi jTilde^i)
... taking Bona-Masso lapse: partial_t alpha = -alpha^2 f K
... and no shift (hence why all the beta terms aren't there)

extra hyperbolic conditioning variables:
a_i = partial_i ln alpha = partial_i alpha / alpha
Phi_i = partial_i phi
dTilde_ijk = 1/2 partial_k gammaTilde_ij

relations:
K = K^i_i = gamma^ij K_ij
ATilde_ij = exp(-4 phi) (K_ij - 1/3 K gamma_ij) = conformal trace-free extrinsic curvature
GammaTilde^i = gammaTilde^jk connTilde^i_jk

how do we recover phi from the state?
Phi_x = partial_x phi
Phi_x,t = partial_x phi,t
phi,t = integral of x over Phi_x,t
phi,t = integral of x over -alpha/6 partial_x K
...for constant alpha...
phi,t = -alpha K / 6
...and that's what the original BSSNOK says
so do like for alpha and keep track of phi and Phi_x separate

1D relations:
	ND: exp(-phi) = psi <=> phi = -ln psi = -ln (gamma^(1/12) = -1/12 ln gamma
phi = -1/12 ln gamma_xx
gammaTilde_xx = exp(-12 phi) gamma_xx			<- isn't this usually exp(-4 phi)?  courtesy of 1/12 * 3 spatial dimensions = 1/4?
...gammaTilde_xx = exp(-ln gamma_xx) gamma_xx 
...gammaTilde_xx = gamma_xx / gamma_xx
...gammaTilde_xx = 1
dTilde_xxx = 1/2 partial_x gammaTilde_xx
...dTilde_xxx = 1/2 partial_x 1
...dTilde_xxx = 0	<- state variable eliminated
GammaTilde^x
...GammaTilde^x = gammaTilde^xx GammaTilde^x_xx
...GammaTilde^x = 1/2 (partial_k gammaTilde_ij + partial_j gammaTilde_ik - partial_i gammaTilde_jk)
...GammaTilde^x = 0

alpha,t = -alpha^2 f K						<- p.165 or so, which does match the paper, and does include a source term
gammaTilde_xx,t = -2 alpha ATilde_xx		<- p.80, looks like this lines up too
phi,t = -alpha K / 6
A_x,t + alpha partial_x Q = 0
...A_x,t + alpha partial_x (f K) = 0
...A_x,t + alpha K f,alpha partial_x alpha + alpha f partial_x K = 0
Phi_x,t + 1/6 alpha partial_x K = 0
dTilde_xxx,t + alpha partial_x ATilde_xx = 0
K,t + alpha exp(-4 phi) / gammaTilde_xx partial_x A_x = 0
ATilde_xx,t + alpha exp(-4 phi) partial_x LambdaTilde^x_xx = 0
GammaTilde^x,t + 4/3 alpha partial_x (K / gammaTilde_xx) = 0
...GammaTilde^x,t + 4/3 alpha / gammaTilde_xx partial_x K - 4/3 alpha K / gammaTilde_xx^2 partial_x gammaTilde_xx = 0

for LambdaTilde^x_xx = (dTilde^x_xx + A_x - GammaTilde_x + 2 Phi_x)^TF
...LambdaTilde^x_xx = (A_x + 2 Phi_x)^TF

... paper says A_x,t + partial_x (alpha f K) = 0 ... so the alpha is on the inside instead of the outside ...
... paper is on Bona-Masso, and the book's entry on BM matches the book's entry on BSSNOK but not the paper's BM ...
the book has no source terms whatsoever, though refers to everything in terms of ~= which neglects some terms (like the source terms maybe?) 

... with 1D simplifications and state variables removed

alpha,t = -alpha^2 f K						<- p.165 or so, which does match the paper, and does include a source term
phi,t = -alpha K / 6
A_x,t + alpha partial_x Q = 0
...A_x,t + alpha partial_x (f K) = 0
...A_x,t + alpha K f,alpha partial_x alpha + alpha f partial_x K = 0
Phi_x,t + alpha/6 partial_x K = 0
K,t + alpha exp(-4 phi) partial_x A_x = 0
ATilde_xx,t + alpha exp(-4 phi) partial_x LambdaTilde^x_xx = 0
...ATilde_xx,t + alpha exp(-4 phi) partial_x (A_x + 2 Phi_x) = 0
...ATilde_xx,t + alpha exp(-4 phi) partial_x A_x + 2 alpha exp(-4 phi) partial_x Phi_x = 0

variables:
alpha, phi, A_x, Phi_x, K, ATilde_xx

	[  alpha  ]   [0,               0, 0,                 0,                   0,       0]     [  alpha  ]   [-alpha^2 f K]
	[   phi   ]   [0,               0, 0,                 0,                   0,       0]     [   phi   ]   [-alpha K / 6]
d_t [   A_x   ]   [alpha K f,alpha, 0, 0,                 0,                   alpha f, 0] d_x [   A_x   ]   [     0      ]
    [  Phi_x  ] + [0,               0, 0,                 0,                   alpha/6, 0]     [  Phi_x  ] = [     0      ]
    [    K    ]   [0,               0, alpha exp(-4 phi), 0,                   0,       0]     [    K    ]   [     0      ]
    [ATilde_xx]   [0,               0, alpha exp(-4 phi), 2 alpha exp(-4 phi), 0,       0]     [ATilde_xx]   [     0      ]
/* [wxMaxima: input   start ] */
load("eigen");
/* [wxMaxima: input   end   ] */
/* [wxMaxima: input   start ] */
assume(alpha>0, f>0);
/* [wxMaxima: input   end   ] */
/* [wxMaxima: input   start ] */
vars : [alpha, phi, A_x, Phi_x, K, ATilde_xx];
/* [wxMaxima: input   end   ] */
/* [wxMaxima: input   start ] */
A : matrix(
[0,0,0,0,0,0],
[0,0,0,0,0,0],
[alpha*K*df_dalpha,0,0,0,alpha*f,0],
[0,0,0,0,1/6*alpha,0],
[0,0,alpha*exp(-4*phi),0,0,0],
[0,0,alpha*exp(-4*phi),2*alpha*exp(-4*phi),0,0]
);
/* [wxMaxima: input   end   ] */
/* [wxMaxima: input   start ] */
determinant(A);
/* [wxMaxima: input   end   ] */
/* [wxMaxima: input   start ] */
'diff(transpose(vars),t) + A . 'diff(transpose(vars),x) =
transpose(matrix([
-alpha^2 * f * K,
-2*alpha*ATilde_xx,
0,0,0,0
]));
/* [wxMaxima: input   end   ] */
/* [wxMaxima: input   start ] */
results:ev(similaritytransform(A),hermitianmatrix=true);
/* [wxMaxima: input   end   ] */
/* [wxMaxima: input   start ] */
eigenvalues : diag_matrix(results[1][1][1], results[1][1][3], results[1][1][3], results[1][1][3], results[1][1][3], results[1][1][2]);
/* [wxMaxima: input   end   ] */
/* [wxMaxima: input   start ] */
results[2][3];
/* [wxMaxima: input   end   ] */
/* [wxMaxima: input   start ] */
eigenvectors : transpose(matrix(
results[2][1][1] * denom(results[2][1][1][3]),
[1,0,0,0,0,0],
results[2][3][1],
[0,0,0,6*sqrt(f)*exp(-2*phi),1,0],
results[2][3][2],
results[2][2][1] * denom(results[2][2][1][3])
));
/* [wxMaxima: input   end   ] */
/* [wxMaxima: input   start ] */
transpose(eigenvectors) . eigenvectors$ ratsimp(%)$ ratsimp(%);
subst([f=-1/3],%)$ratsimp(%);
/* [wxMaxima: input   end   ] */
/* [wxMaxima: input   start ] */
determinant(eigenvectors)$ ratsimp(%);
subst([f=-1/3],%)$ratsimp(%);
/* [wxMaxima: input   end   ] */
/* [wxMaxima: input   start ] */
invert(eigenvectors);
/* [wxMaxima: input   end   ] */
--]]

local class = require 'ext.class'
local Equation = require 'equation'

local BSSNOK1D = class(Equation)
BSSNOK1D.name = 'BSSNOK 1D Hyperbolic'
BSSNOK1D.numStates = 6

function BSSNOK1D:init(args, ...)
	local symmath = require 'symmath'
	local function makesym(field)
		return symmath.clone(assert(args[field], "expected to find field "..field)):simplify() 
	end
	
	local x = assert(args.x)
	
	-- parameters that are symbolic functions -- based on coordinates 
	local stateExprs = table{'alpha', 'gamma_xx', 'K_xx'}:map(function(name)
		return makesym(name), name
	end)

	-- derived functions
	stateExprs.A_x = (stateExprs.alpha:diff(x) / stateExprs.alpha):simplify()
	stateExprs.phi = -symmath.log(stateExprs.gamma_xx)/4
	stateExprs.Phi_x = stateExprs.phi:diff(x):simplify()
	stateExprs.K = (stateExprs.K_xx / stateExprs.gamma_xx):simplify()
	stateExprs.ATilde_xx = symmath.exp(-4 * stateExprs.phi) * (stateExprs.K_xx - stateExprs.K/3 * stateExprs.gamma_xx)
	
	-- convert from symbolic functions to Lua functions
	self.calc = stateExprs:map(function(expr, name)
		return expr:compile{x}, name
	end)
	
	-- parameters that are symbolic functions -- based on f's alpha
	-- NOTICE: assign these *after* mapping stateExprs to calc
	local f = makesym'f'
	local f_param = assert(args.f_param)
	self.calc.f = f:compile{f_param}
	
	local dalpha_f = f:diff(args.f_param):simplify()
	self.calc.dalpha_f = dalpha_f:compile{args.f_param}
end

--phi = -1/(4*n) ln gamma_xx
-- exp(-4n phi) = gamma_xx for n=1
-- volume = sqrt(gamma_xx) = sqrt(exp(-4n phi)) = exp(-2n phi)
-- but let's use n=3 regardless of this being a 1D sim.  I think that's the rule of thumb?
do
	-- TODO fix this.  gamma_xx and A_x are negatived ... 
	-- I think the eigenvalues are negative ...
	-- BIGGER TODO - stop using your own calculated eigenvectors from your own calculated flux. just use the book. that will fix a lot of things.
	local q = function(self,i) return self.qs[i] end
	local alpha = q:_(1)
	local phi = q:_(2)
	local gamma_xx = math.exp:o(-12 * phi)
	local A_x = q:_(3)
	local Phi_x = q:_(4)
	local K = q:_(5)
	local ATilde_xx = q:_(6)
	local volume = alpha * math.sqrt:o(gamma_xx)
	BSSNOK1D:buildGraphInfos{
		{alpha = alpha},
		{A_x = A_x},
		{gamma_xx = gamma_xx},
		{phi = phi},
		{Phi_x = Phi_x},
		{K = K},
		{ATilde_xx = ATilde_xx},
		{volume = volume},
	}
end

function BSSNOK1D:initCell(sim,i)
	local x = sim.xs[i]
	local alpha = self.calc.alpha(x)
	local phi = self.calc.phi(x)
	local A_x = self.calc.A_x(x)
	local Phi_x = self.calc.Phi_x(x)
	local K = self.calc.K(x)
	local ATilde_xx = self.calc.ATilde_xx(x)
	return {alpha, phi, A_x, Phi_x, K, ATilde_xx}
end

function BSSNOK1D:calcInterfaceEigenBasis(sim,i,qL,qR)
	local avgQ = {}
	for j=1,self.numStates do
		avgQ[j] = (qL[j] + qR[j]) / 2
	end
	
	local alpha, phi, A_x, Phi_x, K, ATilde_xx = unpack(avgQ)
	local f = self.calc.f(alpha)
	local dalpha_f = self.calc.dalpha_f(alpha)
	
	local e2p = math.exp(2*phi)
	local ie2p = 1/e2p
	local ie4p = ie2p * ie2p
	local f_1_2 = math.sqrt(f)
	local f_3_2 = f * f_1_2
	local lambda = alpha * f_1_2 * ie2p
	sim.eigenvalues[i] = {-lambda, 0, 0, 0, 0, lambda}
	-- row-major, math-indexed
	sim.fluxMatrix[i] = {
		{0,0,0,0,0,0},
		{0,0,0,0,0,0},
		{alpha * dalpha_f * K, 0, 0, 0, alpha * f, 0},
		{0, 0, 0, 0, alpha/6, 0},
		{0, 0, alpha * ie4p, 0, 0, 0},
		{0, 0, alpha * ie4p, 2 * alpha * ie4p, 0, 0},
	}
	sim.eigenvectors[i] = {
		{0, 1, 0, 0, 0, 0},
		{0, 0, 1, 0, 0, 0},
		{6*f_3_2*e2p, 0, 0, 0, 0, 6*f_3_2*e2p},
		{f_1_2*e2p, 0, 0, 6*f_1_2*ie2p, 0, f_1_2*e2p},
		{-6*f, 0, 0, 1, 0, 6*f},
		{-6*f-2, 0, 0, 0, 1, 6*f+2},
	}
	local tmp1 = e2p/(36*f_3_2)+ie2p/f_1_2
	sim.eigenvectorsInverse[i] = {
		{0, 0, ie2p/(6*f_3_2)*(1 - .5*(f_1_2*e2p*tmp1)), e2p/(72*f_3_2), -1/(12*f), 0},
		{1,0,0,0,0,0},
		{0,1,0,0,0,0},
		{0,0,-e2p/(36*f_3_2),e2p/(6*f_1_2),0,0},
		{0,0,-(12*f+4)*tmp1/(12*f) - (-6*f-2)*ie2p/(6*f_3_2), (12*f+4)*e2p/(72*f_3_2), -(12*f+4)/(12*f), 1},
		{0, 0, tmp1/(12*f), -e2p/(72*f_3_2), 1/(12*f), 0},
	}
end

function BSSNOK1D:calcEigenvaluesFromCons(alpha, phi, ...)
	local f = self.calc.f(alpha)
	local f_1_2 = math.sqrt(f)
	local e2p = math.exp(2*phi)
	local ie2p = 1/e2p
	local lambda = alpha * f_1_2 * ie2p
	return -lambda, 0, 0, 0, 0, lambda
end

function BSSNOK1D:sourceTerm(sim, qs)
	local source = sim:newState()
	for i=1,sim.gridsize do
		local alpha, phi, A_x, Phi_x, K, ATilde_xx = unpack(qs[i])
		local f = self.calc.f(alpha)
		source[i][1] = -alpha * alpha * f * K
		source[i][2] = -alpha * K / 6
	end
	return source
end

return BSSNOK1D
