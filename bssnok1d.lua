--[[
going by Alcubierre's book, with some help from the 5-var system in his paper at http://arxiv.org/pdf/gr-qc/9609015v2.pdf 

p.85 non-hyperbolic description:
d/dt gammaTilde_ij = -2 alpha ATilde_ij
d/dt phi = -1/6 alpha K
d/dt ATilde_ij = exp(-4 phi) ( -D_i D_j alpha + alpha R_ij + 4 pi alpha (gamma_ij (S - rho) - 2 S_ij))^TF + alpha (K ATilde_ij - 2 ATilde_ik ATilde^k_j)
d/dt K = -D_i D^i alpha + alpha (ATilde_ij ATilde^ij + 1/3 K^2) + 4 pi alpha (rho + S)
... taking Bona-Masso lapse: d/dt alpha = -alpha^2 f K
... and no shift

alpha
g_xx
a_x = partial_x ln alpha = partial_x alpha / alpha
Phi_x = partial_x phi
dTilde_xxx = 1/2 partial_x gammaTilde_xx
K = K^x_x
ATilde_xx = exp(-4 phi) (K_xx - 1/3 K gamma_xx) = conformal trace-free extrinsic curvature
GammaTilde^x = gammaTilde^xx connTilde^x_xx

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
phi = -1/(4*n) ln gamma_xx
gammaTilde_xx = exp(-4 phi) gamma_xx 
...gammaTilde_xx = exp(1/n ln gamma_xx) gamma_xx 
...gammaTilde_xx = gamma_xx / gamma_xx^1
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
a_x,t + alpha partial_x Q = 0
...a_x,t + alpha partial_x (f K) = 0
...a_x,t + alpha K f,alpha partial_x alpha + alpha f partial_x K = 0
Phi_x,t + 1/6 alpha partial_x K = 0
dTilde_xxx,t + alpha partial_x ATilde_xx = 0
K,t + alpha exp(-4 phi) / gammaTilde_xx partial_x a_x = 0
ATilde_xx,t + alpha exp(-4 phi) partial_x LambdaTilde^x_xx = 0
GammaTilde^x,t + 4/3 alpha partial_x (K / gammaTilde_xx) = 0
...GammaTilde^x,t + 4/3 alpha / gammaTilde_xx partial_x K - 4/3 alpha K / gammaTilde_xx^2 partial_x gammaTilde_xx = 0

for LambdaTilde^x_xx = (dTilde^x_xx + a_x - GammaTilde_x + 2 Phi_x)^TF
...LambdaTilde^x_xx = (a_x + 2 Phi_x)^TF

... paper says a_x,t + partial_x (alpha f K) = 0 ... so the alpha is on the inside instead of the outside ...
... paper is on Bona-Masso, and the book's entry on BM matches the book's entry on BSSNOK but not the paper's BM ...
the book has no source terms whatsoever, though refers to everything in terms of ~= which neglects some terms (like the source terms maybe?) 

... with 1D simplifications and state variables removed

alpha,t = -alpha^2 f K						<- p.165 or so, which does match the paper, and does include a source term
phi,t = -alpha K / 6
a_x,t + alpha partial_x Q = 0
...a_x,t + alpha partial_x (f K) = 0
...a_x,t + alpha K f,alpha partial_x alpha + alpha f partial_x K = 0
Phi_x,t + alpha/6 partial_x K = 0
K,t + alpha exp(-4 phi) partial_x a_x = 0
ATilde_xx,t + alpha exp(-4 phi) partial_x LambdaTilde^x_xx = 0
...ATilde_xx,t + alpha exp(-4 phi) partial_x (a_x + 2 Phi_x) = 0
...ATilde_xx,t + alpha exp(-4 phi) partial_x a_x + 2 alpha exp(-4 phi) partial_x Phi_x = 0

variables:
alpha, phi, a_x, Phi_x, K, ATilde_xx

	[  alpha  ]   [0,               0, 0,                 0,                   0,       0]     [  alpha  ]   [-alpha^2 f K]
	[   phi   ]   [0,               0, 0,                 0,                   0,       0]     [   phi   ]   [-alpha K / 6]
d_t [   a_x   ]   [alpha K f,alpha, 0, 0,                 0,                   alpha f, 0] d_x [   a_x   ]   [     0      ]
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
vars : [alpha, phi, a_x, Phi_x, K, ATilde_xx];
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
require 'ext'
local Simulation = require 'simulation'
local BSSNOK1DSim = class(Simulation)
BSSNOK1DSim.numStates = 6

function BSSNOK1DSim:init(args, ...)
	BSSNOK1DSim.super.init(self, args, ...)
	local symmath = require 'symmath'
	local function makesym(x) return symmath.clone(assert(x)):simplify() end
	
	local x = assert(args.x)

	local alpha = makesym(args.alpha)
	self.calc_alpha = alpha:compile{x}
	local a_x = (alpha:diff(x) / alpha):simplify()
	self.calc_a_x = a_x:compile{x}
	
	local gamma_xx = makesym(args.g)
	local phi = -symmath.log(gamma_xx)/4
	self.calc_phi = phi:compile{x}
	local Phi_x = phi:diff(x):simplify()
	self.calc_Phi_x = Phi_x:compile{x}

	local h = makesym(args.h)
	local d2x_h = h:diff(x,x):simplify()
	local K_xx = -d2x_h / symmath.sqrt(gamma_xx)
	local K = (K_xx / gamma_xx):simplify()
	self.calc_K = K:compile{x}

	local ATilde_xx = symmath.exp(-4 * phi) * (K_xx - K/3 * gamma_xx)
	self.calc_ATilde_xx = ATilde_xx:compile{x}

	local f = makesym(args.alpha)
	self.calc_f = f:compile{x}
	
	local dalpha_f = f:diff(args.f_param):simplify()
	self.calc_dalpha_f = dalpha_f:compile{args.f_param}

	local get_state = index:bind(self.qs)
	local get_alpha = get_state:index(1)
	local get_phi = get_state:index(2)
	local get_a_x = get_state:index(3)
	local get_Phi_x = get_state:index(4)
	local get_K = get_state:index(5)
	local get_ATilde_xx = get_state:index(6)
	--phi = -1/(4*n) ln gamma_xx
	-- exp(-4n phi) = gamma_xx for n=1
	-- volume = sqrt(gamma_xx) = sqrt(exp(-4n phi)) = exp(-2n phi)
	self.graphInfos = {
		{viewport={0/3, 0/3, 1/3, 1/3}, getter=get_alpha, name='alpha', color={1,0,1}},
		{viewport={0/3, 1/3, 1/3, 1/3}, getter=get_a_x, name='a_x', color={0,1,0}},
		{viewport={1/3, 0/3, 1/3, 1/3}, getter=get_phi, name='phi', color={.5,.5,1}},
		{viewport={1/3, 1/3, 1/3, 1/3}, getter=get_Phi_x, name='Phi_x', color={1,1,0}},
		{viewport={2/3, 0/3, 1/3, 1/3}, getter=get_K, name='K', color={0,1,1}},
		{viewport={2/3, 1/3, 1/3, 1/3}, getter=get_ATilde_xx, name='ATilde_xx', color={0,1,1}},
		{viewport={0/3, 2/3, 1/3, 1/3}, getter=log:compose(index:bind(self.eigenbasisErrors)), name='eigenbasis error', color={1,0,0}, range={-30, 30}},
		{viewport={1/3, 2/3, 1/3, 1/3}, getter=log:compose(index:bind(self.fluxMatrixErrors)), name='reconstruction error', color={1,0,0}, range={-30, 30}},
		{viewport={2/3, 2/3, 1/3, 1/3}, getter=get_alpha * exp:compose(-2 * get_phi), name='volume', color={0,1,1}},
	}
end

function BSSNOK1DSim:initCell(i)
	local x = self.xs[i]
	local alpha = self.calc_alpha(x)
	local phi = self.calc_phi(x)
	local a_x = self.calc_a_x(x)
	local Phi_x = self.calc_Phi_x(x)
	local K = self.calc_K(x)
	local ATilde_xx = self.calc_ATilde_xx(x)
	return {alpha, phi, a_x, Phi_x, K, ATilde_xx}
end

function BSSNOK1DSim:calcInterfaceEigenBasis(i)
	local avgQ = {}
	for j=1,self.numStates do
		avgQ[j] = (self.qs[i-1][j] + self.qs[i][j]) / 2
	end
	
	local alpha, phi, a_x, Phi_x, K, ATilde_xx = unpack(avgQ)
	local f = self.calc_f(alpha)
	local dalpha_f = self.calc_dalpha_f(alpha)
	
	local e2p = exp(2*phi)
	local ie2p = 1/e2p
	local ie4p = ie2p * ie2p
	local f_1_2 = sqrt(f)
	local f_3_2 = f * f_1_2
	local lambda = -alpha * f_1_2 * ie2p
	self.eigenvalues[i] = {-lambda, 0, 0, 0, 0, lambda}
	-- row-major, math-indexed
	self.fluxMatrix[i] = {
		{0,0,0,0,0,0},
		{0,0,0,0,0,0},
		{alpha * dalpha_f * K, 0, 0, 0, alpha * f, 0},
		{0, 0, 0, 0, alpha/6, 0},
		{0, 0, alpha * ie4p, 0, 0, 0},
		{0, 0, alpha * ie4p, 2 * alpha * ie4p, 0, 0},
	}
	self.eigenvectors[i] = {
		{0, 1, 0, 0, 0, 0},
		{0, 0, 1, 0, 0, 0},
		{6*f_3_2*e2p, 0, 0, 0, 0, 6*f_3_2*e2p},
		{f_1_2*e2p, 0, 0, 6*f_1_2*ie2p, 0, f_1_2*e2p},
		{-6*f, 0, 0, 1, 0, 6*f},
		{-6*f-2, 0, 0, 0, 1, 6*f+2},
	}
	local tmp1 = e2p/(36*f_3_2)+ie2p/f_1_2
	self.eigenvectorsInverse[i] = {
		{0, 0, ie2p/(6*f_3_2)*(1 - .5*(f_1_2*e2p*tmp1)), e2p/(72*f_3_2), -1/(12*f), 0},
		{1,0,0,0,0,0},
		{0,1,0,0,0,0},
		{0,0,-e2p/(36*f_3_2),e2p/(6*f_1_2),0,0},
		{0,0,-(12*f+4)*tmp1/(12*f) - (-6*f-2)*ie2p/(6*f_3_2), (12*f+4)*e2p/(72*f_3_2), -(12*f+4)/(12*f), 1},
		{0, 0, tmp1/(12*f), -e2p/(72*f_3_2), 1/(12*f), 0},
	}
end

function BSSNOK1DSim:addSourceToDerivCell(dq_dts, i)
	local alpha, phi, a_x, Phi_x, K, ATilde_xx = unpack(self.qs[i])
	local f = self.calc_f(alpha)
	dq_dts[i][1] = dq_dts[i][1] - alpha * alpha * f * K
	dq_dts[i][2] = dq_dts[i][2] - alpha * K / 6
end

return BSSNOK1DSim

