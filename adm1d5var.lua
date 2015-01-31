--[[
hyperbolic formalism:

state vector: [alpha g A D K]	<- even though alpha and g are already represented by A and D ...
fluxes vector: alpha K * [0, 0, f/g, 1, A/K]
source vector: alpha / g * [-alhpa f K, -2 K g, 0, 0, A D - K^2]

alpha,t + (0),x = -alpha^2 f K / g
g,t + (0),x = -2 alpha K
A,t + (alpha K f / g),x = 0
D,t + (alpha K),x = 0
K,t + (alpha A),x = alpha / g (A D - K^2)

d_t alpha + (0),x = -alpha^2 f K / g
d_t g + (0),x = -2 alpha K
d_t A + K f / g alpha,x + alpha K / g f,x + alpha f / g K,x - alpha K f / g^2 g,x = 0 
d_t D + K alpha,x + alpha K,x = 0
d_t K + A alpha,x + alpha A,x = alpha / g (A D - K^2)

    [alpha]   [   0,            0,          0,   0,      0     ]     [alpha]   [  -alpha^2 f K / g   ]
    [  g  ]   [   0,            0,          0,   0,      0     ]     [  g  ]   [     -2 alpha K      ]
d_t [  A  ] + [f K / g, -alpha f K / g^2,   0,   0, alpha f / g] d_x [  A  ] = [          0          ]
    [  D  ]   [   K,            0,          0,   0,    alpha   ]     [  D  ]   [          0          ]
    [  K  ]   [   A,            0,        alpha, 0,      0     ]     [  K  ]   [alpha / g (A D - K^2)]

eigenvalues of A:
-lambda * (
	-lambda * (
		-lambda * det(
			-lambda		alpha f / g
			alpha		-lamda
		)
	)
) = 0
<=>
lambda = 0 has multiplicity 3,

lambda^2 - alpha^2 f / g = 0 
lambda = +-alpha sqrt(f / g)

for lambda = -alpha sqrt(f/g) the eigenvector is 
	[0, 0, 1, g/f, -sqrt(g/f)] dot state = A + D g/f - K sqrt(g/f)
	= A + D f / g - K * sqrt(f / g)
for lambda = +alpha sqrt(f/g) the eigenvector is 
	[0, 0, 1, g/f, +sqrt(g/f)] dot state = A + D g/f + K sqrt(g/f)
for lambda = 0 the eigenvectors are
[alpha,0,-A,0,-K]
[0,0,0,1,0]
[0,1,0,0,0]

eigenvector matrix is :
[	0			alpha	0	0	0			]
[	0			0		0	1	0			]
[	f/g			-A		0	0	f/g			]
[	1			0		1	0	1			]
[	-sqrt(f/g)	-K		0	0	sqrt(f/g)	]
inverse:
[(sqrt(f)*g^(3/2)*A-f*g*K)/(2*alpha*f^(3/2)*sqrt(g)),0,g/(2*f),0,-sqrt(g)/(2*sqrt(f))]
[1/alpha,0,0,0,0]
[-(g*A)/(alpha*f),0,-g/f,1,0]
[0,1,0,0,0]
[(sqrt(f)*sqrt(g)*K+g*A)/(2*alpha*f),0,g/(2*f),0,sqrt(g)/(2*sqrt(f))]

/* [wxMaxima: input   start ] */
fluxMatrix : matrix(
[0,0,0,0,0],
[0,0,0,0,0],
[f*K/g, -alpha*f*K/g^2, 0,0, alpha*f/g],
[K,0,0,0,alpha],
[A,0,alpha,0,0]);
/* [wxMaxima: input   end   ] */
/* [wxMaxima: input   start ] */
load("eigen");
/* [wxMaxima: input   end   ] */
/* [wxMaxima: input   start ] */
assume(f>0,g>0,alpha>0);
/* [wxMaxima: input   end   ] */
/* [wxMaxima: input   start ] */
result:ev(similaritytransform(fluxMatrix), hermitianmatrix=true);
/* [wxMaxima: input   end   ] */
/* [wxMaxima: input   start ] */
eigenvalueMatrix : diag_matrix(
result[1][1][1],
result[1][1][3], 
result[1][1][3],
result[1][1][3],
result[1][1][2]);
/* [wxMaxima: input   end   ] */
/* [wxMaxima: input   start ] */
transpose(matrix(
result[2][1][1] * sqrt(f^2 + f*g + g^2)/g,
result[2][3][1] * sqrt(K^2 + A^2 + alpha^2),
result[2][3][2],
[0,1,0,0,0],
result[2][2][1] * sqrt(f^2 + f*g + g^2)/g))$
ratsimp(%)$
eigenvectorMatrix : %;
/* [wxMaxima: input   end   ] */
/* [wxMaxima: input   start ] */
invert(eigenvectorMatrix)$
ratsimp(%)$
eigenvectorInverseMatrix : %;
/* [wxMaxima: input   end   ] */
/* [wxMaxima: input   start ] */
eigenvectorMatrix . eigenvectorInverseMatrix$
ratsimp(%);
/* [wxMaxima: input   end   ] */
/* [wxMaxima: input   start ] */
eigenvectorMatrix . eigenvalueMatrix . eigenvectorInverseMatrix $
ratsimp(%);
is(% = fluxMatrix);
fluxMatrix;
/* [wxMaxima: input   end   ] */

--]]

require 'ext'
local Simulation = require 'simulation'

local ADM1D5VarSim = class(Simulation)
	
ADM1D5VarSim.numStates = 5 

function ADM1D5VarSim:init(args, ...)
	Simulation.init(self, args, ...)

	local symmath = require 'symmath'

	local x = assert(args.x)

	local h = symmath.clone(assert(args.h)):simplify()
	self.calc_h = h:compile{x}
	
	local dx_h = h:diff(x):simplify()
	self.calc_dx_h = dx_h:compile{x}
	
	local d2x_h = dx_h:diff(x):simplify()
	self.calc_d2x_h = d2x_h:compile{x}

	local g = symmath.clone(assert(args.g)):simplify()
	self.calc_g = g:compile{x}

	local dx_g = g:diff(x):simplify()
	self.calc_dx_g = dx_g:compile{x}

	local alpha = symmath.clone(assert(args.alpha)):simplify()
	self.calc_alpha = alpha:compile{x}

	local dx_alpha = alpha:diff(x):simplify()
	self.calc_dx_alpha = dx_alpha:compile{x}

	local f = symmath.clone(assert(args.alpha)):simplify()
	self.calc_f = f:compile{assert(args.f_var)}

	self.graphInfos = {
		{viewport={0/3, 0/3, 1/3, 1/3}, getter=index:bind(self.qs):index(1), name='alpha', color={1,0,1}},
		{viewport={0/3, 1/3, 1/3, 1/3}, getter=index:bind(self.qs):index(3), name='A', color={0,1,0}},
		{viewport={1/3, 0/3, 1/3, 1/3}, getter=index:bind(self.qs):index(2), name='g', color={.5,.5,1}},
		{viewport={1/3, 1/3, 1/3, 1/3}, getter=index:bind(self.qs):index(4), name='D', color={1,1,0}},
		{viewport={2/3, 0/3, 1/3, 1/3}, getter=index:bind(self.qs):index(5), name='K', color={0,1,1}},
		{viewport={0/3, 2/3, 1/3, 1/3}, getter=log:compose(index:bind(self.eigenbasisErrors)), name='error', color={1,0,0}, range={-30, 30}},
		{viewport={1/3, 2/3, 1/3, 1/3}, getter=log:compose(index:bind(self.fluxMatrixErrors)), name='error', color={1,0,0}, range={-30, 30}},
	}
end
	
function ADM1D5VarSim:initCell(i)
	local x = self.xs[i]
	local alpha = self.calc_alpha(x)
	local g = self.calc_g(x)
	local A = self.calc_dx_alpha(x) / self.calc_alpha(x)
	local D = 1/2 * self.calc_dx_g(x)
	local K = -self.calc_d2x_h(x) / sqrt(self.calc_g(x))
	return {alpha, g, A, D, K}
end

function ADM1D5VarSim:calcInterfaceEigenBasis(i)
	local avgQ = {}
	for j=1,self.numStates do
		avgQ[j] = (self.qs[i-1][j] + self.qs[i][j]) / 2
	end
	
	local alpha, g, A, D, K = unpack(avgQ)
	local x = self.ixs[i]
	local f = self.calc_f(alpha)
	local lambda = alpha * sqrt(f / g)		
	self.eigenvalues[i] = {-lambda, 0, 0, 0, lambda}
	-- row-major, math-indexed
	self.fluxMatrix[i] = {
		{0,0,0,0,0},
		{0,0,0,0,0},
		{f*K/g, -alpha*f*K/g^2, 0,0, alpha*f/g},
		{K,0,0,0,alpha},
		{A,0,alpha,0,0},
	}
	self.eigenvectors[i] = {
		{0,			alpha,	0,	0,	0			},	-- alpha
		{0,			0,		0,	1,	0			},	-- g
		{f/g,		-A,		0,	0,	f/g			},	-- A
		{1,			0,		1,	0,	1			},	-- D
		{-sqrt(f/g),-K,		0,	0,	sqrt(f/g)	},	-- K
	}
	self.eigenvectorsInverse[i] = {
		{(g * A / f - K * sqrt(g / f)) / (2 * alpha), 0, g / (2 * f), 0, -.5 * sqrt(g / f)}, 
		{1 / alpha, 0, 0, 0, 0}, 
		{-(g * A) / (alpha * f), 0, -g / f, 1, 0}, 
		{0, 1, 0, 0, 0}, 
		{(g * A / f + K * sqrt(g / f)) / (2 * alpha), 0, g / (2 * f), 0, .5 * sqrt(g / f)}, 
	}
	-- note that because we have zero eigenvalues that the eigendecomposition cannot reconstruct the flux matrix
end

function ADM1D5VarSim:addSourceToDerivCell(i)
	local alpha, g, A, D, K = unpack(self.qs[i])
	local f = self.calc_f(alpha)
	self.dq_dts[i][1] = self.dq_dts[i][1] - alpha * alpha * f * K / g
	self.dq_dts[i][2] = self.dq_dts[i][2] - 2 * alpha * K
	self.dq_dts[i][5] = self.dq_dts[i][5] + alpha * (A * D - K * K) / g
end

return ADM1D5VarSim

