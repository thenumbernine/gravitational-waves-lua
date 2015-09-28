local class = require 'ext.class'
local table = require 'ext.table'

local ADM3to5Var = require 'adm1d3to5var'
local RoeBackwardEulerLinear = require 'roe_backwardeuler_linear'

local ADM3to5VarRoeBackwardEulerLinear = class(RoeBackwardEulerLinear)

ADM3to5VarRoeBackwardEulerLinear.name = ADM3to5VarRoeBackwardEulerLinear

function ADM3to5VarRoeBackwardEulerLinear:init(args)
	args = table(args)
	args.equation = ADM3to5Var(args)
	ADM3to5VarRoeBackwardEulerLinear.super.init(self, args)
end

return ADM3to5VarRoeBackwardEulerLinear
