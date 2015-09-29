local class = require 'ext.class'
local table = require 'ext.table'

local ADM1D3to5Var = require 'adm1d3to5var'

local Roe = require 'roe'

local ADM1D3to5VarRoe = class(Roe)

function ADM1D3to5VarRoe:init(args)
	args = table(args)
	args.equation = ADM1D3to5Var(args)
	ADM1D3to5VarRoe.super.init(self, args)
end

return ADM1D3to5VarRoe
