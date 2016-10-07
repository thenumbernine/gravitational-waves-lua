local class = require 'ext.class'
local roe = require 'roe'
local RoeNewtonKrylov = class(Roe)

function RoeNewtonKrylov:init(args)
	RoeNewtonKrylov.super.init(self, args)
	self.linearSolver = args.linearSolver or require 'linearsolvers'.conjres
	self.errorLogging = args.errorLogging
end

function 
