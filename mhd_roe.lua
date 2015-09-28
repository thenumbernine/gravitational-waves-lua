local class = require 'ext.class'
local Roe = require 'roe'
local MHD = require 'mhd'
local MHDRoe = class(Roe)
MHDRoe.equation = MHD()
return MHDRoe
