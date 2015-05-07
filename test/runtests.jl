
using PyPlot
include("../src/LDPCsim.jl")
using LDPCsim


ldpcSim(1000,[-2.0:1:10],[50])
savefig("test/ldpcVsSNR_M50rHalf.png") # run from root directory

ldpcSim(1000,[3.5],[4,8,16,32,64,128])
savefig("test/ldpcVsM_snr3p5_rHalf.png") # run from root directory
