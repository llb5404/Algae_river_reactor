# the Revise package allows you to rerun code after making variable/function/method changes without re-compiling everything
using Revise
using Debugger
#https://pyco2sys.readthedocs.io/en/latest/co2sys_nd/#results

break_on(:error)

_revise_mode__ = :eval

include("AlgaeRiverReactor.jl")
using .AlgaeRiverReactor

Revise.track("AlgaeRiverReactor.jl")
Revise.track("LoadParameters.jl")
Revise.track("PDE_AlgaeBiomass.jl")
Revise.track("PDE_Temperature.jl")
Revise.track("PDE_CO2.jl")
Revise.track("MakePlots.jl")

@show Revise.watched_files

AlgaeRiverReactor.Run()
