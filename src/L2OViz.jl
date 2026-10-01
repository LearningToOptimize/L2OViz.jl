module L2OViz

using CairoMakie
using Statistics

include("select_entries.jl")
include("_utils.jl")
include("plot_variable.jl")
include("plot_graph_variable.jl")
include("plot_scale_histograms.jl")
include("animate_variable.jl")
include("animate_graph_variable.jl")

export plot_variable, plot_graph_variable, plot_scale_histograms,
       animate_variable, animate_graph_variable

end
