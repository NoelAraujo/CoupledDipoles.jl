# using CoupledDipoles, Revise
using GLMakie;
Makie.inline!(true);

using ColorSchemes

mycg = cgrad(myGrad3, 10, categorical=true)
@time let
    Random.seed!(123)
    n = 15
    x, y, color = rand(n), rand(n), rand(n)
    cmaps = [mycg, :inferno, :cool, :viridis, :thermal]

    function FigGridHeatSharedCbarH(myData, titles)
        fig = Figure(resolution = (500, 600), fontsize = 20)
        g_1 = fig[1, 1] = GridLayout()
        g_2 = fig[1, 2] = GridLayout()
        
        minmax =  (0, 100)
        println(minmax)
        axes = []
        nData = length(myData)
        for (i, title) in zip(1:nData, titles)
            ax = Axis(
                g_1[i, 1],
                title = title,
                xlabel = L"\rho/k^3",
                ylabel = L"\Delta/\Gamma",
            )
            pnts = heatmap!(myData[i], colormap = cmaps[1], colorrange = minmax)
            ax.xticks = [1, 10]
            ax.yticks = [1, 10]
            cbar = Colorbar(
                fig,
                pnts,
                minorticksvisible=true,
                flipaxis = true,
                tickwidth = 2,
                ticklabelsize = 14,
                height = Relative(0.9),
            )
            cbar.ticks = [0, 50, 100]
            g_2[2:3, 1] = cbar
            push!(axes, ax)
        end

        # retain the xlabel only of the last figure
        [axes[j].xlabelvisible = false for j = 1:(nData-1)]

        Label(g_2[1, 1, Bottom()], L"P_{\omega}", textsize = 30, padding = (0, -10, -10, 0))

        rowsize!(g_2, 1, Relative(1 / 30))
        colgap!(fig.layout, 10)
        rowgap!(g_1, 0)
        rowgap!(g_2, 2)
        fig
    end

    myData = [33rand(10, 10) .+ 66, 33rand(10, 10) .+ 33, 33*rand(10, 10)]
    titles = ["super", "sub", "loc"]
    fig = FigGridHeatSharedCbarH(myData, titles)
    fig
end
