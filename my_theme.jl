# ref: https://ai.googleblog.com/2019/08/turbo-improved-rainbow-colormap-for.html
turboColorMap = cgrad(["#331b3d", "#4d6edf", "#3eb9e9", "#78ed9b","#9de94f","#ebce4c","#ce3b2a","#781e14"], [0.51, 1.5, 2.53, 3.49, 4.49, 5.52, 6.51, 7.5]./8)
blueColorMap = cgrad(["#00517b", "#0080a1", "#00aec6", :white])
vik10ColorMap = cgrad(["#001261", "#033e7d", "#1e6f9c", "#71a7c3", "#c8dce6", "#e9cdbc", "#d29674", "#bd6533", "#8a2706", "#590008"], 15, categorical = false)
# RGYB = cgrad(["#d73e2c", "#e94f31", "#f07860",  "#468f5a", "#5dba74", "#71e092", "#e39a30", "#f0af34", "#fac762", "#0171c1", "#0091f8", "#34a8fb"])
RGYB = cgrad(["#f07860",   "#71e092", "#fac762",  "#34a8fb"])
goodRed = RGBf(255/255,0.0,98/255)
backgroundcolor = "#06001A"
my_colors = [:white, "#5BD9DE", "#F75E85", "#F7A85E",  "#5EF793", "#b59edc"]

my_theme = Theme(
    backgroundcolor = "#06001A",
    colormap = turboColorMap,
	Axis = (
        backgroundcolor = "#06001A",
		xgridcolor = :grey, xlabelcolor = :white,
		ygridcolor = :grey, ylabelcolor = :white,
		xtickcolor = :white, xticklabelcolor = :white,
		ytickcolor = :white, yticklabelcolor = :white,
        titlecolor=:white,
        topspinecolor=:white, bottomspinecolor=:white,
        leftspinecolor=:white, rightspinecolor=:white,
        xgridwidth=0.2,ygridwidth=0.2,
        xminortickcolor=:white,
        yminortickcolor=:white,
        yminorticks=IntervalsBetween(9),
        xminorticks=IntervalsBetween(9),
        xlabelsize=20, ylabelsize=25,
        xticklabelsize=20, yticklabelsize=20,
	),
	Axis3 = (
		xgridcolor = :grey, xlabelcolor = :white,
		ygridcolor = :grey, ylabelcolor = :white,
		zgridcolor = :grey, zlabelcolor = :white,
		xtickcolor = :white, xticklabelcolor = :white,
		ytickcolor = :white, yticklabelcolor = :white,
        ztickcolor = :white, zticklabelcolor = :white,
        xspinecolor_1=:white,xspinecolor_2=:white,xspinecolor_3=:white,
        yspinecolor_1=:white,yspinecolor_2=:white,yspinecolor_3=:white,
        zspinecolor_1=:white,zspinecolor_2=:white,zspinecolor_3=:white,
        titlecolor=:white,
        xlabelsize=20, ylabelsize=25,
        xticklabelsize=20, yticklabelsize=20,
	),
	Colorbar = (labelcolor = :white, tickcolor = :white, minortickcolor = :white,
		rightspinecolor = :white, ticklabelcolor = :white),
    palette = (color = cgrad(my_colors),
                linestyle=[:solid, :dash, :dot,:dashdot], marker = [:circle, :diamond, :rect]),

    Lines = (cycle = Cycle([:color, :linestyle], covary = false), linewidth=5,),
    Scatter = (cycle = Cycle([:marker, :color], covary = true),),

    Legend = (bgcolor="#06001A", framecolor=:white, labelcolor=:white, titlecolor=:white),
    # Barplot = (label_color = :white) # does not work ?
)


myWhite_theme = Theme(
    Axis = (
        xgridwidth=0.2,ygridwidth=0.2,
        yminorticks=IntervalsBetween(9),
        xminorticks=IntervalsBetween(9),
        xlabelsize=20, ylabelsize=25,
        xticklabelsize=20, yticklabelsize=20,
	),
    colormap = vik10ColorMap,
    palette = (color = cgrad(["#684c41", "#4fa3a5", "#fdae38", "#f75435","#45c5f7"]),
        linestyle=[:solid, :dash, :dot, :dashdot],  marker = [:circle, :diamond, :rect]),
    Lines = (cycle = Cycle([:color, :linestyle], covary = false), linewidth=5,),
    Scatter = (cycle = Cycle([:marker, :color], covary = true),),
)