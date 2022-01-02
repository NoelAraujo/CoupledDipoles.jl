using PrettyNumbers
function get_my_yticks(numero)
    io = IOBuffer()
    pretty_number(io, numero, show_significand = false)
    nice_number = io |> take! |> String
    return nice_number
end

function add_box_inset(
    fig;
    left = 100,
    right = 250,
    bottom = 200,
    top = 300,
    bgcolor = :white,
)
    inset_box = Axis(
        fig,
        bbox = BBox(left, right, bottom, top),
        xlabelsize = 20,
        ylabelsize = 20,
        xlabel = L"|r - r_{cm}|",
        ylabel = L"log_{10} \; |\psi|^2",
        xticklabelsize = 12,
        yticklabelsize = 12,
        backgroundcolor = bgcolor,
    )
    # bring content upfront
    translate!(inset_box.scene, 0, 0, 10)
    elements = keys(inset_box.elements)
    filtered = filter(ele -> ele != :xaxis && ele != :yaxis, elements)
    foreach(ele -> translate!(inset_box.elements[ele], 0, 0, 9), filtered)
    return inset_box
end

myColor1 = RGBf(1 / 255, 148 / 255, 154 / 255) # Teal
myColor2 = RGBf(0, 67 / 255, 105 / 255) # Navy_Blue
myColor3 = RGBf(219 / 255, 31 / 255, 72 / 255) # 

Fig1_theme = Theme(
    Axis = (
        xlabelsize = 30,
        ylabelsize = 35,
        xticksize = 10,
        yticksize = 10,
        xtickalign = 1,
        ytickalign = 1,
        xlabelpadding = -5,
        xticks = LinearTicks(3),
    ),
    palette = (
        color = [(myColor1, 1.0), (myColor2, 1.0), (myColor3, 1.0)],
        marker = [:utriangle, :circle, :rect],
    ),
    Scatter = (cycle = Cycle([:marker, :color], covary = true),),
)

Fig2_theme = Theme(
    Axis = (
        xlabelsize = 30,
        ylabelsize = 35,
        xticksize = 10,
        yticksize = 10,
        xtickalign = 1,
        ytickalign = 1,
        xlabelpadding = -5,
        xticks = LinearTicks(3),
    ),
    palette = (
        color = [(myColor3, 1.0), (myColor2, 1.0), (myColor1, 1.0)],
        marker = [:rect, :circle, :utriangle],
    ),
    Scatter = (cycle = Cycle([:marker, :color], covary = true),),
)


spatialProfile_theme = Theme(
    Axis = (
        xlabelsize = 30,
        ylabelsize = 35,
        xticksize = 10,
        yticksize = 10,
        xtickalign = 1,
        ytickalign = 1,
        xlabelpadding = -5,
        yscale = log10,
        xlabel = L"|r - r_{cm}|",
        ylabel = L"|\psi|^2",
    ),
    palette = (color = [myColor1, myColor2, myColor3], marker = [:utriangle, :circle, :rect]),
    Scatter = (cycle = Cycle([:marker, :color], covary = true),),
)
