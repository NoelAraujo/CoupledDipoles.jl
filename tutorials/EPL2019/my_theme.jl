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

function convert_hex_rhb(hex_color)
    r = parse(Int, hex_color[1:2], base=16)/255
    g = parse(Int, hex_color[3:4], base=16)/255
    b = parse(Int, hex_color[5:6], base=16)/255
    return RGBf(r,g,b)
end


myHex1 = ["03071E","370617","6A040F","9D0208","D00000","DC2F02","E85D04","F48C06","FAA307","FFBA08"]
myGrad1 = cgrad(convert_hex_rhb.(myHex1))
myColor4 = convert_hex_rhb("01FA70")

## More Options
# myHex2 = ["00072D", "001C55", "0A2472", "0E6BA8", "A6E1FA"]
# myGrad2 = cgrad(convert_hex_rhb.(myHex2))
# myColor4 = convert_hex_rhb("F4B000")

# myHex3 = ["001219", "9B2226", "EE9B00" ]
# myGrad3 = cgrad(convert_hex_rhb.(myHex3))
# myColor4 = convert_hex_rhb("01FA70")

myColors = cgrad(myGrad1, 3, categorical=true, rev=false)
myColor1 = myColors[1]
myColor2 = myColors[2]
myColor3 = myColors[3]


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
