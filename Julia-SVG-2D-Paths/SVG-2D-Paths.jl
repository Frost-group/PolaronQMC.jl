# Feynman 1D path-integral Paths drawn with SVG Paths!

# Paths in Paths!
# See here for spec / tutorial: https://developer.mozilla.org/en/docs/Web/SVG/Tutorial/Paths

svgf = open("path.svg", "w")

tSlices = 100 # Also used as size of SVG canvas

@printf(
    svgf,
    "<svg width=\"%d\" height=\"%d\" xmlns=\"http://www.w3.org/2000/svg\">\n",
    tSlices,
    tSlices
)

#//  <path d="M10 10"/>

function printstartpath(x, y)
    @printf(svgf, "<path d=\"M %d %d ", x, y)
end

function printsegment(dx, dy)
    @printf(svgf, " l %d %d", dx, dy)
    #    @printf("printsegment: %d %d\n",dx,dy)
end

function printendpath()
    colours = ["#7fc97f", "#beaed4", "#fdc086", "#ffff99", "#386cb0", "#f0027f"] # Colourbrewer 6 colours qualitative
    colour = rand(colours)

    @printf(svgf, "  \" stroke=\"%s\" fill=\"none\" stroke-width=\"0.5\"  \/>", colour)
end

function drawpath()
    printstartpath(tSlices / 2.0, 0)

    dt = 1
    sx = 0
    x1 = 5 # This target, but it tends to make the paths look rather artificial!

    # Really need to compose this random vector so paths connect from [x0,t0]->[x1,t1]
    for t = 0:dt:tSlices
        dx = rand([1, -1]) # Forces either a (1,1) or (1,-1) move; zizzag space-time paths

        if (tSlices - t - 1) < abs(sx - x1) # Force a reconnect in the space-time paths
            dx = -(sx - x1) / abs(sx - x1) # Head in this direction only.
        end
        printsegment(dx, dt)
        sx += dx # keep track of where we are, to reconnect space-time paths
    end

    printendpath()
end

function drawpathPermutations()
    printstartpath(tSlices / 2.0, 0)

    # x1 is an offset for the bottom of the path. I can't remember why I added this!
    x1 = 10

    dxs = shuffle([fill(1, x1 + div(tSlices, 2)); fill(-1, -x1 + div(tSlices, 2))])
    # Shuffle set of [1,...,1] and [-1,...-1] arrays together, generating a random set of dxs which sum to zero
    # i.e. we have a path that connects from x0,t0 to x1,t1, by definition

    dt = 1
    for dx in dxs
        printsegment(dx, dt)
    end

    printendpath()
end

# A few random paths...
Paths = 50
for i = 1:Paths
    drawpathPermutations()
end

@printf(svgf, "\n</svg>\n")
@printf(svgf, "<svg width=\"200\" height=\"200\" xmlns=\"http://www.w3.org/2000/svg\">\n")
@printf(svgf, "<svg width=\"100\" height=\"100\" xmlns=\"http://www.w3.org/2000/svg\">\n")

#//  <path d="M10 10"/>

function printstartpath(x, y)
    @printf(svgf, "<path d=\"M %d %d ", x, y)
end

function printsegment(dx, dy)
    @printf(svgf, " l %d %d", dx, dy)
    #    @printf("printsegment: %d %d\n",dx,dy)
end

function printendpath()
    colours = ["#7fc97f", "#beaed4", "#fdc086", "#ffff99", "#386cb0", "#f0027f"] # Colourbrewer 6 colours qualitative
    colour = rand(colours)

    @printf(svgf, "  \" stroke=\"%s\" fill=\"none\" stroke-width=\"0.5\"  \/>", colour)
end

function drawpath()
    printstartpath(tSlices / 2.0, 0)

    dt = 1
    sx = 0
    # Really need to compose this random vector so paths connect from [x0,t0]->[x1,t1]
    for t = 0:dt:tSlices
        dx = rand([1, -1]) # Forces either a (1,1) or (1,-1) move; zizzag space-time paths

        if (tSlices - t - 1) < abs(sx - x1) # Force a reconnect in the space-time paths
            dx = -(sx - x1) / abs(sx - x1) # Head in this direction only.
        end
        printsegment(dx, dt)
        sx += dx # keep track of where we are, to reconnect space-time paths
    end

    printendpath()
end

function drawpathPermutations()
    printstartpath(tSlices / 2.0, 0)

    dxs = shuffle([fill(1, div(tSlices, 2)); fill(-1, div(tSlices, 2))])
    # Shuffle set of [1,...,1] and [-1,...-1] arrays together, generating a random set of dxs which sum to zero
    # i.e. we have a path that connects from x0,t0 to x1,t1, by definition

    dt = 1
    for dx in dxs
        printsegment(dx, dt)
    end

    printendpath()
end

# A few random paths...
Paths = 50
for i = 1:Paths
    drawpathPermutations()
end

@printf(svgf, "\n</svg>\n")
close(svgf)
