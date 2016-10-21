
# Paths in Paths!
# See here for spec / tutorial: https://developer.mozilla.org/en/docs/Web/SVG/Tutorial/Paths

svgf=open("path.svg","w")

@printf(svgf,"<svg width=\"200\" height=\"200\" xmlns=\"http://www.w3.org/2000/svg\">\n")

#//  <path d="M10 10"/>

function printstartpath(x,y)
    @printf(svgf,"<path d=\"M %f %f ",x,y)
end

function printsegment(dx,dy)
    @printf(svgf,"<path d=\"M l %f %f",dx,dy)
end

function printendpath()
    @printf(svgf,"\"\\>")
end

function drawpath()
    printstartpath(0,0)
    printsegment(1,1)
    printsegment(1,0)
    printsegment(1,1)
    printendpath()
end

drawpath()

@printf(svgf,"\n</svg>\n")
close(svgf)
