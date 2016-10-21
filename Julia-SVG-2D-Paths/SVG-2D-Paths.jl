
# Paths in Paths!
# See here for spec / tutorial: https://developer.mozilla.org/en/docs/Web/SVG/Tutorial/Paths

svgf=open("path.svg","w")

@printf(svgf,"<svg width=\"200\" height=\"200\" xmlns=\"http://www.w3.org/2000/svg\">\n")

#//  <path d="M10 10"/>

function printsegment(dx,dy)
    @printf(svgf,"<path d=\"M l %f %f",dx,dy)
end

printsegment(10,10)

function drawpath()
    printsegment(10,10)
end

close(svgf)
