
# Paths in Paths!
# See here for spec / tutorial: https://developer.mozilla.org/en/docs/Web/SVG/Tutorial/Paths

svgf=open("path.svg","w")

@printf(svgf,"<svg width=\"200\" height=\"200\" xmlns=\"http://www.w3.org/2000/svg\">\n")

#//  <path d="M10 10"/>

function printstartpath(x,y)
    @printf(svgf,"<path d=\"M %d %d ",x,y)
end

function printsegment(dx,dy)
    @printf(svgf," l %d %d",dx,dy)
    @printf("printsegment: %d %d\n",dx,dy)
end

function printendpath()
    @printf(svgf,"\" stroke=\"black\" fill=\"transparent\" stroke-width=\"0.5\"  \\>")
end

function drawpath()
    printstartpath(0,0)
    
    dx=1
    for x in 0:1:100
        printsegment(dx,dx*rand([1,-1]) )
    end
    
    printendpath()
end

drawpath()

@printf(svgf,"\n</svg>\n")
close(svgf)
