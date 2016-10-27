
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
#    @printf("printsegment: %d %d\n",dx,dy)
end

function printendpath()
    @printf(svgf,"  \" stroke=\"black\" fill=\"transparent\" fill-opacity=\"0.0\" stroke-width=\"0.5\"  \/>")
end

function drawpath()
    printstartpath(0,100)
    
    dx=1
    # Really need to compose this random vector so paths connect from [x0,t0]->[x1,t1]
    for x in 0:dx:200
        printsegment(dx,dx*rand([1,-1]) )
    end
    
    printendpath()
end

# A few random paths...
for i in 1:10
    drawpath()
end

@printf(svgf,"\n</svg>\n")
close(svgf)
