# Feynman 1D path-integral Paths drawn with SVG Paths!
# See here for spec / tutorial: https://developer.mozilla.org/en/docs/Web/SVG/Tutorial/Paths

svgf=open("path.svg","w")

tSlices=100 # Also used as size of SVG canvas

@printf(svgf,"<svg width=\"%d\" height=\"%d\" xmlns=\"http://www.w3.org/2000/svg\">\n",tSlices,tSlices)

#//  <path d="M10 10"/>

function printstartpath(x,y)
    @printf(svgf,"<path d=\"M %d %d ",x,y)
end

function printsegment(dx,dy)
    @printf(svgf," l %d %d",dx,dy)
#    @printf("printsegment: %d %d\n",dx,dy)
end

function printendpath()
    colours=["#7fc97f","#beaed4","#fdc086","#ffff99","#386cb0","#f0027f"] # Colourbrewer 6 colours qualitative
    colour=rand(colours)

    @printf(svgf,"  \" stroke=\"%s\" fill=\"none\" stroke-width=\"0.5\"  \/>",colour)
end

function drawpath()
    printstartpath(tSlices/2.0,0)
    
    dt=1
    sx=0
    # Really need to compose this random vector so paths connect from [x0,t0]->[x1,t1]
    for t in 0:dt:tSlices
        dx=rand([1,-1]) # Forces either a (1,1) or (1,-1) move; zizzag space-time paths
        
        if (tSlices-t-1)<abs(sx) # Force a reconnect in the space-time paths
            dx=-sx/abs(sx) # Head in this direction only.
        end
        printsegment(dx,dt )
        sx+=dx # keep track of where we are, to reconnect space-time paths
    end
    
    printendpath()
end

# A few random paths...
Paths=50
for i in 1:Paths
    drawpath()
end

@printf(svgf,"\n</svg>\n")
close(svgf)
