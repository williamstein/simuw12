def random_function():
    import random
    f = 1
    x = var('x')
    for i in range(5):
        f += random.choice([-2,-1,1,2])*random.choice([sin(x), cos(x)]).subs(x=random.choice([-3,-2,-1,1,2,3])*x+random.randint(0,2))
    return f.expand()

def exercise(n):
    xmin=-1; xmax=2*pi
    while True:
        f = random_function()
        g = f.derivative()
        pg = g.plot(xmin, xmax, color='lightgreen')
        if max(abs(pg.ymax()), abs(pg.ymin())) < xmax:
            break
    s = "Below is a plot of $$f(x)=%s.$$  The {\em \color{red}derivative} of $f(x)$ is the function whose value at $x$ is the {\\em slope} of the graph of $f$ at $x$.  Plot the derivative of $f(x)$ by sketching the tangent line to the graph at maybe 10 points, and at each point, plot the slope of that line, then connect your points (it's a good idea to include all points at which the derivative is 0).  There is enough space vertically to fit the derivative.  {\em After} you finish carefully plotting the derivative, enter $f(x)$ into Sage, and plot {\color{blue}\\verb|f.derivative()|} to check your work.\n"%latex(f)
    s += "\\begin{center}\\includegraphics{functions/%s.pdf}\\end{center}\\newpage\n\n"%n
    pf = f.plot(xmin,xmax,thickness=3,axes=True,frame=True,gridlines='minor',aspect_ratio=1)
    pf.ymin(min(pf.ymin(),pg.ymin()))
    pf.ymax(max(pf.ymax(),pg.ymax()))
    pf.save('functions/%s.pdf'%n)
    open('functions/%s.tex'%n,'w').write(s)
    #(pf+pg).save('functions/%s-back.pdf'%n)
