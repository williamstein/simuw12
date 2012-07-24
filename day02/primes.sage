colors = ['red', 'blue', 'green', 'purple', 'cyan']

for p in primes(200):
    print r"\tikz \node[draw,circle]{\color{" + colors[p%len(colors)] + "}%s};\n"%p
