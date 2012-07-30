"""
Functions related to the congruent number problem for elliptic curves.

AUTHOR:
    -- 2005-10-17: William Stein
    -- 2006-03-12: William Stein; updated for SAGE 1.1.0

EXAMPLES:

    sage: attach "cong.sage"
    sage: cong_number_sets(101)
    ([], [])
    sage: is_conj_cong_number(101)
    True
    sage: is_cong_number(101)
    True
    sage: print conj_cong_number_list(30)
    [5, 6, 7, 13, 14, 15, 20, 21, 22, 23, 24, 28, 29, 30]
    sage: cong_number_curve(101)
    Elliptic Curve defined by y^2  = x^3 - 10201*x over Rational Field
    sage: cong_curve_gens(101)
    Entering qsieve::search: y^2 = 1x^3 + 0x^2 + -10201x^1 + 0
     ...
    sage: P = _[0]
    sage: cong_triangle(P)
    (-2332125199241651236502079933/936040998152650,
     267980280100/44538033219,
     6996375597724953709506260201/2808122994457950)    
    sage: cong_triangle(2*P)
    (194078525650704658678027874375600681957178024918822355502897426822473572163657494987784312896632352861060981280801/51931304811041127502104551232937598672828704386020947660827418999486010525512364380759400727151544100,
    499788470246347274586979219879071925647657369364415186761776176129573731800/9246622194528447038907717193458873919847209114593519978984797638591935971364011099399,
    -194078525650704658678027874375600681957178024939120541557027444504370638033509208129865614081887258432093181280801/51931304811041127502104551232937598672828704386020947660827418999486010525512364380759400727151544100)
    """

def cong_number_sets(n):
    """
    Given a positive integer n, returns the two sets appearing in the
    conjectural criterion for when a number is congruent.
    """
    n = ZZ(n)
    n = ZZ(prod([p for p, e in n.factor() if e%2 == 1]))
    
    if n % 2 == 0:   # even case
        E = []       # with c even
        O = []       # with c odd
        a_bound = (sqrt(n/8) + 1).floor()
        c_bound = (sqrt(n)/4 + 1).floor()
        half_n = n//2
        for c in range(-c_bound, c_bound+1):
            c_square = c^2
            for a in range(-a_bound, a_bound+1):
                a_square = a^2
                z = half_n - 4*a_square - 8*c_square
                try:
                    b = z.sqrt(extend=False)
                    if b.denominator() == 1:
                        b = b.numerator()
                        assert 4*a^2 + b^2 + 8*c^2 == n/2
                        if c % 2 == 0:
                            E.append((a,b,c))
                        else:
                            O.append((a,b,c))
                except ValueError:
                    pass
    else:
        E = []       # with c even
        O = []       # with c odd
        a_bound = floor(sqrt(n/2)+1)
        c_bound = floor(sqrt(n/8) + 1)
        for c in range(-c_bound, c_bound+1):
            c_square = c^2
            for a in range(-a_bound, a_bound+1):
                a_square = a^2
                z = n - 2*a_square - 8*c_square
                try:
                    b = z.sqrt(extend=False)
                    if b.denominator() == 1:
                        b = b.numerator()
                        assert 2*a^2 + b^2 + 8*c^2 == n
                        if c % 2 == 0:
                            E.append((a,b,c))
                        else:
                            O.append((a,b,c))
                except ValueError:
                    pass
    return E, O

def is_conj_cong_number(n):
    """
    Returns True if n is conjecturally a congruent number, according
    to the Birch and Swinnerton-Dyer conjecture.
    """
    E, O = cong_number_sets(n)
    return len(E) == len(O)

def is_cong_number(n, verbose=True):
    """
    Returns True if n is provably a congruent number.  This computes
    the rank of the corresponding elliptic curve.
    """
    return len(cong_curve_gens(n, verbose=verbose)) > 0

def conj_cong_number_list(bound):
    """
    Returns the conjectural congruent numbers up to a given bound.
    """
    return [n for n in range(1,bound+1) if is_conj_cong_number(n)]

def cong_number_curve(n):
    """
    Returns the elliptic curve corresponding to the integer n.  This
    curve has positive rank if and only if n is a congruent number.
    """
    return EllipticCurve([-n^2,0])

def cong_curve_gens(n, verbose=True):
    """
    Returns generators for the Mordell-Weil group of the elliptic
    curve corresponding to n.
    """
    E = cong_number_curve(n)
    return E.gens(verbose=verbose, proof=False, descent_second_limit=25)

def cong_triangle(P):
    """
    Returns the triangle corresponding to the point P on a congruent
    number elliptic curve.  The point P must have nonzero y coordinate.
    """
    m = P.curve().a_invariants()[3]
    n = (-m).sqrt()
    (x, y) = (P[0],P[1])
    if y == 0:
        raise ArithmeticError, "point does not define a rational right triangle"
    T = ((n^2-x^2)/y, -2*x*n/y, (n^2+x^2)/y)
    assert (T[0]^2 + T[1]^2 == T[2]^2), "bug in program."
    return T

def point_from_triangle(a,b,c):
    """
    Return the point on $y^2=x^3-n^2x$ corresponding to the
    right triangle with side lengths $a,b,c$.
    """
    n = a*b/2
    E = EllipticCurve([-n^2, 0])
    return E([-n*b/(a+c), 2*n^2/(a+c)])

def rational_n(x):
    """
    Returns a rational number.  For any rational Q, there exists an
    integer x such that Q = rational_n(x).  Negative x values result
    in negative Q.
    """
    if x == 0:
      return QQ(0)
    if x < 0:
      s = -1
      x = -x
    else:
      s = 1
    r = floor((sqrt(1+8*(x-1))-1)/2)
    t = (r^2 + r)/2
    n = x - t
    d = r - n + 2
    return s*n/d


def rational_triangles(bnd):
    for i in range(2,bnd):
        t = rational_n(i)
        x = (t^2-1)/(t^2+1)
        y = (t^3)/(t^2+1)
        if not ( x == -1 or x == 0):
            print "%s,%s"%(x,y)


def L(n, prec):
    """
    Compute $L(E_n,1)$ using prec $a_n$.

    EXAMPLES:
        sage: L(1,10)
        0.65551389243851388
        sage: L(34,272)
        -0.00000021288041361683671        
    """
    E = EllipticCurve([-n^2,0])
    v = E.anlist(prec)
    N = E.conductor()
    sqrtN = sqrt(N)
    lval = 2*sum(v[n]/n*exp(-2*pi*n/sqrtN) for n in range(1,prec))
    return lval


def theta(prec, c=1, base_ring=ZZ):
    """
    Return the theta series Theta(q^c) to precision O(q^prec) over the
    given base_ring.
    """
    R.<q> = PowerSeriesRing(base_ring)
    v = [0]*prec
    j = int(sqrt(prec))
    c = int(c)
    for n in range(1,j):
        m = n*n
        if m*c >= prec:
            break
        v[m*c] = 2
    v[0] = 1
    return R(v, prec)

def tunnell_forms(prec=30, base_ring=ZZ):
    """
    Return the two q-expansions appearing in Tunnell's theorem.
    """
    t=cputime()
    T = theta(prec, 1, base_ring)
    T2 = theta(prec, 2, base_ring)
    T4 = theta(prec, 4, base_ring)
    T8 = theta(prec, 8, base_ring)
    T32 = theta(prec, 32, base_ring)
    print "time to construct series", cputime(t)
    t = cputime()
    f = T*(2*T32 - T8)
    f1 = f*T2
    f2 = f*T4
    print "time to multiply out", cputime(t)
    return f1, f2

def congruent_numbers(nmax):
    """
    Return all congruent numbers (according to Tunnell's
    criterion, so conjectural) that are less than nmax.
    """
    f1, f2 = tunnell_forms(nmax, base_ring=ZZ)
    v1 = [n for n in range(nmax) if n%2==1 and f1[n]==0 and is_squarefree(n)]
    v2 = [n for n in range(nmax) if n%2==0 and f2[int(n/2)]==0 and is_squarefree(n)]
    v = v1 + v2
    # now add in multiples by squares.
    s = []
    for a in v:
        for j in range(2,int(sqrt(nmax/a))+1):
            n = a*j*j
            if n < nmax:
                s.append(a*j*j)
    t = v + s
    t.sort()
    return t
    
    

table1000 = [
5, 6, 7, 13, 14, 15, 20, 21, 22, 23, 24, 28, 29, 30, 31, 34, 37, 38,
39, 41, 45, 46, 47, 52, 53, 54, 55, 56, 60, 61, 62, 63, 65, 69, 70,
71, 77, 78, 79, 80, 84, 85, 86, 87, 88, 92, 93, 94, 95, 96, 101, 102,
103, 109, 110, 111, 112, 116, 117, 118, 119, 120, 124, 125, 126, 127,
133, 134, 135, 136, 137, 138, 141, 142, 143, 145, 148, 149, 150, 151,
152, 154, 156, 157, 158, 159, 161, 164, 165, 166, 167, 173, 174, 175,
180, 181, 182, 183, 184, 188, 189, 190, 191, 194, 197, 198, 199, 205,
206, 207, 208, 210, 212, 213, 214, 215, 216, 219, 220, 221, 222, 223,
224, 226, 229, 230, 231, 237, 238, 239, 240, 244, 245, 246, 247, 248,
252, 253, 254, 255, 257, 260, 261, 262, 263, 265, 269, 270, 271, 276,
277, 278, 279, 280, 284, 285, 286, 287, 291, 293, 294, 295, 299, 301,
302, 303, 306, 308, 309, 310, 311, 312, 313, 316, 317, 318, 319, 320,
323, 325, 326, 327, 330, 333, 334, 335, 336, 340, 341, 342, 343, 344,
348, 349, 350, 351, 352, 353, 357, 358, 359, 365, 366, 367, 368, 369,
371, 372, 373, 374, 375, 376, 380, 381, 382, 383, 384, 386, 389, 390,
391, 395, 397, 398, 399, 404, 405, 406, 407, 408, 410, 412, 413, 414,
415, 421, 422, 423, 426, 429, 430, 431, 434, 436, 437, 438, 439, 440,
442, 444, 445, 446, 447, 448, 453, 454, 455, 457, 461, 462, 463, 464,
465, 468, 469, 470, 471, 472, 476, 477, 478, 479, 480, 485, 486, 487,
493, 494, 495, 496, 500, 501, 502, 503, 504, 505, 508, 509, 510, 511,
514, 517, 518, 519, 525, 526, 527, 532, 533, 534, 535, 536, 540, 541,
542, 543, 544, 546, 548, 549, 550, 551, 552, 557, 558, 559, 561, 564,
565, 566, 567, 568, 572, 573, 574, 575, 580, 581, 582, 583, 585, 589,
590, 591, 592, 596, 597, 598, 599, 600, 602, 604, 605, 606, 607, 608,
609, 613, 614, 615, 616, 621, 622, 623, 624, 628, 629, 630, 631, 632,
636, 637, 638, 639, 644, 645, 646, 647, 651, 653, 654, 655, 656, 658,
660, 661, 662, 663, 664, 668, 669, 670, 671, 674, 677, 678, 679, 685,
686, 687, 689, 692, 693, 694, 695, 696, 700, 701, 702, 703, 709, 710,
711, 717, 718, 719, 720, 721, 723, 724, 725, 726, 727, 728, 731, 732,
733, 734, 735, 736, 741, 742, 743, 749, 750, 751, 752, 756, 757, 758,
759, 760, 761, 764, 765, 766, 767, 773, 774, 775, 776, 777, 781, 782,
783, 788, 789, 790, 791, 792, 793, 796, 797, 798, 799, 805, 806, 807,
813, 814, 815, 820, 821, 822, 823, 824, 828, 829, 830, 831, 832, 837,
838, 839, 840, 845, 846, 847, 848, 850, 852, 853, 854, 855, 856, 860,
861, 862, 863, 864, 866, 869, 870, 871, 876, 877, 878, 879, 880, 884,
885, 886, 887, 888, 889, 890, 892, 893, 894, 895, 896, 901, 902, 903,
904, 905, 909, 910, 911, 915, 916, 917, 918, 919, 920, 924, 925, 926,
927, 933, 934, 935, 941, 942, 943, 948, 949, 950, 951, 952, 956, 957,
958, 959, 960, 965, 966, 967, 973, 974, 975, 976, 980, 981, 982, 983,
984, 985, 987, 988, 989, 990, 991, 992, 995, 997, 998, 999]





def random_elliptic_curve(p):
    """
    Construct and return a random elliptic curver over the finite
    field of order p.
    """
    p = ZZ(p)
    if not is_prime(p):
        raise ValueError, "p (=%s) must be a prime integer."%p
    F = FiniteField(p)
    while True:
        try:
            return EllipticCurve(F, [F.random_element(), F.random_element()])
        except ArithmeticError:
            pass
    return E
    
