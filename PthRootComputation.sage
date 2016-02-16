import timeit
# Field Characteristic
p = 5
F, x = GF(p)['x'].objgen()
Q = GF(p)

# Fixed Element to compute roots
_const_element = x**203 + x**202 + 4*x**201 + x**199 + x**197 \
+ 4*x**195 + x**194 + x**193 + x**154 + 3*x**152 + 3*x**151 \
+ 4*x**150 + 3*x**149 + x**147 + 4*x**146 + x**145 + x**102 \
+ 4*x**100 + x**99 + 2*x**98 + 3*x**97 + 4*x**52 + x**51 + x**50 \
+ 2*x**49 + x**4 + x**2 + x + 2

#Irreducible binomial over F5
_irr_bin = x**256 + 2
#Irreducible trinomial over F5
_irr_tri = x**271 + 2*x**121 + 1
#Irreducible friendly tetranomial over F5
_irr_tet = x**271 + x**196 + 4*x**121 + 1
#irreducible friendly pentanomial over F5
_irr_pen = x**271 + 2*x**241 + 2*x**211 + 4*x**181 + 2
#Irreducible friendly 40-nomial
_irr_40 = x**271 + x**266 + x**261 + x**256 + x**251 + x**246\
 + x**241 + x**236 + x**231 + x**226 + x**221 + x**216 + x**211\
 + x**206 + x**201 + x**196 + x**191 + x**186 + x**181 + x**176\
 + x**171 + x**166 + x**161 + x**156 + x**151 + x**146 + x**141\
 + x**136 + x**131 + x**126 + x**121 + x**116 + x**111 + x**106\
 + x**101 + x**96 + x**91 + 3*x**86 + 4*x**81 + 3
#Irreducible equally spaced heptanomial over F5
_irr_eql_hept = x**294 + x**245 + x**196 + x**147 + x**98 + x**49 + 1
_hept_consts = [x**549, x**69, x**275, x**138]
#Irreducible equally spaced 17-nomial over F5
_irr_eql_big = x**272 + x**255 + x**238 + x**221 \
+ x**204 + x**187 + x**170 + x**153 + x**136 + x**119 \
+ x**102 + x**85 + x**68 + x**51 + x**34 + x**17 + 1
_big_consts = [x**58, x**116, x**463, x**232]

def printHw(consts):
	for i  in range(0, p-1):
		print "x"+str(i+1)+"="+str(consts[i])+", "

def breakElement(element):
	E = [0]*p
	for exp, coef in enumerate(element.coefficients(sparse=False)):
		i = exp%p
		E[i] = E[i] + coef*x**((exp-i)/p)
	return E

def compute_friendly_constants(poly, hw=False):
	const = [0]*(p-1)
	j = [0]*(p-1)
	t = [0]*(p-1)
	aux = poly.coefficients(sparse=False)
	a0 = Q(aux[0])
	r = (len(aux) - 1)%p
	f = sum(c*x**i for i, c in enumerate(aux[r::p]))
	u = int(Q(1)/(p-r))
	for i in range(1,p):
		j[i-1] = (u*i)%p
		t[i-1] = (u*i - j[i-1])/p
	for i in range(1,p):
		a0jip = Q(1)/( (-a0)**j[i-1])
		const[i-1] =	a0jip*x**( i*t[r-1] - r*t[i-1] + i )*f**j[i-1]
	if hw :
		printHw(const)
	return const

def pthRoot(C, consts):
	return C[0] + sum( consts[i-1]*C[i] for i in range(1,p))

def pthRootWithReduction(C, consts, poly):
	return (C[0] + sum( consts[i-1]*C[i] for i in range(1,p)))

def benchmark(C, consts):
	t = timeit.Timer(lambda: pthRoot(C, consts)).timeit(number = 100000)/100000
	return t

def benchmarkWithReduction(C, consts, poly):
	t = timeit.Timer(lambda: pthRootWithReduction(C, consts, poly)).timeit(number =100000)/100000
	return t

def test_friendly(poly, msg="Testing friendly"):
	print "***\n" + msg
	C = breakElement(_const_element)
	consts = compute_friendly_constants(poly)
	test = (pthRoot(C, consts)**p)%poly
	assert(_const_element/test == 1)
	print benchmark(C, consts)

def test_equally(poly, consts, msg="Testing equally"):
	print "***\n" + msg
	C = breakElement(_const_element)
	test = (pthRootWithReduction(C, consts, poly)**p)%poly
	assert(_const_element/test == 1)
	print benchmarkWithReduction(C, consts, poly)

test_friendly(_irr_bin, "Irreducible Binomial")
test_friendly(_irr_tri, "Irreducible Trinomial")
test_friendly(_irr_tet, "Irreducible Tetranomial")
test_friendly(_irr_pen, "Irreducible Pentanomial")
test_friendly(_irr_40, "Irreducible 40-nomial")
test_equally(_irr_eql_hept, _hept_consts, "Irreducible heptanomial")
test_equally(_irr_eql_big, _big_consts, "Irreducible 17-nomial")
