from itertools import product, combinations, imap
# Field characteristic.
p = 3
# Hamming weight.
k = 4

#Generator of combination of friendly exponents
def friendly_exp(m):
    max_d = (m-1)/p + 1
    for dist in combinations(range(1, max_d), k-2):
        yield [m] + [m - d*p for d in dist] + [0]
#Generator of combination of improved friendly exponents
def improved_exp(m):
    max_d = ((m-1)/p + 1) / (k-2)
    for d in range(1, max_d):
        yield [m] + [m-i*d*p for i in range(1, k-1)] + [0]

#Checks if polynomial is irreducible
x = GF(p)['x'].gen()
def has_irr(friendly):
    for coefs in product(range(1, p), repeat=k):
        poly = sum(coef * x**exp for coef, exp \
        in zip(coefs, friendly))
        if poly.is_irreducible():
            return poly
            
#Program Loop
for m in range(10, 550):
    if m in Primes():
        has_friendly = any(imap(has_irr, friendly_exp(m)))
        has_improved = any(imap(has_irr, improved_exp(m)))
        print m, has_friendly, has_improved
