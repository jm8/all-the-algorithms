from random import randint

# returns true for prime numbers and some composite numbers
# returns false for most composite numbers
def fermat_primality_test(p, iterations = 3):
    if p < 2: return False
    if p == 3 or p == 2: return True

    for i in range(iterations):
        a = randint(2, p - 1)
        if pow(a, p - 1, p) != 1:
            return False
    return True