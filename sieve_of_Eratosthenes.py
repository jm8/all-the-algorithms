import math

# returns a list of all numbers less than or equal to n
def sieve(n):
    l = [True] * (n + 1)
    l[0] = False
    l[1] = False
    for i in range(2, math.ceil(math.sqrt(n))):
        if l[i]:
            for j in range(i * i, n + 1, i):
                l[j] = False
    output = []
    for i in range(0, len(l)):
        if l[i]:
            output.append(i)
    return output