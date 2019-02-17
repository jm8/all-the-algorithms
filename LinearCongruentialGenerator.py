class LinearCongruentialGenerator:
    def __init__(self, a, c, m, seed = 0):
        self.a = a
        self.c = c
        self.m = m
        self.current = seed
    def new(self):
        self.current *= self.a
        self.current += self.c
        self.current %= self.m
        return self.current