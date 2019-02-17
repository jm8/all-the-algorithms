# Blum Blum Shub RNG
class BlumBlumShub:
    def __init__(self, seed, p, q):
        self.p = p
        self.q = q
        self.current = seed
    def new(self):
        self.current = (self.current**2) % (self.p * self.q)
        return self.current
