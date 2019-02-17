# Lagged Fibonacci Generator
class LaggedFibonacci:
    def __init__(self, seed, j, k, m):
        self.j =j
        self.k = k
        self.m = m
        self.current = seed

    def fibonacci(self, x):
        if x < 2: # fibonacci sequence starts at 0
            return 1
        if x < len(self.LUT):
            return self.LUT[x]
        return self.fibonacci(x- 1) + self.fibonacci(x - 2)

    def fibonacci_preload(self, amount):
        f1 = 1 
        f2 = 1
        self.LUT = [1, 1]
        for i in range(amount):
            temp = f2
            f2 += f1
            f1 = temp
            self.LUT.append(f2)

    def new(self):
        self.current += 1
        return (self.fibonacci(self.current - self.j) + self.fibonacci(self.current - self.k)) % self.m
r = LaggedFibonacci(12, 7, 10, 100)