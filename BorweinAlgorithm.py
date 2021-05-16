import math
def BorweinQuadratic(iterations):
	a = math.sqrt(2)
	b = 0
	p = 2 + math.sqrt(2)
	for i in range(iterations):
		sqrta = math.sqrt(a)
		b = ((1 + b) * sqrta) / (a + b)
		a = (sqrta + (1 / sqrta)) / 2
		p = ((1 + a) * p * b) / (1 + b)
	return p
