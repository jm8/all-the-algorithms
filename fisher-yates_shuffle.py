from random import randint

# improved version of Fisher-Yates shuffle by Durstenfeld
def shuffle(array):
    for i in range(len(array) - 1, 0, -1):
        j = randint(0, i)
        # swap array elements at i and j
        t = array[i]
        array[i] = array[j]
        array[j] = t
    return array
