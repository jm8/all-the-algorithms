from random import shuffle
def is_sorted(array):
    for i in range(len(array) - 1):
        if array[i] > array[i+1]:
            return False
    return True
def bogo_sort(array):
    i = 0
    while not is_sorted(array):
        shuffle(array)
        i += 1
    print(i)
    return array