from math import floor
def stooge_sort(array, i = 0, j = - 1):
    j = len(array) - 1 if j == -1 else j
    if array[i] > array[j]:
        temp = array[i]
        array[i] = array[j]
        array[j] = temp
    if j - i + 1 > 2:
        t = floor((j - i + 1) / 3)
        stooge_sort(array, i, j - t)
        stooge_sort(array, i + t, j)
        stooge_sort(array, i, j - t)
    return array
print(stooge_sort([1254, 156, 14652, 124, 12, 63, 2116, 146, 254, 1245, 742]))