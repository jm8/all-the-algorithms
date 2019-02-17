from math import floor
def comb_sort(array, shrink = 1.3):
    gap = len(array)
    sorted = False
    while(not sorted):
        gap = max([floor(gap / shrink), 1])
        sorted = True if gap == 1 else False
        for i in range(len(array) - gap):
            if(array[i] > array[i+gap]):
                temp = array[i]
                array[i] = array[i+gap]
                array[i+gap] = temp
    return array
a = [98, 154, 1, 45, 14, 1576, 19, 10, 91, 2]
a = comb_sort(a)
print(a)