def partition(array, low, high):
    pivot = array[high]
    i = low
    for j in range(low, high):
        if array[j] < pivot:
            temp = array[j]
            array[j] = array[i]
            array[i] = temp
            i += 1
    temp = array[i]
    array[i] = array[high]
    array[high] = temp
    return i

def quick_sort(array, low = 0, high = -1):
    high = len(array) - 1 if high == -1 else high
    if low < high:
        p = partition(array, low, high)
        quick_sort(array, low, p-1)
        quick_sort(array, p+1, high)
    return a
    
a = [98, 154, 1, 45, 14, 1576, 19, 10, 91, 2]
print(quick_sort(a))