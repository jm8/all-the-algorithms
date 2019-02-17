def bubble_sort(a):
    sorted = False
    while(not sorted):
        sorted = True
        for i in range(len(a) - 1):
            print(i)
            if a[i] > a[i+1]:
                temp = a[i]
                a[i] = a[i+1]
                a[i+1] = temp
                sorted = False
    return a