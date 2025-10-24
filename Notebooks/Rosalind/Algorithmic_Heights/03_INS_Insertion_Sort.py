# https://rosalind.info/problems/ins/

def insertion_sort(arr: list) -> int:
    count = 0
    for i in range(1, len(arr)):
        key = arr[i]
        j = i - 1
        while j >= 0 and key < arr[j]:
            arr[j+1] = arr[j]
            j -= 1
            count += 1
        arr[j+1] = key
    return count

arr = [6, 10, 4, 5, 1, 2]
print(insertion_sort(arr))
print(arr)
