# https://rosalind.info/problems/bins/

n = 5
m = 6
A = [10, 20, 30, 40, 50]
B = [40, 10, 35, 15, 40, 20]

def binary_search_loop(arr: list, n: int):
    low = 0
    high = len(arr) - 1
    while low <= high:
        mid = (low + high) // 2
        if arr[mid] == n:
            return mid + 1
        elif arr[mid] < n:
            low = mid + 1
        else:
            high = mid - 1
    return -1

res = [str(binary_search_loop(A, k)) for k in B]
print(" ".join(res))
