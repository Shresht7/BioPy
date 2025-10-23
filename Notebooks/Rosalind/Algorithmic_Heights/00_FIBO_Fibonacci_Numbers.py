# https://rosalind.info/problems/fibo/

def fibonacci_recursion(n: int) -> int:
    """Returns the nth number in the Fibonacci sequence"""
    if n == 0 or n == 1:
        return n
    return fibonacci_recursion(n - 1) + fibonacci_recursion(n - 2)

print(fibonacci_recursion(6))

def fibonacci_for(n: int) -> int:
    """Returns the nth number in the Fibonacci sequence"""
    if n == 0 or n == 1:
        return n
    a, b = 0, 1
    for _ in range(n - 1):
        a, b = b, a + b
    return b

print(fibonacci_for(6))
