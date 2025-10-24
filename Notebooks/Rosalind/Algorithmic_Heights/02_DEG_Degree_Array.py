# https://rosalind.info/problems/deg/

input_str = """
    6 7
    1 2
    2 3
    6 3
    5 6
    2 5
    2 4
    4 1
"""

def parse_input(s: str) -> list[tuple[int, int]]:
    ret = []
    for line in input_str.splitlines():
        line = line.strip()
        if not line:
            continue
        tuple = [int(x) for x in line.split()]
        ret.append(tuple)
    return ret

edges = parse_input(input_str)

def problem(edges: list[tuple[int, int]]) -> list[int]:
    ret = [0 for _ in range(len(edges))]
    for [x, y] in edges:
        if not ret[x]:
            ret[x] = 1
        else:
            ret[x] += 1
        if not ret[y]:
            ret[y] = 1
        else:
            ret[y] += 1
    return ret[1:]

solution = problem(edges)
print(solution)
