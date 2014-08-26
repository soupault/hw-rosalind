__author__ = 'egor'


def bin_search(lst, elem, offset=0):
    if len(lst) == 0:
        return None
    mid_idx = len(lst) // 2
    curr_elem = lst[mid_idx]
    if elem == curr_elem:
        return offset + mid_idx + 1
    elif elem < curr_elem:
        return bin_search(lst[:mid_idx], elem, offset)
    elif elem > curr_elem:
        return bin_search(lst[mid_idx + 1:], elem, offset + 1 + mid_idx)


def insertion_sort(A):
    n = len(A)
    for i in range(1, n):
        k = i
        while k > 0 and A[k] < A[k - 1]:
            A[k], A[k - 1] = A[k - 1], A[k]
            k -= 1
    return A


def merge_sort(A):
    if len(A) == 1:
        return A
    elif len(A) > 2:
        mid_idx = len(A) // 2
        return merge_sorted_arrays(merge_sort(A[:mid_idx]),
                                   merge_sort(A[mid_idx:]))
    else:
        return insertion_sort(A)


def merge_sorted_arrays(array_a, array_b):
    idx_a, idx_b = 0, 0
    result = []

    while True:
        if array_a[idx_a] < array_b[idx_b]:
            result.append(array_a[idx_a])
            idx_a += 1
        else:
            result.append(array_b[idx_b])
            idx_b += 1

        if idx_a >= len(array_a):
            result.extend(array_b[idx_b:])
            break
        elif idx_b >= len(array_b):
            result.extend(array_a[idx_a:])
            break
    return result