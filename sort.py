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


def quick_sort_list(s):
    """
    Perform QUICK SORT procedure on a list of numbers
    :param s: List of number to be sorted
    :return:
    """
    def choose_pivot(a, lo, hi):
        pivot_idx = (hi - lo) // 2 + lo
        return pivot_idx

    def partition(a, lo, hi):
        pivot_idx = choose_pivot(a, lo, hi)
        pivot_val = a[pivot_idx]
        # Swap elements
        a[hi], a[pivot_idx] = a[pivot_idx], a[hi]

        store_idx = lo
        # Compare elements to pivot
        for i in range(lo, hi):
            if a[i] < pivot_val:
                a[i], a[store_idx] = a[store_idx], a[i]
                store_idx += 1
        # Move pivot to its final place
        a[hi], a[store_idx] = a[store_idx], a[hi]
        return store_idx

    def quick_sort(a, lo, hi):
        p = partition(a, lo, hi)
        if lo < p - 1:
            quick_sort(a, lo, p - 1)
        if p + 1 < hi:
            quick_sort(a, p + 1, hi)

    return quick_sort(s, 0, len(s)-1)


def topological_sort(adj_m):
    """
    Performs TOPOLOGICAL SORT on Directed Acyclic Graph (DAG)
    :param adj_m: Graph in the form of adjacency matrix.
                  Indexes should start from 0.
    :return: Sorted list of vertices
    """
    import numpy as np

    order = []

    used_mask = np.zeros((1, np.size(adj_m, axis=0)))
    sum_in = np.transpose(np.sum(adj_m, axis=0))

    while not np.all(used_mask):
        root_idxs = np.where((sum_in == 0) & np.logical_not(used_mask))[1]
        for idx in root_idxs:
            used_mask[0, idx] = 1
            sum_in -= adj_m[idx, :]
            order.append(idx)
    return order


def heap_sort(a):
    """
    Performs HEAP SORTING
    :param a: Correct HEAP (list) to be sorted in the form of vector
    :return: List with all elements from the heap sorted
    """
    def sift_down(h, idx):
        curr_idx = idx
        child_l = 2 * curr_idx + 1
        child_r = 2 * curr_idx + 2

        while True:
            if child_l >= len(h):
                break
            elif child_r >= len(h):
                if h[curr_idx] < h[child_l]:
                    h[curr_idx], h[child_l] = h[child_l], h[curr_idx]
                    curr_idx = child_l
                else:
                    break
            else:
                less_than_l = h[curr_idx] < h[child_l]
                less_than_r = h[curr_idx] < h[child_r]
                if not less_than_l and not less_than_r:
                    break
                elif less_than_l and less_than_r:
                    tmp_child = child_l if h[child_l] > h[child_r] else child_r
                elif less_than_l:
                    tmp_child = child_l
                else:
                    tmp_child = child_r

                h[curr_idx], h[tmp_child] = h[tmp_child], h[curr_idx]
                curr_idx = tmp_child

            child_l = 2 * curr_idx + 1
            child_r = 2 * curr_idx + 2

    b = a.copy()
    s = []

    while True:
        if len(b) == 1:
            s.append(b[0])
            break
        s.append(b[0])
        b[0] = b[-1]
        b = b[:-1]
        sift_down(b, 0)
    return reverse(s)