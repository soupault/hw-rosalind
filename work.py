from sort import *
from misc import *

file_in = "input_small.txt"
# file_in = "input_large.txt"
file_out = "result.txt"

if __name__ == "__main__":
    # ---------------------------- Get data ----------------------------------
    with open(file_in) as f:
        lines = [e.strip() for e in f.readlines()]

    # Input is an INTEGER
    # n = int(lines[0])
    # Input is an ARRAY
    # A = [int(e) for e in lines[1].split()]
    # Input is a SET OF ARRAYS
    # data = []
    # for l in lines[1:]:
    #     data.append([int(e) for e in l.split()])

    n = int(lines[0])
    A = [int(e) for e in lines[1].split()]

    # print(n, "\n", A, "\n")
    # quit()
    # -------------------------- Main routine --------------------------------
    result = []

    # print(result)
    # quit()
    # --------------------------- Write data ---------------------------------
    with open(file_out, 'w') as f:
        # Result is a NUMBER
        # f.writelines(str(result))
        # Result is a LIST
        # f.writelines([str(e) + " " for e in result])
        # Result is a SET OF LISTS
        # for line in result:
        #     f.writelines([str(e) + " " for e in line] + ["\n"])
