from sort import *
from misc import *

# file_in = "input_small.txt"
file_in = "input_large.txt"
file_out = "result.txt"

if __name__ == "__main__":
    # ---------------------------- Get data ----------------------------------
    with open(file_in) as f:
        lines = [e.strip() for e in f.readlines()]

    # --- Input is an INTEGER
    # n = int(lines[0])
    # --- Input is an ARRAY
    # A = [int(e) for e in lines[1].split()]
    # --- Input is a SET OF ARRAYS
    # data = []
    # for l in lines[1:]:
    #     data.append([int(e) for e in l.split()])
    # --- Input is a STRING
    # s = lines[0]
    # --- Input is an ARRAY OF STRINGS
    # s, t = [l for l in lines]

    # print(s, "\n", t, "\n")
    # quit()
    # -------------------------- Main routine --------------------------------
    result = []

    # print(result)
    # quit()
    # --------------------------- Write data ---------------------------------
    with open(file_out, 'w') as f:
        pass
        # --- Result is a NUMBER
        # f.writelines([str(result) + "\n"])
        # --- Result is a LIST
        # f.writelines([str(e) + " " for e in result])
        # --- Result is a SET OF LISTS
        # for line in result:
        #     f.writelines([str(e) + " " for e in line] + ["\n"])
        # --- Result is a STRING
        # f.writelines(result)
        # --- Result is a SET OF STRINGS
        # f.writelines([line + "\n" for line in result])
