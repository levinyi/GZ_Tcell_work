import sys



with open() as f:
    for line in f:
        line = line.rstrip("\n")
        if line.startswith("#"):
            print(line)
            continue
        a = line.split("/")
