import random
from sys import argv


if len(argv) > 2:
    print("Incorrect command line arguments")
    exit(-1)

try:
    range_limit = int(argv[1])
except:
    print("Incorrect command line arguments")
    exit(-1)



random.seed()

array = []
for i in range(1, range_limit + 1):
    array.append(i)

random.shuffle(array)

with open("data.txt", "w") as f:
    for ael in array:
        s = str(ael) + " "
        f.write(s)
    f.write("\n")

print("OK")
