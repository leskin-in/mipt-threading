import random

random.seed()

array = []
for i in range(1, 2500001):
    array.append(i)

random.shuffle(array)

with open("data.txt", "w") as f:
    for ael in array:
        s = str(ael) + " "
        f.write(s)
    f.write("\n")

print("OK")
