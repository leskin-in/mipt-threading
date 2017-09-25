with open("data.txt") as f:
    data = f.read()

lines = data.split("\n")

for i in range(1, len(lines)):
    words = lines[i].split(" ")
    for j in range(len(words) - 1):
        if(words[i] > words[i + 1]):
            print("Found error in line", i)
            break

print("Finished")
