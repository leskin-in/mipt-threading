with open("data.txt") as f:
    data = f.read()

lines = data.split("\n")

for i in range(1, len(lines)):
    words = lines[i].split(" ")
    for j in range(len(words) - 2):
        if(int(words[j]) > int(words[j + 1])):
            print("Found error in line", i)
            print(words[j], ">", words[j+1])
            break

print("Finished")
