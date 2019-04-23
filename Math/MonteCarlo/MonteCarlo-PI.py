import random

total = 5000000
in_count = 0

for i in range(total):
    x=random.random() #0-1
    y=random.random() #0-1

    dis=(x**2+y**2)**0.5

    if dis <= 1: # inside the circle
        in_count+=1

print('Pi= ', 4*in_count/total)
