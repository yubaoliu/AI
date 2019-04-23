# Calculate the square of interest area in a logo image using Mont carlo
from PIL import Image
import random

img = Image.open('../../Assets/Images/apple-logo.png')

count = 200000
in_count = 0

for i in range(count):
    x = random.randint(0, img.width-1)
    y = random.randint(0, img.height-1)

    color = img.getpixel((x, y))
#    print(color)
    if color == 1:
        in_count += 1


print('image square= ', img.width*img.height)
# Estimate the count number in interested area
print('Estimated Square: ', int(img.width*img.height*in_count/count))

# Real square
real_count = 0
for x in range(img.width):
    for y in range(img.height):
        if 1 == img.getpixel((x,y)):
            real_count += 1

print('Real Square: ', real_count)
