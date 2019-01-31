import cv2
import numpy as np

img = cv2.imread('../../Assets/Images/lena.jpg', 1)

imgInfo = img.shape
height = imgInfo[0]
width = imgInfo[1]

gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
cv2.imshow('src', gray)

count = np.zeros(256, np.float)

for i in range(0, height):
    for j in range(0, width):
        pixel = gray[i, j]
        index = int(pixel)
        count[index] = count[index]+1

for i in range(0, 256):
    count[i] = count[i]/(height*width)

#计算累积概率
sum1 = float(0)
for i in range(0, 256):
    sum1 = sum1 + count[i]
    count[i] = sum1

#print(count)

#计算映射表
map1 = np.zeros(256, np.uint16)
for i in range(0, 256):
    map1[i] = np.uint16(count[i]*255)

for i in range(0, height):
    for j in range(0, width):
        pixel = gray[i,j]
        gray[i, j] = map1[pixel]

cv2.imshow('dst', gray)
cv2.waitKey(0)