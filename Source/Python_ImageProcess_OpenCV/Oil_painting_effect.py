import cv2
import numpy as np

img = cv2.imread('../../Assets/Images/dog.jpeg', 1)

imgInfo = img.shape
height = imgInfo[0]
width = imgInfo[1]

cv2.imshow('img', img)

gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)

dst = np.zeros((height, width, 3), np.uint8)
for i in range(4, height-4):
    for j in range(4, width-4):
        array1 = np.zeros(8, np.uint8) # 8个灰度等级的像素个数
        for m in range(-4, 4):
            for n in range( -4, 4):
                p1 = int(gray[i+m, j+n]/32)
                array1[p1] = array1[p1]+1

        #计算哪个段所含像素最多
        currentMax = array1[0]
        l = 0 #当前处于哪个块中
        for k in range(0, 8):
            if currentMax < array1[k]:
                currentMax = array1[k]
                l = k

       #Calculate the mean value

        count = 0
        newB = 0
        newG = 0
        newR = 0
        for m in range(-4, 4):
            for n in range(-4, 4):
                if gray[i+m, j+n] >= (l*32) and gray[i+m, j+n] <= ( (l+1)*32):
                    (b, g, r) = img[i+m, j+n]
                    newB = newB + b
                    newG = newG + g
                    newR = newR + r
                    count = count + 1

        dst[i, j] = ( int(newB/count), int(newG/count), int(newR/count))

cv2.imshow('dst',dst)
cv2.waitKey(0)
