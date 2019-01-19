
import cv2
import numpy as np

img = cv2.imread('../../Assets/Images/flower-white.jpeg', 1)

imgInfo = img.shape
height = imgInfo[0]
width = imgInfo[1]

dst = np.zeros((height, width, 3), np.uint8)

for i in range(0, height):
    for j in range(0, width):
        (b, g, r) = img[i, j]
        dst[i, j] = (255 - b, 255-g, 255-r)

cv2.imshow('src', img)
cv2.imshow('dst', dst)
cv2.waitKey(0)