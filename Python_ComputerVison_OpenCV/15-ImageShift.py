import cv2
import numpy as np

img = cv2.imread('../../Assets/Images/lena.jpg', 1)

cv2.imshow("src",img)
imgInfo = img.shape

height = imgInfo[0]
width = imgInfo[1]


matShift = np.float32([[1, 0, 100], [0, 1, 100]]) # 2 3

dst = cv2.warpAffine(img, matShift, (height, width)) # p1: original data; p2: shift mat; p3: image info

cv2.imshow('dst', dst)

cv2.waitKey(0)