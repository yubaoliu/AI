import cv2
import numpy as np

img = cv2.imread('../../Assets/Images/lena.jpg', 1)

cv2.imshow("src",img)
imgInfo = img.shape

height = imgInfo[0]
width = imgInfo[1]


matShift = np.float32([[0.5, 0, 0], [0, 0.5, 0]]) # 2 3

dst = cv2.warpAffine(img, matShift, (int(height/2), int(width/2) )) # p1: original data; p2: shift mat; p3: image info

cv2.imshow('dst', dst)

cv2.waitKey(0)