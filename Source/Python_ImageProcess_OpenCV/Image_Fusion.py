import cv2
import numpy as np
import random

img0 = cv2.imread('../../Assets/Images/flower0.jpeg', 1)
img1 = cv2.imread('../../Assets/Images/dog.jpeg', 1)

imgInfo = img0.shape
height = imgInfo[0]
width = imgInfo[1]

# ROI
roiH = int(height-50)
roiW = int(width-50)

img0ROI = img0[0:roiH, 0:roiW]
img1ROI = img1[0:roiH, 0:roiW]


dst = np.zeros((roiH, roiW, 3), np.uint8)
dst = cv2.addWeighted(img0ROI, 0.5, img1ROI, 0.5, 0)

cv2.imshow('img0', img0)
cv2.imshow('img1', img1)
cv2.imshow('dst', dst)
cv2.waitKey(0)