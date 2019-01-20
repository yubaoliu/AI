import cv2
import numpy as np

img = cv2.imread('../../Assets/Images/senn01.jpeg', 1)

gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
cv2.imshow('src', gray)
dst = cv2.equalizeHist(gray)
cv2.imshow('dst', dst)
cv2.waitKey(0)
