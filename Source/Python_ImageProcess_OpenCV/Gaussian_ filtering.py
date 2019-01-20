import cv2

img = cv2.imread('../../Assets/Images/LenaWithNoise.png', 1)
cv2.imshow('src', img)


dst = cv2.GaussianBlur(img, (5, 5), 1, 5)


cv2.imshow('dst', dst)
cv2.waitKey(0)