
import cv2

img0 = cv2.imread('../../Assets/Images/flower-white.jpeg', 1)
img1 = cv2.imread('../../Assets/Images/flower-white.jpeg', 0)

print(img0.shape)
print(img1.shape)

cv2.imshow('src', img1)
cv2.waitKey(0)

