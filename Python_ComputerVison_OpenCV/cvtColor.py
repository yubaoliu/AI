
import cv2

img = cv2.imread('../../Assets/Images/flower-white.jpeg', 1)

print(img.shape)

dst = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)

cv2.imshow('dst', dst)
cv2.waitKey(0)