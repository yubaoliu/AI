import cv2
import numpy as np

img = cv2.imread('../Assets/Images/flower2.jpeg', 1)
for i in range(200, 300):
    img[i, 200] = (255, 255, 255)
    img[i, 200-1] = (255, 255, 255)
    img[i, 200+1] = (255, 255, 255)

for i in range(150, 250):
    img[250, i] = (255, 255, 255)
    img[250+1, i] = (255, 255, 255)
    img[250-1, i] = (255, 255, 255)

#cv2.imwrite('damaged.jpg', img)
cv2.imshow('damaged image', img)

imgInfo = img.shape
height = imgInfo[0]
width = imgInfo[1]

paint = np.zeros((height, width, 1), np.uint8)

for i in range(200, 300):
    paint[i, 200] = 255
    paint[i, 200-1] = 255
    paint[i, 200+1] = 255

for i in range(150, 250):
    paint[250, i] = 255
    paint[250+1, i] = 255
    paint[250-1, i] = 255
cv2.imshow('paint', paint)

imgDst = cv2.inpaint(img, paint, 3, cv2.INPAINT_TELEA)

cv2.imshow('Repaired image', imgDst)
cv2.waitKey(0)
