# 100 -> 200 x
# 100 -> 300 y

import cv2

img = cv2.imread('../Assets/Images/lena.jpg', 1)
dst = img[100:200, 100:300]
cv2.imshow('dst', dst)
cv2.waitKey(0)