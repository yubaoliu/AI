import cv2
import numpy as np

newImageInfo = (500, 500, 3)
dst = np.zeros(newImageInfo, np.uint8)

#line segment
cv2.line(dst, (100, 100), (400, 400), (0, 0, 255))

# line width
cv2.line(dst, (100, 200), (400, 200), (0, 255, 255), 20)
# line style
cv2.line(dst, (100, 300), (400, 300), (0, 255, 0), 20, cv2.LINE_AA)
# triangle
cv2.line(dst, (200, 150), (50, 250), (25, 100, 255))
cv2.line(dst, (50, 250), (400, 380), (25, 100, 255))
cv2.line(dst, (400, 380), (200, 150), (25, 100, 255))

cv2.imshow('dst', dst)
cv2.waitKey(0)

