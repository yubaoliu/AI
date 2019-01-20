import cv2
import numpy as np

newImageInfo = (500, 500, 3)
dst = np.zeros(newImageInfo, np.uint8)

cv2.rectangle(dst, (50, 100), (200, 300), (255, 0, 0), -1) #-1: fill; >0: line width
cv2.circle(dst, (250,250), (50), (0, 255, 0), 2)
cv2.ellipse(dst, (256, 256), (150, 100), 0, 0, 180, (255, 255, 0), -1)

points = np.array([ [150, 50], [140, 140], [200, 170], [250, 250], [150, 50]], np.int32)
print("shape", points.shape)
print(points)

points = points.reshape((-1, 1, 2))
print("inverse: ", points.shape)
print(points)

cv2.polylines(dst, [points], True, (0, 255, 255))

cv2.imshow('dst', dst)
cv2.waitKey(0)