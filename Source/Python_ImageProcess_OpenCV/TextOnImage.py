import cv2

img = cv2.imread('../../Assets/Images/flower0.jpeg', 1)
font = cv2.FONT_HERSHEY_SIMPLEX
cv2.rectangle(img, (50, 20),(250, 170), (0, 255, 0), 3)
cv2.putText(img, 'this is flower', (10, 100), font, 1, (200, 100, 255), 2, cv2.LINE_AA)

cv2.imshow('src',img)
cv2.waitKey(0)
