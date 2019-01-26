#1. Load video

import cv2
cap = cv2.VideoCapture("../../Assets/Videos/panda.mp4")
isOpened = cap.isOpened()
print(isOpened)
fps = cap.get(cv2.CAP_PROP_FPS)
width = int(cap.get(cv2.CAP_PROP_FRAME_WIDTH))
height = int(cap.get(cv2.CAP_PROP_FRAME_HEIGHT))
print('FPS: ', fps, '\t Width:', width, ' \t Height:',  height)
i = 0

while(isOpened):
    if i == 50:
        break;
    else:
        i = i +1
    (flag, frame) = cap.read()
    fileName = '../../Assets/Output/'+ 'image' + str(i)+'.jpg'
    print(fileName)
    if flag == True:
        cv2.imwrite(fileName, frame, [cv2.IMWRITE_JPEG_QUALITY, 100])
print('end!')