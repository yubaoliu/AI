import cv2

filePath = '../../Assets/Output/'
img = cv2.imread(filePath+'image1.jpg')
imgInfo = img.shape
size = (imgInfo[1], imgInfo[0])
print(size)
fourcc = cv2.VideoWriter_fourcc(*'mp4v')
videoWrite = cv2.VideoWriter(filePath+'images2Video.mp4', fourcc, 5, size)

for i in range(1, 51):
    fileName = 'image' + str(i)+'.jpg'
    img = cv2.imread(filePath+fileName)
    videoWrite.write(img)
print('end')