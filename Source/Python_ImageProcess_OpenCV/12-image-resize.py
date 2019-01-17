
import cv2
img = cv2.imread('../../Assets/Images/lena.jpg', 1)

imgInfo = img.shape
print(imgInfo)

height = imgInfo[0]
width = imgInfo[1]
mode = imgInfo[2]

#same Scale放大 縮小 等比例
dstHeight = int(height*0.5)
dstWidth = int(width*0.5)

# 插值：最近臨域插值， 雙線性插值 像素關系重採樣 立方插值
dst = cv2.resize(img, (dstWidth, dstHeight))

cv2.imshow('image', dst)

cv2.waitKey(0)


