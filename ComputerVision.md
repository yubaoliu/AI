---
title: "Computer Vision Study Notes"
data: 2018-12-08
---
# Prepare
## Configure Development Environment
1. Anaconda (All work should be developed in the Virtual environment)
2. OpenCV
3. Tensorflow


# Source Code
## Tensorflow
### HelloWorld.ipynb
Function: test tensorflow and OpenCV
```Python
# Test tensorflow
import tensorflow as tf  #run this statement to enable 'Tab' hint function

hello = tf.constant('hello tensor flow') #create a string
sess = tf.Session() # create a session
print(sess.run(hello))

#Test Opencv
import cv2
print('Hello Opencv')

```

# Study Resources
1. [人工智能——机器视觉及图像识别](https://www.bilibili.com/video/av33208345/?p=1)


# Reference
