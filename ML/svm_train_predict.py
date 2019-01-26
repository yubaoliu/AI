import cv2
import numpy as np
import matplotlib.pyplot as plt

# prepare the data, [height, weight]
rand1 = np.array([[155, 48], [159, 50], [164, 53], [168, 56], [172, 60]])
rand2 = np.array([[152, 53], [156, 55], [160, 56], [172, 64], [174, 65]])

plt.figure('SVM Training and Prediction')
plt.scatter(rand1[:, 0], rand1[:, 1], color='blue', marker='x', label='Negative')
plt.scatter(rand2[:, 0], rand2[:, 1], color='red', marker='o', label='Positive')

# Label, 0 -> girl; 1 -> boy, Supervised Learning
label = np.array([[0], [0], [0], [0], [0], [1], [1], [1], [1], [1]])

# data
data = np.vstack((rand1, rand2))
data = np.array(data, dtype='float32')

svm = cv2.ml.SVM_create()

svm.setType(cv2.ml.SVM_C_SVC)
svm.setKernel(cv2.ml.SVM_LINEAR)
svm.setC(0.01)

result = svm.train(data, cv2.ml.ROW_SAMPLE, label)

pt_data = np.vstack([[167, 55], [162, 57]]) #0 -> gril; 1->boy
pt_data = np.array(pt_data, dtype='float32')
print(pt_data)

(part1, part2) = svm.predict(pt_data)

print(part1, part2)
plt.show()

