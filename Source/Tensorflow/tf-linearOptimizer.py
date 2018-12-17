import tensorflow as tf
import numpy as np

# Generate 100 random points
x_data=np.random.rand(100)

y_data=x_data*0.1+0.2

#Liniear Module
b=tf.Variable(0.)
k=tf.Variable(0.)
y=k*x_data+b

loss=tf.reduce_mean(tf.square(y_data-y))    #measured value-predicted value

optimizer=tf.train.GradientDescentOptimizer(0.2)

init=tf.global_variables_initializer()
train=optimizer.minimize(loss)
with tf.Session() as sess:
    sess.run(init)
    for step in range(201):
        sess.run(train)
        if step%20==0:
            print(step,sess.run([k,b]))
